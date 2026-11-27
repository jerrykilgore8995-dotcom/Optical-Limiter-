#pragma once

#include "FullDeck.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace optical {

inline constexpr double fullScalePeakVolts = 17.370675219607563; // +24 dBu for a full-scale sine

class Engine {
public:
    bool prepare(double sampleRate, const Settings& controls) {
        inverse = 2 * sampleRate;
        m = MNA(unknownCount, fixedCount);
        std::copy(dcBias.begin(), dcBias.end(), m.v.begin());
        previous = older = trialStart = refinementStart = m.v;
        settings = {};
        lastInput = lastDelta = 0;
        failures = refined = 0;
        stampDeck(m, settings, 0, 0);
        if (!m.matrix.prepare() || !newton(0, 0, 50))
            return false;
        m.acceptCompanions(true);
        previous = older = m.v;
        settings = controls;
        if (!newton(0, 0, 70))
            return false;
        m.acceptCompanions(true);
        previous = older = m.v;
        m.matrix.resize(unknownCount);
        stampDeck(m, settings, 0, inverse);
        savedCompanions = m.companions;
        return m.matrix.prepare();
    }

    void configure(const Settings& next) {
        if (settings.gdb != next.gdb || settings.ff != next.ff || settings.ratio != next.ratio
                || settings.soft != next.soft || settings.sens != next.sens || settings.hpf != next.hpf
                || settings.att != next.att || settings.rel != next.rel || settings.makeup != next.makeup) {
            settings = next;
            lastDelta = 0;
        }
    }

    double step(double input) {
        const double delta = input - lastInput;
        const double predictor = std::abs(lastDelta) > 1e-9
            ? std::clamp(delta / lastDelta, -1.0, 1.25) : 0.0;
        for (int i = 0; i < unknownCount; ++i)
            m.v[i] = previous[i] + predictor * (previous[i] - older[i]);
        if (newton(input, inverse, 10)) {
            m.acceptCompanions(false);
        } else {
            const auto savedInput = lastInput;
            savedCompanions = m.companions;
            m.v = previous;
            refinementStart = previous;
            // Refine only the interval that could not converge. Keep the
            // accepted substeps and their trapezoidal histories, and roll
            // the entire host sample back if refinement runs out of budget.
            constexpr int divisions = 1024;
            int position = 0;
            int stride = divisions / 2;
            constexpr int maxAttempts = 120;
            for (int attempt = 0; position < divisions && attempt < maxAttempts && stride > 0; ++attempt) {
                stride = std::min(stride, divisions - position);
                const double fraction = static_cast<double>(position + stride) / divisions;
                const double coefficient = inverse * divisions / stride;
                if (newton(savedInput + delta * fraction, coefficient, 14)) {
                    m.acceptCompanions(false);
                    refinementStart = m.v;
                    position += stride;
                    stride = std::min(stride * 2, divisions / 2);
                } else {
                    m.v = refinementStart;
                    stride /= 2;
                }
            }
            if (position != divisions) {
                m.v = previous;
                m.companions = savedCompanions;
                ++failures;
                return voltage(outputNode);
            }
            ++refined;
        }
        older.swap(previous);
        previous = m.v;
        lastDelta = delta;
        lastInput = input;
        return voltage(outputNode);
    }

    double voltage(int node) const { return m.voltage(node); }
    double reductionDb() const {
        return 20 * std::log10((1 + 10000 * cg(voltage(cellNode))) / (1 + 10000 * 2e-8));
    }
    unsigned failedSteps() const { return failures; }
    double residual() const { return lastError; }
    int worstEquation() const { return worst; }
    unsigned refinedSteps() const { return refined; }
    int solverOperations() const { return m.matrix.operationCount(); }

private:
    void stamp(double input, double coefficient) {
        m.matrix.clear();
        projectBehavior(m, settings);
        stampDeck(m, settings, input, coefficient);
    }

    double merit() {
        double error = 0;
        double peak = 0;
        for (int i = 0; i < unknownCount; ++i) {
            const double value = m.matrix.rhs[i] * (i < nodeCount ? 1.0 : 0.00015);
            if (!std::isfinite(value))
                return std::numeric_limits<double>::infinity();
            error += value * value;
            if (std::abs(value) > peak) {
                peak = std::abs(value);
                worst = i;
            }
        }
        return error;
    }

    bool newton(double input, double coefficient, int limit) {
        stamp(input, coefficient);
        lastError = merit();
        for (int iter = 0; iter < limit; ++iter) {
            if (lastError < 1e-18) return true;
            if (!m.matrix.solve()) {
                stamp(input, coefficient);
                if (!m.matrix.solvePivoted()) return false;
            }
            trialStart = m.v;
            const double initial = lastError;
            for (double scale = 1; scale >= 1. / 8192; scale *= .5) {
                for (int i = 0; i < unknownCount; ++i)
                    m.v[i] = trialStart[i] + scale * m.matrix.x[i];
                stamp(input, coefficient);
                lastError = merit();
                if (lastError < initial) break;
            }
            if (lastError >= initial) break;
        }
        return lastError < 1e-15;
    }

    MNA m{unknownCount, fixedCount};
    Settings settings;
    std::vector<double> previous, older, trialStart, refinementStart;
    std::vector<MNA::Companion> savedCompanions;
    double inverse = 0, lastInput = 0, lastDelta = 0;
    double lastError = 0;
    unsigned failures = 0, refined = 0;
    int worst = 0;
};

} // namespace optical
