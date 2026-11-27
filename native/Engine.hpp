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
            savedCompanions = m.companions;
            // Uniform substeps keep the feedback peak detector's history
            // consistent across a burst-to-silence transition.
            const bool solved = settings.ff > .5
                ? refineAdaptive(delta) || refineUniform(delta)
                : refineUniform(delta) || refineAdaptive(delta);
            if (!solved) {
                restoreSample();
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
    void restoreSample() {
        m.v = previous;
        m.companions = savedCompanions;
    }

    bool refineUniform(double delta) {
        for (int steps : {2, 4, 8, 16, 32}) {
            restoreSample();
            bool complete = true;
            for (int i = 1; i <= steps; ++i) {
                const double fraction = static_cast<double>(i) / steps;
                if (!newton(lastInput + delta * fraction, inverse * steps, 16)) {
                    complete = false;
                    break;
                }
                m.acceptCompanions(false);
            }
            if (complete) return true;
        }
        return false;
    }

    bool refineAdaptive(double delta) {
        restoreSample();
        refinementStart = previous;
        // Keep each accepted trapezoidal history and refine only the
        // interval that could not converge.
        constexpr int divisions = 1024;
        constexpr int maxAttempts = 120;
        int position = 0;
        int stride = divisions / 2;
        for (int attempt = 0; position < divisions && attempt < maxAttempts && stride > 0; ++attempt) {
            stride = std::min(stride, divisions - position);
            const double fraction = static_cast<double>(position + stride) / divisions;
            const double coefficient = inverse * divisions / stride;
            if (newton(lastInput + delta * fraction, coefficient, 14)) {
                m.acceptCompanions(false);
                refinementStart = m.v;
                position += stride;
                stride = std::min(stride * 2, divisions / 2);
            } else {
                m.v = refinementStart;
                stride /= 2;
            }
        }
        return position == divisions;
    }

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
