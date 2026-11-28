#pragma once

#include "Detector.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace optical {

inline constexpr double fullScalePeakVolts = 17.370675219607563; // +24 dBu for a full-scale sine

class Engine {
public:
    // Only accepted history is needed to continue processing on a
    // separately prepared solver. Matrix scratch stays thread-local.
    struct State {
        std::array<double, unknownCount + fixedCount> solution, previous, older;
        std::array<MNA::Companion, companionCount> companions;
        Detector detector;
        Settings settings;
        double lastInput, lastDelta;
        unsigned failures, refined;
    };

    void save(State& state) const {
        std::copy(m.v.begin(), m.v.end(), state.solution.begin());
        std::copy(previous.begin(), previous.end(), state.previous.begin());
        std::copy(older.begin(), older.end(), state.older.begin());
        std::copy(m.companions.begin(), m.companions.end(), state.companions.begin());
        state.detector = detector;
        state.settings = settings;
        state.lastInput = lastInput;
        state.lastDelta = lastDelta;
        state.failures = failures;
        state.refined = refined;
    }

    void restore(const State& state) {
        std::copy(state.solution.begin(), state.solution.end(), m.v.begin());
        std::copy(state.previous.begin(), state.previous.end(), previous.begin());
        std::copy(state.older.begin(), state.older.end(), older.begin());
        std::copy(state.companions.begin(), state.companions.end(), m.companions.begin());
        detector = state.detector;
        settings = state.settings;
        lastInput = state.lastInput;
        lastDelta = state.lastDelta;
        failures = state.failures;
        refined = state.refined;
    }

    bool prepare(double sampleRate, const Settings& controls) {
        inverse = 2 * sampleRate;
        m = MNA(unknownCount, fixedCount);
        m.companions.resize(companionCount);
        std::copy(dcBias.begin(), dcBias.end(), m.v.begin());
        previous = older = trialStart = refinementStart = m.v;
        settings = {};
        detector = {};
        detector.configure(settings);
        lastInput = lastDelta = 0;
        failures = refined = 0;
        fixedStamps = LinearStamp(unknownCount, unknownCount + fixedCount);
        stampFixed(fixedStamps);
        fixedStamps.build();
        if (!prepareMatrix(0) || !newton(0, 0, 50))
            return false;
        accept(true);
        previous = older = m.v;
        settings = controls;
        if (settings.ff <= .5) settings.soft = 0;
        detector.configure(settings);
        if (!newton(0, 0, 70))
            return false;
        accept(true);
        previous = older = m.v;
        if (!prepareMatrix(inverse)) return false;
        subdivisionMatrix = std::move(m.matrix);
        if (!prepareMatrix<true>(inverse)) return false;
        savedCompanions = m.companions;
        return true;
    }

    void configure(Settings next) {
        if (next.ff <= .5) next.soft = 0;
        if (settings.gdb != next.gdb || settings.ff != next.ff || settings.ratio != next.ratio
                || settings.soft != next.soft || settings.sens != next.sens || settings.hpf != next.hpf
                || settings.att != next.att || settings.rel != next.rel || settings.makeup != next.makeup) {
            settings = next;
            detector.configure(settings);
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
            accept(false);
        } else {
            savedCompanions = m.companions;
            savedDetector = detector.state;
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

    double voltage(int node) const {
        if (node == timingNode) return detector.state.timing.value;
        if (node == detectorNode) return detector.state.level;
        return m.voltage(node);
    }
    double reductionDb() const {
        return 20 * std::log10((1 + 10000 * cg(voltage(cellNode))) / (1 + 10000 * 2e-8));
    }
    unsigned failedSteps() const { return failures; }
    double residual() const { return lastError; }
    int worstEquation() const {
        int row = 0;
        double peak = 0;
        for (int i = 0; i < unknownCount; ++i) {
            const double value = std::abs(m.matrix.rhs[i]) * (i < nodeCount ? 1. : .00015);
            if (value > peak) {
                peak = value;
                row = i;
            }
        }
        return row;
    }
    unsigned refinedSteps() const { return refined; }
    int solverOperations() const { return m.matrix.operationCount(); }

private:
    void accept(bool dc) {
        m.acceptCompanions(dc);
        detector.accept(m.voltage(preNode), m.voltage(grNode));
    }

    void restoreSample() {
        m.v = previous;
        m.companions = savedCompanions;
        detector.state = savedDetector;
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
                accept(false);
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
                accept(false);
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

    template<bool FixedCaps = false>
    bool prepareMatrix(double coefficient) {
        m.matrix.resize(unknownCount);
        setKnownTerminals(m, settings, 0);
        detector.begin(coefficient);
        const auto control = detector.control(m.voltage(preNode), m.voltage(grNode));
        projectBehavior(m, settings, control.target);
        fixedStamps.addJacobian(m.matrix);
        if constexpr (FixedCaps) stampCapacitorJacobian(m, coefficient);
        m.matrix.beginDynamicStamp();
        stampDynamic<FixedCaps>(m, settings, coefficient, control.target, control.preSlope, control.grSlope);
        if (!m.matrix.prepare()) return false;
        m.matrix.clear();
        fixedStamps.addJacobian(m.matrix);
        if constexpr (FixedCaps) stampCapacitorJacobian(m, coefficient);
        m.matrix.saveLinearEntries();
        m.matrix.cachePreparedFactorization();
        return true;
    }

    template<bool FixedCaps>
    void stamp(double input, double coefficient, bool raw = false) {
        setKnownTerminals(m, settings, input);
        const auto control = detector.control(m.voltage(preNode), m.voltage(grNode));
        projectBehavior(m, settings, control.target);
        fixedStamps.apply(m.matrix, m.v.data(), raw);
        stampDynamic<FixedCaps>(m, settings, coefficient, control.target, control.preSlope, control.grSlope);
    }

    double merit() {
        double error = 0;
        for (int i = 0; i < nodeCount; ++i) {
            const double value = m.matrix.rhs[i];
            error += value * value;
        }
        for (int i = nodeCount; i < unknownCount; ++i) {
            const double value = m.matrix.rhs[i] * .00015;
            error += value * value;
        }
        return std::isfinite(error) ? error : std::numeric_limits<double>::infinity();
    }

    bool newton(double input, double coefficient, int limit) {
        if (coefficient == inverse) return solveStep<true>(input, coefficient, limit);
        if (coefficient == 0.) return solveStep<false>(input, coefficient, limit);
        std::swap(m.matrix, subdivisionMatrix);
        const bool solved = solveStep<false>(input, coefficient, limit);
        std::swap(m.matrix, subdivisionMatrix);
        return solved;
    }

    template<bool FixedCaps>
    bool solveStep(double input, double coefficient, int limit) {
        detector.begin(coefficient);
        stamp<FixedCaps>(input, coefficient);
        lastError = merit();
        for (int iter = 0; iter < limit; ++iter) {
            if (lastError < 1e-18) return true;
            if (!m.matrix.solve()) {
                stamp<FixedCaps>(input, coefficient, true);
                if (!m.matrix.solvePivoted()) return false;
            }
            trialStart = m.v;
            const double initial = lastError;
            for (double scale = 1; scale >= 1. / 8192; scale *= .5) {
                for (int i = 0; i < unknownCount; ++i)
                    m.v[i] = trialStart[i] + scale * m.matrix.x[i];
                // A single nodal residual can already exceed the full
                // merit limit. Keep the last trial exact for the fallback.
                if (coefficient > 0 && scale > 1. / 8192
                        && rejectTrial(m, fixedStamps, coefficient, initial * (1. + 1e-10))) continue;
                stamp<FixedCaps>(input, coefficient);
                lastError = merit();
                if (lastError < initial) break;
            }
            if (lastError >= initial) break;
        }
        return lastError < 1e-15;
    }

    MNA m{unknownCount, fixedCount};
    Sparse subdivisionMatrix;
    Settings settings;
    Detector detector;
    Detector::State savedDetector;
    LinearStamp fixedStamps{0, 0};
    std::vector<double> previous, older, trialStart, refinementStart;
    std::vector<MNA::Companion> savedCompanions;
    double inverse = 0, lastInput = 0, lastDelta = 0;
    double lastError = 0;
    unsigned failures = 0, refined = 0;
};

} // namespace optical
