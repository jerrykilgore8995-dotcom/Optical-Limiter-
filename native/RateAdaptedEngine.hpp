#pragma once

#include "Engine.hpp"

#include <array>
#include <cstddef>

namespace optical {

class RateAdaptedEngine {
    static constexpr std::size_t halfDelay = 48;
    static constexpr std::size_t historySize = 128;
    static constexpr std::size_t mask = historySize - 1;
    using History = std::array<double, historySize>;

public:
    struct State {
        Engine::State circuit;
        History input, evenOutput, oddOutput;
        std::size_t index;
    };

    void save(State& state) const {
        circuit.save(state.circuit);
        state.input = input;
        state.evenOutput = evenOutput;
        state.oddOutput = oddOutput;
        state.index = index;
    }

    void restore(const State& state) {
        circuit.restore(state.circuit);
        input = state.input;
        evenOutput = state.evenOutput;
        oddOutput = state.oddOutput;
        index = state.index;
    }

    void process(float* samples, int frames, bool dry = false) {
        for (int i = 0; i < frames; ++i) {
            const double inputSample = std::isfinite(samples[i]) ? static_cast<double>(samples[i]) : 0.;
            samples[i] = static_cast<float>(step(inputSample * fullScalePeakVolts, dry) / fullScalePeakVolts);
        }
    }

    static constexpr int latencyForRate(double hostRate) {
        return hostRate < 88200.0 ? 2 * halfDelay : 0;
    }

    bool prepare(double hostRate, const Settings& settings) {
        input.fill(0.);
        evenOutput.fill(0.);
        oddOutput.fill(0.);
        index = 0;
        oversample = latencyForRate(hostRate) != 0;
        return circuit.prepare(hostRate * (oversample ? 2. : 1.), settings);
    }

    void configure(const Settings& settings) { circuit.configure(settings); }
    double reductionDb() const { return circuit.reductionDb(); }
    unsigned failedSteps() const { return circuit.failedSteps(); }

    double step(double sample, bool dry = false) {
        if (!oversample)
            return dry ? sample : circuit.step(sample);

        index = (index + 1) & mask;
        input[index] = sample;

        // The even phase of the half-band interpolator is a pure delay.
        const double even = input[(index - halfDelay) & mask];
        const double odd = 2. * oddPhase(input, index);
        evenOutput[index] = dry ? even : circuit.step(even);
        oddOutput[index] = dry ? odd : circuit.step(odd);

        // Take the even phase of the decimator. Its odd samples start one
        // host frame earlier than the newest even sample.
        if (dry) return input[(index - 2 * halfDelay) & mask];
        return .5 * evenOutput[(index - halfDelay) & mask] + oddPhase(oddOutput, index - 1);
    }

private:
    static double oddPhase(const History& history, std::size_t newest) {
        double result = 0.;
        for (std::size_t i = 0; i < taps.size(); ++i) {
            const std::size_t j = 2 * taps.size() - 1 - i;
            result += taps[i] * (history[(newest - i) & mask] + history[(newest - j) & mask]);
        }
        return result;
    }

    // 193-tap symmetric half-band FIR from a Kaiser-windowed sinc
    // with beta 10.5. The even phase is zero except for the center (0.5).
    // Listed are the first 48 odd taps, normalized to unity DC gain.
    static constexpr std::array<double, halfDelay> taps{
        -1.2271660068992396e-06, 2.6954361726209467e-06, -5.004589300925929e-06, 8.4486249282110871e-06,
        -1.3386663180585624e-05, 2.0250980682805254e-05, -2.9555167459297196e-05, 4.1902308407037097e-05,
        -5.7993097310532623e-05, 7.8633802120652686e-05, -0.00010474401845213469, 0.00013736417521929622,
        -0.00017766279346478073, 0.00022694354846009628, -0.00028665224819907645, 0.00035838392120073081,
        -0.00044389030670407789, 0.00054508816577797002, -0.00066406898935552836, 0.00080311087822555716,
        -0.00096469362403124798, 0.0011515183486282882, -0.0013665334896622242, 0.00161296949377737,
        -0.0018943853558627485, 0.0022147312135930917, -0.0025784327084517841, 0.0029905049688408831,
        -0.0034567071890417821, 0.0039837533984735333, -0.0045796020006122857, 0.0052538574467587132,
        -0.0060183344532344276, 0.0068878627922010294, -0.0078814567274757656, 0.0090240523692067271,
        -0.01034915742740594, 0.011903020111808328, -0.013751434693888625, 0.015991352876434082,
        -0.01877178376085609, 0.022333982429077033, -0.027095472662905774, 0.033846085953352277,
        -0.044280044198002522, 0.062805001801177179, -0.10558718139349087, 0.31813788967983248
    };

    Engine circuit;
    History input{}, evenOutput{}, oddOutput{};
    std::size_t index = 0;
    bool oversample = false;
};

} // namespace optical
