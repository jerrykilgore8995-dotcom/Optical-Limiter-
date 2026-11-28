#pragma once
#include <algorithm>
#include <cmath>
#include <iterator>

namespace optical {

class ControlRipple {
public:
    template<class V> V output(V timing, double inverse) const {
        if (inverse <= 0. || wet == 0.) return timing;
        return timing + (V{mean} - timing) * wet;
    }

    double frequency() const { return period > 0. ? 1. / period : 0.; }

    bool accept(double audio, double timing, double inverse, double loopRatio = 1., bool peak = false, double attack = .01) {
        if (!(inverse > 0.)) {
            *this = {};
            return false;
        }
        const double dt = 2. / inverse;
        const double slope = timing - previousTiming;
        elapsed += dt;
        // Recover quickly through the original detector on peaks that exceed
        // the recent cycle and when a low-frequency tone stops altogether.
        const double absAudio = std::abs(audio);
        if (absAudio >= .13 * referencePeak) belowSignal = 0.;
        else belowSignal += dt;
        const double prevCyclePeak = cyclePeak;
        cyclePeak = std::max(cyclePeak, absAudio);
        if (referencePeak > 1e-5 && (absAudio > 1.18 * referencePeak
                || belowSignal > .13 * period)) {
            stableCycles = 0;
            stableTime = 0.;
        }

        const bool crossed = audio >= 0. && previousAudio < 0.;
        if (crossed) {
            const double fraction = -previousAudio / (audio - previousAudio);
            const double completed = elapsed - dt * (1. - fraction);
            const double crossingTiming = previousTiming + fraction * slope;
            integral += dt * fraction * (previousTiming + crossingTiming) * .5;
            const double current = integral / completed;
            const bool inBand = completed >= 1. / 22000. && completed <= 1. / 15.;
            const bool stable = inBand && std::abs(completed - period) < .025 * completed
                && std::abs(prevCyclePeak - referencePeak) < .12 * prevCyclePeak;
            if (stable && prevCyclePeak > 1e-5) {
                stableCycles = std::min(stableCycles + 1, 4);
                stableTime += completed;
            } else {
                stableCycles = 0;
                stableTime = 0.;
            }
            if (stableCycles == 2 || wet == 0.)
                mean = current;
            else if (stableCycles > 2) {
                const double correction = settledCorrection(1. / completed, loopRatio > 1., peak);
                const double active = std::clamp(mean, 0., 1.);
                const double near10k = std::max(0., 1. - std::abs(1. / completed - 10000.) / 1500.);
                const double quiet = std::clamp((.22 - mean) / .14, 0., 1.)
                                   * std::clamp(mean / .035, 0., 1.);
                const double nearKnee = peak ? near10k * quiet * (loopRatio > 1. ? -.358 : -.375) : 0.;
                mean += (current - mean) * (.65 / loopRatio) + .065 * (correction * active + nearKnee);
            }
            period = inBand ? completed : 0.;
            referencePeak = prevCyclePeak;
            elapsed = dt * (1. - fraction);
            cyclePeak = absAudio;
            integral = elapsed * (crossingTiming + timing) * .5;
        } else {
            integral += dt * (previousTiming + timing) * .5;
        }
        previousAudio = audio;
        previousTiming = timing;
        const bool locked = stableCycles >= 2 && stableTime >= 5. * attack
                         && elapsed < 1.3 * period;
        wet += std::min(1., dt * (locked ? 400. : 1200.)) * (double(locked) - wet);
        if (!locked && wet < 1e-6) wet = 0.;
        return crossed;
    }

private:
    // Effective steady attenuation of the reference includes a small
    // contribution from its periodic detector modulation. Match its
    // fundamental gain without reproducing that modulation.
    static double settledCorrection(double frequency, bool feedback, bool peak) {
        constexpr double f[]{15., 20., 40., 50., 100., 200., 8000., 10000., 12000.};
        constexpr double ffAvg[]{.015, .016, .026, .027, .016, 0., 0., 0., 0.};
        constexpr double ffPeak[]{.045, .046, .033, .028, .004, 0., 0., -.024, 0.};
        constexpr double fbAvg[]{.15, .135, .067, .049, .035, 0., 0., 0., 0.};
        constexpr double fbPeak[]{.04, .055, .150, .115, .06, 0., 0., -.095, 0.};
        const auto& map = feedback ? (peak ? fbPeak : fbAvg) : (peak ? ffPeak : ffAvg);
        const auto* high = std::lower_bound(std::begin(f), std::end(f), frequency);
        if (high == std::begin(f)) return map[0];
        if (high == std::end(f)) return 0.;
        const auto index = high - std::begin(f);
        const double blend = (frequency - high[-1]) / (*high - high[-1]);
        return map[index - 1] + blend * (map[index] - map[index - 1]);
    }

    double previousAudio = 0., previousTiming = 0., elapsed = 0., integral = 0.;
    double cyclePeak = 0., referencePeak = 0., belowSignal = 0., period = 0.;
    double mean = 0., wet = 0., stableTime = 0.;
    int stableCycles = 0;
};

} // namespace optical
