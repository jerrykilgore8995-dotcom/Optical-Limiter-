#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace optical {

// Suppress periodic CV ripple from the low-frequency detector. Predictable
// gain modulation at twice the audio frequency generates odd audio harmonics.
class ControlRipple {
public:
    template<class V>
    V output(V timing, double inverse) const {
        if (wet == 0. || inverse != configuredRate) return timing;
        V filtered = timing;
        for (const auto& notch : notches) filtered = notch.output(filtered);
        return timing + (filtered - timing) * wet;
    }

    void accept(double audio, double timing, double inverse) {
        if (inverse <= 0.) {
            *this = {};
            return;
        }
        const double dt = 2. / inverse;
        elapsed += dt;
        if (audio >= 0. && previousAudio < 0.) {
            const double fraction = audio / (audio - previousAudio);
            const double sinceCrossing = dt * fraction;
            const double measuredPeriod = elapsed - sinceCrossing;
            elapsed = sinceCrossing;
            const bool inBand = measuredPeriod >= 1. / 300. && measuredPeriod <= 1. / 45.;
            const bool stable = inBand && std::abs(measuredPeriod - period) < .025 * measuredPeriod;
            stableCycles = stable ? std::min(stableCycles + 1, 3) : 0;
            period = inBand ? measuredPeriod : 0.;
        }
        previousAudio = audio;
        const bool locked = stableCycles >= 2 && elapsed < 1.5 * period;
        if (!locked && (elapsed >= 1.5 * period || inverse != configuredRate)) wet = 0.;
        if (locked && (inverse != configuredRate || std::abs(configuredPeriod - period) > .0005 * period)) {
            const bool resetHistory = wet == 0. || inverse != configuredRate;
            for (std::size_t i = 0; i < notches.size(); ++i) {
                const double angle = 4. * std::numbers::pi * (i + 1) * dt / period;
                notches[i].tune(angle, timing, resetHistory);
            }
            configuredRate = inverse;
            configuredPeriod = period;
        }
        if (inverse != configuredRate) return;
        double filtered = timing;
        for (auto& notch : notches) filtered = notch.accept(filtered);
        // Leave timing transients on the original path when the filtered
        // envelope departs materially from the existing control voltage.
        const double wanted = locked
            ? std::clamp((.075 - std::abs(filtered - timing)) / .04, 0., 1.) : 0.;
        wet += std::min(1., dt * 350.) * (wanted - wet);
    }

private:
    struct Notch {
        double b0 = 1., b1 = 0., a2 = 0.;
        double x1 = 0., x2 = 0., y1 = 0., y2 = 0.;
        void tune(double angle, double timing, bool reset) {
            constexpr double quality = 1.4;
            const double width = std::sin(angle) / (2. * quality);
            b0 = 1. / (1. + width);
            b1 = -2. * std::cos(angle) * b0;
            a2 = (1. - width) * b0;
            if (reset) x1 = x2 = y1 = y2 = timing;
        }
        template<class V>
        V output(V x) const { return x * b0 + b0 * x2 + b1 * (x1 - y1) - a2 * y2; }
        double accept(double x) {
            const double y = output(x);
            x2 = x1; x1 = x;
            y2 = y1; y1 = y;
            return y;
        }
    };

    std::array<Notch, 3> notches{};
    double previousAudio = 0., elapsed = 0., period = 0.;
    double configuredRate = 0., configuredPeriod = 0., wet = 0.;
    int stableCycles = 0;
};

} // namespace optical
