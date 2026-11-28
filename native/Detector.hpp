#pragma once

#include "FullDeck.hpp"
#include "ControlRipple.hpp"
#include "PeakCalibration.hpp"

#include <type_traits>

namespace optical {

// A scalar with the derivative with respect to the selected sidechain input.
// Piecewise expressions follow the same branch as the SPICE model.
struct DetectorVoltage {
    double value = 0., slope = 0.;
    DetectorVoltage() = default;
    DetectorVoltage(double x, double dx = 0.) : value(x), slope(dx) {}
};
inline DetectorVoltage operator+(DetectorVoltage a, DetectorVoltage b) {
    return {a.value + b.value, a.slope + b.slope};
}
inline DetectorVoltage operator-(DetectorVoltage a, DetectorVoltage b) {
    return {a.value - b.value, a.slope - b.slope};
}
inline DetectorVoltage operator*(DetectorVoltage a, double b) {
    return {a.value * b, a.slope * b};
}
inline DetectorVoltage operator*(double a, DetectorVoltage b) { return b * a; }
inline DetectorVoltage operator*(DetectorVoltage a, DetectorVoltage b) {
    return {a.value * b.value, a.slope * b.value + a.value * b.slope};
}
inline DetectorVoltage operator/(DetectorVoltage a, double b) { return a * (1. / b); }
inline DetectorVoltage min(DetectorVoltage a, DetectorVoltage b) { return a.value < b.value ? a : b; }
inline DetectorVoltage max(DetectorVoltage a, DetectorVoltage b) { return a.value > b.value ? a : b; }
inline DetectorVoltage abs(DetectorVoltage a) { return a.value < 0 ? a * -1. : a; }

class Detector {
public:
    struct Result { double target, preSlope, grSlope; };
    template<class V> struct StateVariable { V value{}, slope{}; };
    using History = StateVariable<double>;
    struct State {
        History hp, hf, avg, sub, over, ratio, timing;
        double peak = 0., peakGain = 1.;
        double quietTime = 0., onsetTime = 0.;
        ControlRipple ripple;
        CellDistortion cellDistortion;
        double level = 1e-9;
    };

    void configure(const Settings& next) { settings = next; tunedRate = -1.; }

    void begin(double inverse) {
        using namespace detector_model;
        if (tunedRate != inverse) {
            tunedRate = inverse;
            const auto gain = settings.ff > .5 ? 1. - 1. / settings.ratio : settings.ratio - 1.;
            ratioGain = gain;
            constexpr double discharge = 1. / (rp * cp);
            peakDecay = inverse > 0. ? (inverse - discharge) / (inverse + discharge) : 0.;
            hp.tune(inverse, 6.283185307 * std::max(10., settings.hpf));
            hf.tune(inverse, 1. / (rhf * cfb));
            avg.tune(inverse, 1. / (ra * ca));
            sub.tune(inverse, 1. / (rsub * csub));
            over.tune(inverse, 1. / (rif * cif));
            ratio.tune(inverse, 1. / (10000. * gain * crrf));
            release.tune(inverse, 1. / settings.rel);
            attack.tune(inverse, 1. / settings.att);
            const double top = std::max(.01, 50000. * (1. - settings.sens / 30.));
            const double bottom = std::max(.01, 50000. * settings.sens / 30.);
            threshold = 3. * bottom / (top + bottom);
        }
        hp.begin(state.hp, inverse); hf.begin(state.hf, inverse);
        avg.begin(state.avg, inverse);
        sub.begin(state.sub, inverse); over.begin(state.over, inverse);
        ratio.begin(state.ratio, inverse); attack.begin(state.timing, inverse);
        release.begin(state.timing, inverse);
        if (settings.sens <= 0.) darkTarget = targetFor(timeFor(commandFor(ratio(0.)))).value;
    }

    Result control(double pre, double gr) const {
        if (settings.sens <= 0.) return {darkTarget, 0., 0.};
        const auto sample = evaluate<false>(DetectorVoltage{selected(pre, gr), 1.});
        const auto target = targetFor(state.ripple.output(sample.timing, tunedRate));
        return {target.value, settings.ff * target.slope, (1. - settings.ff) * target.slope};
    }

    void accept(double pre, double gr, double cellVolts = 0., double cellConductance = 2e-8) {
        const auto sample = evaluate<true>(selected(pre, gr));
        save(state.hp, sample.hp); save(state.hf, sample.hf);
        if (tunedRate > 0.) {
            const double dt = 2. / tunedRate;
            if (std::abs(sample.hf) < .0001) {
                state.quietTime = std::min(.01, state.quietTime + dt);
            } else {
                if (state.quietTime >= .002) state.onsetTime = .00012;
                state.quietTime = 0.;
            }
            state.onsetTime = std::max(0., state.onsetTime - dt);
        }
        save(state.avg, sample.avg);
        state.peak = sample.peak;
        save(state.sub, sample.sub);
        save(state.over, sample.over); save(state.ratio, sample.ratio);
        save(state.timing, sample.timing);
        const bool crossing = state.ripple.accept(sample.hf, sample.timing, tunedRate,
                                                  settings.ff > .5 ? 1. : settings.ratio,
                                                  settings.ratio >= 8., settings.att);
        const double frequency = state.ripple.frequency();
        if (crossing && settings.ratio >= 8.) state.peakGain = peakCalibration(frequency, state.peak);
        state.cellDistortion.accept(cellVolts, cellConductance, tunedRate, crossing, frequency);
        state.level = sample.level;
    }

    State state;

private:
    struct Pole {
        double gain = 0., invDenominator = 0., offset = 0.;
        void tune(double inverse, double rate) {
            invDenominator = 1. / (inverse + rate);
            gain = rate * invDenominator;
        }
        void begin(const History& h, double inverse) {
            offset = inverse == 0 ? 0. : (inverse * h.value + h.slope) * invDenominator;
        }
        template<class V> V operator()(V input) const { return input * gain + offset; }
    };

    template<class V> struct Sample {
        V hp, hf, avg, avgo;
        V sub, over, ratio, timing, level;
        V peak;
    };
    static double value(double v) { return v; }
    static double value(DetectorVoltage v) { return v.value; }
    static double log10(double v) { return std::log10(v); }
    static DetectorVoltage log10(DetectorVoltage v) {
        return {std::log10(v.value), v.slope / (2.302585092994046 * v.value)};
    }
    double selected(double pre, double gr) const { return settings.ff * pre + (1. - settings.ff) * gr; }

    template<class V> V commandFor(V rrout) const {
        return min(V{4.}, max(V{0.}, rrout * -1.));
    }

    template<class V> V timeFor(V command) const {
        const auto normal = release(command);
        return value(command) > value(normal) ? attack(command) : normal;
    }

    static DetectorVoltage targetFor(DetectorVoltage value) {
        const auto timing = min(4., max(0., value));
        constexpr double scale = 1.151292546497023;
        if (timing.value <= .000173700424)
            return {2e-8, 0.};
        const double e = std::exp(scale * timing.value);
        return min(.005, max(2e-8, DetectorVoltage{(e - 1.) / 10000., scale * e * timing.slope / 10000.}));
    }

    template<bool AdvanceAll, class V> Sample<V> evaluate(V input) const {
        using namespace detector_model;
        Sample<V> r{};
        r.hp = hp(input);
        r.hf = hf((settings.hpf > 0 ? input - r.hp : input) * -1.);
        const auto rectified = max(r.hf, r.hf * -1.);
        if (AdvanceAll || settings.ratio < 8.) {
            r.avg = avg(rectified);
            r.avgo = r.avg * eavg1;
        }
        if (AdvanceAll || settings.ratio >= 8.) {
            // A sharp edge after silence charges the reference peak
            // detector above its steady envelope. Match that initial
            // charge with a short onset gain on the ideal rectifier.
            const double onsetGain = settings.ff > .5 && (state.quietTime >= .002
                                                || state.onsetTime > 0.) ? 1.58 : 1.;
            r.peak = max(rectified * onsetGain, V{state.peak * peakDecay});
        }
        r.level = max(settings.ratio < 8. ? r.avgo : r.peak * state.peakGain, V{1e-9});
        r.sub = sub(log10(max(V{.001}, min(V{15.}, r.level))) * -2. - threshold);
        r.over = over(r.sub * -1.);
        const auto above = max(V{-15.}, min(V{15.}, r.over)) * 10.;
        auto knee = max(V{0.}, above);
        if (settings.sens <= 0.) knee = V{0.};
        else if (settings.soft > .5 && value(above) > -5. && value(above) < 5.) {
            const auto x = above + 5.;
            if constexpr (std::is_same_v<V, double>) knee = x * x / 20.;
            else knee = {x.value * x.value / 20., x.value * x.slope / 10.};
        }
        r.ratio = ratio(knee * (-.1 * ratioGain));
        r.timing = timeFor(commandFor(r.ratio));
        return r;
    }

    void save(History& h, double value) const {
        h.slope = tunedRate == 0. ? 0. : tunedRate * (value - h.value) - h.slope;
        h.value = value;
    }

    Settings settings;
    Pole hp, hf, avg, sub, over, ratio, attack, release;
    double peakDecay = 0., tunedRate = -1., threshold = 0., ratioGain = 0., darkTarget = 2e-8;
};

} // namespace optical
