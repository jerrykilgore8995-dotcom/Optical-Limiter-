#pragma once

#include "FullDeck.hpp"
#include "ControlRipple.hpp"

#include <array>
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
    template<class V> struct PeakState {
        StateVariable<V> half, rect, position, velocity, hold;
    };
    struct State {
        History hp, hf, half, rect, avg, avgo, sub, over, ratio, timing;
        PeakState<double> peak;
        ControlRipple ripple;
        std::array<double, 7> hfPast{};
        double hfRate = 0.;
        std::size_t hfValid = 0;
        double level = 1e-9;
    };

    void configure(const Settings& next) { settings = next; tunedRate = -1.; }

    void begin(double inverse) {
        using namespace detector_model;
        if (tunedRate != inverse) {
            tunedRate = inverse;
            const auto gain = settings.ff > .5 ? 1. - 1. / settings.ratio : settings.ratio - 1.;
            ratioGain = gain;
            hp.tune(inverse, 6.283185307 * std::max(10., settings.hpf));
            hf.tune(inverse, 1. / (rhf * cfb));
            half.tune(inverse, 1. / (rhalf * chb));
            rect.tune(inverse, 1. / (rrect * crect));
            avg.tune(inverse, 1. / (ra * ca));
            avgo.tune(inverse, 1. / (raf * caf));
            sub.tune(inverse, 1. / (rsub * csub));
            over.tune(inverse, 1. / (rif * cif));
            ratio.tune(inverse, 1. / (10000. * gain * crrf));
            release.tune(inverse, 1. / settings.rel);
            attack.tune(inverse, 1. / settings.att);
            const double top = std::max(.01, 50000. * (1. - settings.sens / 30.));
            const double bottom = std::max(.01, 50000. * settings.sens / 30.);
            threshold = 3. * bottom / (top + bottom);
            peak.tune(inverse * peakSteps);
        }
        hp.begin(state.hp, inverse); hf.begin(state.hf, inverse);
        half.begin(state.half, inverse); rect.begin(state.rect, inverse);
        avg.begin(state.avg, inverse); avgo.begin(state.avgo, inverse);
        sub.begin(state.sub, inverse); over.begin(state.over, inverse);
        ratio.begin(state.ratio, inverse); attack.begin(state.timing, inverse);
        release.begin(state.timing, inverse);
        peak.begin(state, inverse);
        if (settings.sens <= 0.) darkTarget = targetFor(timeFor(commandFor(ratio(0.)))).value;
    }

    Result control(double pre, double gr) const {
        if (settings.sens <= 0.) return {darkTarget, 0., 0.};
        const auto sample = evaluate<false>(DetectorVoltage{selected(pre, gr), 1.});
        const auto target = targetFor(state.ripple.output(sample.timing, tunedRate));
        return {target.value, settings.ff * target.slope, (1. - settings.ff) * target.slope};
    }

    void accept(double pre, double gr) {
        const auto sample = evaluate<true>(selected(pre, gr));
        if (state.hfRate != tunedRate) state.hfValid = 0;
        state.hfRate = tunedRate;
        std::copy_backward(state.hfPast.begin(), state.hfPast.end() - 1, state.hfPast.end());
        state.hfPast[0]=sample.hf;
        state.hfValid = std::min(state.hfPast.size(), state.hfValid + 1);
        save(state.hp, sample.hp); save(state.hf, sample.hf);
        save(state.half, sample.half); save(state.rect, sample.rect);
        save(state.avg, sample.avg); save(state.avgo, sample.avgo);
        state.peak = sample.peak;
        save(state.sub, sample.sub);
        save(state.over, sample.over); save(state.ratio, sample.ratio);
        save(state.timing, sample.timing);
        state.ripple.accept(sample.hf, sample.timing, tunedRate);
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
        V hp, hf, half, rect, avg, avgo;
        V sub, over, ratio, timing, level;
        PeakState<V> peak;
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
        r.half = half(max(r.hf, V{0.}));
        r.rect = rect(r.half * 2. - r.hf);
        if (AdvanceAll || settings.ratio < 8.) {
            r.avg = avg(r.rect);
            r.avgo = avgo(r.avg * eavg1);
        }
        if (AdvanceAll || settings.ratio >= 8.) {
            r.peak = peak.advance(state.peak, r.hf);
        }
        r.level = max(settings.ratio < 8. ? r.avgo : r.peak.hold.value, V{1e-9});
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

    // Resolve the short charge pulse on smaller intervals while keeping the
    // circuit's audio rate, interpolator and trapezoidal companions intact.
    static constexpr int peakSteps = 4;
    class PeakIntegrator {
    public:
        void tune(double inverse) {
            using namespace detector_model;
            rate = inverse;
            pole(1. / (rhalf * chb), halfGain, halfHistory);
            pole(1. / (rrect * crect), rectGain, rectHistory);
            constexpr double discharge = 1. / (rp * cp);
            constexpr double charging = 1. / (8. * cp);
            constexpr double driveLimit = 1. / (100. * cp);
            freeGain = 1. / (inverse + discharge);
            chargeHistory = 1. / (inverse + discharge + charging);
            chargeGain = charging * chargeHistory;
            limitHistory = 1. / (inverse + discharge + driveLimit);
            limitBias = 13.5 * driveLimit * limitHistory;
            if (inverse > 0.) {
                positionCoupling = (gpkint1 / cpkfb) / inverse;
                velocityHistory = 1. / (inverse + gpkdamp1 / cpkvel + (gpkint1 / cpkvel) * positionCoupling);
                velocityGain = gpkint1 / cpkvel * velocityHistory;
                invRate = 1. / inverse;
            }
        }

        void begin(const State& state, double mainRate) {
            // Causal fractional-delay FIR at quarter-sample positions. Fit through
            // 0.47 pi rad/sample with unity DC gain and the intended DC delay.
            static constexpr double weights[peakSteps - 1][8] = {
                {0.087044829961725939, 1.2177424344667742, -0.56808518142735165, 0.47088371753083397, -0.3344693796132695, 0.18051674528844114, -0.066442592896427968, 0.012817917320156981},
                {0.2614342092785068, 1.2232409791429912, -0.95424666841387562, 0.8567513280241239, -0.62970078738126078, 0.34668231057488869, -0.1293934744134787, 0.025249894381531784},
                {0.55552474059080603, 0.86924117290026393, -0.87470710236641158, 0.83535386375161547, -0.63217698297898051, 0.3543177097555506, -0.13396380258539459, 0.026430258051734802}
            };
            if (state.hfRate == mainRate && state.hfValid == state.hfPast.size()) {
                for (int i = 0; i < peakSteps - 1; ++i) {
                    double offset = 0.;
                    for (std::size_t j = 0; j < state.hfPast.size(); ++j)
                        offset += weights[i][j + 1] * state.hfPast[j];
                    interpolation[i] = {offset, weights[i][0]};
                }
                return;
            }
            // A step-size change invalidates the uniform FIR history.
            for (int i = 0; i < peakSteps - 1; ++i) {
                const double f = static_cast<double>(i + 1) / peakSteps;
                const double h01 = f * f * (3. - 2. * f);
                const double h11 = f * f * (f - 1.);
                const double tangent = mainRate == 0. ? 0. : state.hf.slope * 2. / mainRate;
                const double gain = h01 + 2. * h11;
                interpolation[i] = {state.hf.value * (1. - gain) + tangent * f * (1. - f), gain};
            }
        }

        template<class V>
        PeakState<V> advance(const PeakState<double>& history, V input) const {
            PeakState<V> state{copy<V>(history.half), copy<V>(history.rect), copy<V>(history.position),
                               copy<V>(history.velocity), copy<V>(history.hold)};
            if (rate > 0.) {
                for (const auto& entry : interpolation)
                    tick(state, input * entry.slope + entry.value);
            }
            tick(state, input);
            return state;
        }

    private:
        void pole(double g, double& input, double& history) const {
            history = 1. / (rate + g);
            input = g * history;
        }
        template<class V> static StateVariable<V> copy(const History& h) { return {V{h.value}, V{h.slope}}; }
        template<class V> V past(const StateVariable<V>& state) const {
            return rate == 0. ? V{0.} : state.value * rate + state.slope;
        }
        template<class V> void update(StateVariable<V>& s, V next) const {
            s.slope = rate == 0. ? V{0.} : (next - s.value) * rate - s.slope;
            s.value = next;
        }
        template<class V> void tick(PeakState<V>& s, V hf) const {
            using namespace detector_model;
            update(s.half, max(hf, V{0.}) * halfGain + past(s.half) * halfHistory);
            update(s.rect, (s.half.value * 2. - hf) * rectGain + past(s.rect) * rectHistory);
            if (rate > 0.) {
                const auto previousPosition = past(s.position) * invRate;
                update(s.velocity, (hf - previousPosition) * velocityGain + past(s.velocity) * velocityHistory);
                update(s.position, previousPosition + s.velocity.value * positionCoupling);
            } else {
                update(s.position, hf);
                update(s.velocity, V{0.});
            }
            const auto edge = hf * epkedge1 + s.position.value * epkedge2 + s.velocity.value * epkedge3;
            const auto target = max(bcrest(s.rect.value, hf, s.position.value), bpkdrive(bpkshift(edge, hf)));
            const auto bias = past(s.hold);
            update(s.hold, max(bias * freeGain,
                               min(bias * limitHistory + limitBias, target * chargeGain + bias * chargeHistory)));
        }

        std::array<History, peakSteps - 1> interpolation{};
        double rate = 0., invRate = 0.;
        double halfGain = 0., halfHistory = 0., rectGain = 0., rectHistory = 0.;
        double positionCoupling = 0., velocityHistory = 0., velocityGain = 0.;
        double freeGain = 0., chargeHistory = 0., chargeGain = 0., limitHistory = 0., limitBias = 0.;
    };

    Settings settings;
    Pole hp, hf, half, rect, avg, avgo, sub, over, ratio, attack, release;
    PeakIntegrator peak;
    double tunedRate = -1., threshold = 0., ratioGain = 0., darkTarget = 2e-8;
};

} // namespace optical
