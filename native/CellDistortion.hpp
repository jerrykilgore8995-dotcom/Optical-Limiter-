#pragma once

#include <algorithm>
#include <cmath>

namespace optical {

struct CellCurrentFit {
    double peak = 0.;
    double second = 0.;
    double fourth = 0.;

    double evenAc(double voltage) const {
        if (peak <= 0.) return 0.;
        const double x = std::clamp(voltage / peak, -1.3, 1.3);
        const double square = x * x;
        return peak * (second * (2. * square - 1.)
             + fourth * (8. * square * (square - 1.) + 1.));
    }
};

// Propagate the fitted nonlinear current through the linearized feedback
// and output networks. The main audio solve retains its original odd
// distortion, including at levels where an even fit would mix into K3.
struct EvenResponse {
    double previous = 0., feedback = 0., dc = 0., output = 0.;

    double step(const CellCurrentFit& fit, double volts, double conductance,
                double inverseRate, double makeupGain) {
        const double parallel = conductance + .0001;
        constexpr double cellEvenGain = 2.;
        const double drive = -cellEvenGain * (conductance / parallel) * fit.evenAc(volts);
        const double rc = inverseRate * 100e-12 / parallel;
        feedback += (drive + previous - 2. * feedback) / (1. + rc);
        previous = drive;
        constexpr double outputTime = 47e-6 / (1. / 20000. + 1. / 10047.);
        dc += (feedback - dc) / (1. + .5 * inverseRate * outputTime);
        output = makeupGain * (feedback - dc) * (10000. / 10047.);
        return output;
    }
};

// Fit the cell current at a given AC terminal voltage and illumination.
class CellDistortion {
public:
    void accept(double volts, double conductance, double inverse, bool crossing, double frequency) {
        if (inverse <= 0.) {
            *this = {};
            return;
        }
        const double dt = 2. / inverse;
        const double square = volts * volts;
        energy += dt * (previousSquare + square) * .5;
        previousSquare = square;
        elapsed += dt;
        if ((!crossing || frequency <= 0.) && elapsed < .075) return;

        const double rms = std::sqrt(energy / elapsed);
        const double blend = std::min(1., elapsed / (elapsed + .02));
        voltageRms = ready ? voltageRms + blend * (rms - voltageRms) : rms;
        ready = true;
        energy = elapsed = 0.;
        if (voltageRms < 1e-5) {
            fit = {};
            return;
        }

        const double g = std::clamp(conductance, 2e-8, .005);
        const double logResistance = std::clamp(std::log10(1. / (3000. * g)), 0., 1.);
        const double x = std::min(voltageRms, 10.);
        // The breakpoints and slopes come from the two log-log traces in
        // PerkinElmer Figure 3(a). These values are amplitude ratios.
        const double at3k = std::max(2e-5, .000198 * x * x);
        const double at30k = std::max(4e-5, .000682 * x * x);
        double total = at3k * std::exp(logResistance * std::log(at30k / at3k));
        // The requested 10 kHz K2 reference applies on the deeply lit
        // side outside the catalog's resistance range.
        const double highTone = std::max(0., 1. - std::abs(frequency - 10000.) / 1000.);
        const double wellLit = std::clamp((3000. * g - 1.) / .35, 0., 1.);
        const double feedbackWeight = g / (g + .0001);
        total += highTone * wellLit * (std::max(total, .0001 / feedbackWeight) - total);
        total *= std::min(1., voltageRms / .01);
        const double fourthRatio = frequency > 5000. ? .1 * std::clamp((5500. - frequency) / 500., 0., 1.)
            : .2 - .1 * std::clamp((frequency - 1000.) / 4000., 0., 1.);
        fit.peak = std::sqrt(2.) * voltageRms;
        fit.second = total / std::sqrt(1. + fourthRatio * fourthRatio);
        fit.fourth = fit.second * fourthRatio;
    }

    CellCurrentFit fit;

private:
    double previousSquare = 0., energy = 0., elapsed = 0., voltageRms = 0.;
    bool ready = false;
};

} // namespace optical
