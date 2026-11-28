#include "RateAdaptedEngine.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string_view>

namespace {

bool checkDriveThreshold() {
    constexpr double rate = 44100.;
    struct Case {
        double inputDbu, trimDb, sensitivity, expectedGr;
        double ratio = 2., frequency = 1000.;
    };
    constexpr Case cases[]{
        {4., -18., 1., 2.98}, {4., 0., 1., 2.98}, {4., 18., 1., 2.98},
        {-24., 18., 29., 2.98}, {24., -18., -17.6, 2.96},
        {4., 18., optical::thresholdOffSensitivity, 0.},
        {4., 0., 1., 5.14, 8., 5000.}, {4., 18., 1., 5.14, 8., 5000.},
    };
    for (const auto& test : cases) {
        optical::Settings settings;
        settings.ff = 1.;
        settings.soft = 0.;
        settings.ratio = test.ratio;
        settings.gdb = test.trimDb;
        settings.sens = test.sensitivity;
        optical::RateAdaptedEngine engine;
        if (!engine.prepare(rate, settings)) return false;
        const double peak = std::numbers::sqrt2 * .775 * std::pow(10., test.inputDbu / 20.);
        for (int i = 0; i < static_cast<int>(rate * .75); ++i)
            engine.step(peak * std::sin(2. * std::numbers::pi * test.frequency * i / rate));
        const double reduction = engine.reductionDb();
        if (engine.failedSteps() || std::abs(reduction - test.expectedGr) > .15) {
            std::cerr << "Threshold check at " << test.inputDbu << " dBu and " << test.trimDb
                      << " dB trim returned " << reduction << " dB reduction, expected "
                      << test.expectedGr << '\n';
            return false;
        }
    }
    return true;
}

bool checkControlRipple() {
    constexpr double rate = 88200.;
    constexpr double inverse = 2. * rate;
    for (const double frequency : {20., 40., 50., 100., 1000.}) {
        optical::ControlRipple ripple;
        double total = 0., energy = 0.;
        int count = 0;
        for (int i = 0; i < static_cast<int>(rate * .6); ++i) {
            const double phase = 2. * std::numbers::pi * frequency * i / rate;
            double timing = .8;
            for (int harmonic = 2; harmonic <= 18; harmonic += 2)
                timing += .02 * std::cos(harmonic * phase) / harmonic;
            const double output = ripple.output(timing, inverse);
            if (i >= static_cast<int>(rate * .45)) {
                const double offset = output - .8;
                total += offset;
                energy += offset * offset;
                ++count;
            }
            ripple.accept(std::sin(phase), timing, inverse);
        }
        const double offset = total / count;
        if (std::abs(offset) > .01 || energy / count - offset * offset > 1e-10) return false;
    }
    return true;
}

bool checkCellFit() {
    constexpr double rate = 88200.;
    constexpr int count = 8820;
    constexpr double peak = .5 * std::numbers::sqrt2;
    for (const double resistance : {3000., 30000.}) {
        optical::CellDistortion model;
        double k2 = 0., k4 = 0.;
        double previous = 0.;
        for (int i = 0; i < count * 4; ++i) {
            const double phase = 2. * std::numbers::pi * 1000. * i / rate;
            const double volts = peak * std::cos(phase);
            if (i >= count * 3) {
                const double correction = model.fit.evenAc(volts);
                k2 += correction * std::cos(2. * phase);
                k4 += correction * std::cos(4. * phase);
            }
            model.accept(volts, 1. / resistance, 2. * rate, previous < 0. && volts >= 0., 1000.);
            previous = volts;
        }
        const double total = 2. * std::hypot(k2, k4) / (count * peak);
        const double target = resistance == 3000. ? .0000495 : .0001705;
        if (std::abs(total / target - 1.) > .01 || std::abs(k4 / k2 - .2) > .002) return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--drive-threshold") return checkDriveThreshold() ? 0 : 1;
    if (argc == 2 && std::string_view(argv[1]) == "--ripple") return checkControlRipple() ? 0 : 1;
    if (argc == 2 && std::string_view(argv[1]) == "--cell-fit") return checkCellFit() ? 0 : 1;
    if (argc < 8 || argc > 10 || (argc == 10 && std::string_view(argv[9]) != "burst")) {
        std::cerr << "Usage: OpticalEngineTest FF soft ratio sensitivity input_dBu seconds sample_rate [frequency] [burst]\n";
        return 2;
    }
    optical::Settings s;
    s.ff = std::atof(argv[1]);
    s.soft = std::atof(argv[2]);
    s.ratio = std::atof(argv[3]);
    s.sens = std::atof(argv[4]);
    if (s.sens == 0.) s.sens = optical::thresholdOffSensitivity;
    const double amplitude = std::sqrt(2.) * .775 * std::pow(10., std::atof(argv[5]) / 20.);
    const double seconds = std::atof(argv[6]);
    const double rate = std::atof(argv[7]);
    const double frequency = argc >= 9 ? std::atof(argv[8]) : 1000;
    const bool burst = argc == 10;
    optical::Engine engine;
    if (!engine.prepare(rate, s)) {
        std::cerr << "Operating-point initialization failed at " << optical::unknownNames[engine.worstEquation()]
                  << ", squared residual " << engine.residual() << '\n';
        return 1;
    }
    std::cout << std::setprecision(13) << "time,pre,gr,out,cell,det,timing\n";
    for (int i = 0; i < static_cast<int>(seconds * rate); ++i) {
        const double t = i / rate;
        const double level = burst && (t < .025 || t >= .205) ? 0. : amplitude;
        engine.step(level * std::sin(2 * std::numbers::pi * frequency * t));
        std::cout << t << ',' << engine.voltage(optical::preNode) << ',' << engine.voltage(optical::grNode)
                  << ',' << engine.voltage(optical::outputNode) << ',' << engine.voltage(optical::cellNode)
                  << ',' << engine.voltage(optical::detectorNode) << ',' << engine.voltage(optical::timingNode) << '\n';
    }
    std::cerr << "Missed steps: " << engine.failedSteps() << " of " << static_cast<int>(seconds * rate) << '\n';
    return engine.failedSteps() != 0;
}
