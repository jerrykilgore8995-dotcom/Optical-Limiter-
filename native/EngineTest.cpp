#include "Engine.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <string_view>

namespace {

bool checkControlRipple() {
    constexpr double rate = 88200.;
    constexpr double inverse = 2. * rate;
    for (const double frequency : {100., 1000.}) {
        optical::ControlRipple ripple;
        double error = 0.;
        int count = 0;
        for (int i = 0; i < static_cast<int>(rate * .25); ++i) {
            const double phase = 2. * std::numbers::pi * frequency * i / rate;
            const double timing = .8 + .01 * std::cos(2. * phase)
                + .005 * std::cos(4. * phase) + .0025 * std::cos(6. * phase);
            const double output = ripple.output(timing, inverse);
            if (frequency == 1000. && output != timing) return false;
            if (i >= static_cast<int>(rate * .15)) {
                error += (output - .8) * (output - .8);
                ++count;
            }
            ripple.accept(std::sin(phase), timing, inverse);
        }
        if (frequency == 100. && std::sqrt(error / count) > .00001) return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string_view(argv[1]) == "--ripple") return checkControlRipple() ? 0 : 1;
    if (argc < 8 || argc > 10 || (argc == 10 && std::string_view(argv[9]) != "burst")) {
        std::cerr << "Usage: OpticalEngineTest FF soft ratio sensitivity input_dBu seconds sample_rate [frequency] [burst]\n";
        return 2;
    }
    optical::Settings s;
    s.ff = std::atof(argv[1]);
    s.soft = std::atof(argv[2]);
    s.ratio = std::atof(argv[3]);
    s.sens = std::atof(argv[4]);
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
