#include "Engine.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numbers>

int main(int argc, char** argv) {
    if (argc != 8 && argc != 9) {
        std::cerr << "Usage: OpticalEngineTest FF soft ratio sensitivity input_dBu seconds sample_rate [frequency]\n";
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
    const double frequency = argc == 9 ? std::atof(argv[8]) : 1000;
    optical::Engine engine;
    if (!engine.prepare(rate, s)) {
        std::cerr << "Operating-point initialization failed at " << optical::unknownNames[engine.worstEquation()]
                  << ", squared residual " << engine.residual() << '\n';
        return 1;
    }
    std::cout << std::setprecision(13) << "time,pre,gr,out,cell,det,timing\n";
    for (int i = 0; i < static_cast<int>(seconds * rate); ++i) {
        const double t = i / rate;
        engine.step(amplitude * std::sin(2 * std::numbers::pi * frequency * t));
        std::cout << t << ',' << engine.voltage(optical::preNode) << ',' << engine.voltage(optical::grNode)
                  << ',' << engine.voltage(optical::outputNode) << ',' << engine.voltage(optical::cellNode)
                  << ',' << engine.voltage(optical::detectorNode) << ',' << engine.voltage(optical::timingNode) << '\n';
    }
    std::cerr << "Missed steps: " << engine.failedSteps() << " of " << static_cast<int>(seconds * rate) << '\n';
    return engine.failedSteps() != 0;
}
