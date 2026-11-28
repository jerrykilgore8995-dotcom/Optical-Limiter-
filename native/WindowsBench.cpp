#include "RateAdaptedEngine.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <numbers>

int main() {
    constexpr double rate = 44100;
    constexpr int blockSize = 512;
    constexpr int blocks = 2;
    for (double frequency : {0., 1000., 20000.}) {
        auto settings = optical::Settings{};
        if (frequency == 20000.) settings.sens = 0;
        const double amplitude = optical::fullScalePeakVolts * (frequency == 20000. ? .007943282347 : .10);
        auto sample = [=](int i) { return amplitude * std::sin(2 * std::numbers::pi * frequency * i / rate); };
        optical::RateAdaptedEngine engine;
        if (!engine.prepare(rate, settings)) return 1;
        for (int i = 0; i < blockSize; ++i) engine.step(sample(i));
        optical::benchStamps = optical::benchSolves = optical::benchPivoted = 0;
        const auto start = std::chrono::steady_clock::now();
        double energy = 0;
        for (int i = blockSize; i < blockSize * (blocks+1); ++i) {
            const auto output = engine.step(sample(i));
            energy += output * output;
        }
        const auto elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        std::cout << frequency << " Hz " << elapsed / blocks << " ms/512"
                  << ", stamps " << optical::benchStamps / blocks
                  << ", solves " << optical::benchSolves / blocks << ", pivoted " << optical::benchPivoted / blocks
                  << ", missed " << engine.failedSteps() << ", RMS " << std::sqrt(energy / (blocks * blockSize)) << std::endl;
    }
}
