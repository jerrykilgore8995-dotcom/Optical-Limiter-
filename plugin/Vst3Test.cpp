#include <juce_audio_processors/juce_audio_processors.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <memory>
#include <numbers>
#include <vector>

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter();

namespace {

bool checkHighFrequencyResponse(juce::AudioPluginInstance& instance) {
    auto& parameters = instance.getParameters();
    const auto threshold = std::find_if(parameters.begin(), parameters.end(), [](const auto* p) {
        return p->getName(128) == "Threshold Sensitivity";
    });
    if (threshold == parameters.end()) return false;
    (*threshold)->setValueNotifyingHost(0.f);

    constexpr double rate = 44100.;
    constexpr int blockSize = 512;
    constexpr int settle = 4410;
    constexpr int measure = 4410;
    constexpr double amplitude = .007943282347; // -42 dBFS peak, equivalent to -18 dBu.
    std::array<double, 2> gain{};
    int band = 0;
    for (const double frequency : {1000., 20000.}) {
        instance.prepareToPlay(rate, blockSize);
        if (instance.getLatencySamples() != 96) return false;
        juce::AudioBuffer<float> audio(2, blockSize);
        juce::MidiBuffer midi;
        std::complex<double> output{};
        std::vector<double> blockTimes;
        for (int frame = 0; frame < settle + measure; frame += blockSize) {
            for (int i = 0; i < blockSize; ++i) {
                const double angle = 2 * std::numbers::pi * frequency * (frame + i) / rate;
                audio.setSample(0, i, static_cast<float>(amplitude * std::sin(angle)));
                audio.setSample(1, i, 0.f);
            }
            const auto start = std::chrono::steady_clock::now();
            instance.processBlock(audio, midi);
            if (frame >= settle)
                blockTimes.push_back(std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - start).count());
            for (int i = 0; i < blockSize && frame + i < settle + measure; ++i) {
                if (frame + i < settle) continue;
                const double angle = 2 * std::numbers::pi * frequency * (frame + i) / rate;
                output += double(audio.getSample(0, i)) * std::complex<double>{std::cos(angle), -std::sin(angle)};
            }
        }
        gain[band++] = 20 * std::log10(2 * std::abs(output) / (measure * amplitude));
        std::sort(blockTimes.begin(), blockTimes.end());
        std::cout << "Threshold-off " << frequency << " Hz, one driven channel, median "
                  << blockTimes[blockTimes.size() / 2] << " ms/512, max "
                  << blockTimes.back() << " ms/512\n";
        instance.releaseResources();
    }
    const double rolloff = gain[1] - gain[0];
    std::cout << "Packaged VST3 through gain at 1 kHz " << gain[0] << " dB, at 20 kHz "
              << gain[1] << " dB, relative response " << rolloff << " dB, latency 96 samples\n";
    return std::abs(rolloff + .06927) < .5 && std::abs(gain[0] + .15097) < .25;
}

bool checkBypass(juce::AudioProcessor& processor, bool useHostBypass, double rate) {
    constexpr int impulse = 9;
    constexpr int blockSize = 73;
    processor.prepareToPlay(rate, blockSize);
    const int latency = rate < 88200. ? 96 : 0;
    if (processor.getLatencySamples() != latency) return false;
    auto* bypass = processor.getBypassParameter();
    if (!bypass) return false;
    bypass->setValueNotifyingHost(useHostBypass ? 0.f : 1.f);
    juce::AudioBuffer<float> audio(2, blockSize);
    juce::MidiBuffer midi;
    for (int block = 0; block < 3; ++block) {
        audio.clear();
        if (block == 0) audio.setSample(0, impulse, .25f);
        if (useHostBypass) processor.processBlockBypassed(audio, midi);
        else processor.processBlock(audio, midi);
        for (int i = 0; i < blockSize; ++i) {
            const float expected = block * blockSize + i == impulse + latency ? .25f : 0.f;
            if (audio.getSample(0, i) != expected || audio.getSample(1, i) != 0.f) return false;
        }
    }
    bypass->setValueNotifyingHost(0.f);
    processor.releaseResources();
    return true;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    juce::ScopedJuceInitialiser_GUI gui;
    juce::AudioPluginFormatManager formats;
    auto format = std::make_unique<juce::VST3PluginFormat>();
    juce::OwnedArray<juce::PluginDescription> descriptions;
    format->findAllTypesForFile(descriptions, juce::String::fromUTF8(argv[1]));
    formats.addFormat(format.release());
    if (descriptions.size() != 1) {
        std::cerr << "VST3 scan did not return one plug-in\n";
        return 1;
    }
    juce::String error;
    auto instance = formats.createPluginInstance(*descriptions[0], 44100, 512, error);
    if (!instance) {
        std::cerr << error << '\n';
        return 1;
    }
    for (int channels : {1, 2}) {
        const auto bus = channels == 1 ? juce::AudioChannelSet::mono() : juce::AudioChannelSet::stereo();
        juce::AudioProcessor::BusesLayout layout;
        layout.inputBuses.add(bus);
        layout.outputBuses.add(bus);
        if (!instance->setBusesLayout(layout)) return 1;
        instance->prepareToPlay(44100, 512);
        juce::AudioBuffer<float> audio(channels, 512);
        juce::MidiBuffer midi;
        std::array<double, 2> outEnergy{}, inEnergy{};
        std::array<double, 10> times{};
        constexpr int warmup = 2;
        for (int block = 0; block < warmup + static_cast<int>(times.size()); ++block) {
            const bool rightDriven = channels == 2 && block >= warmup;
            for (int i = 0; i < 512; ++i) {
                const auto input = static_cast<float>(.10 * std::sin(2 * std::numbers::pi * 1000
                                                                    * (block * 512 + i) / 44100.));
                audio.setSample(0, i, input);
                if (channels == 2) {
                    const auto right = rightDriven ? static_cast<float>(.08 * std::sin(
                        2 * std::numbers::pi * 733 * (block * 512 + i) / 44100.)) : 0.f;
                    audio.setSample(1, i, right);
                    if (block == warmup + static_cast<int>(times.size()) - 1) inEnergy[1] += right * right;
                }
                if (block == warmup + static_cast<int>(times.size()) - 1) inEnergy[0] += input * input;
            }
            const auto start = std::chrono::steady_clock::now();
            instance->processBlock(audio, midi);
            if (block >= warmup)
                times[static_cast<std::size_t>(block - warmup)] =
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - start).count();
            for (int i = 0; i < 512; ++i) {
                for (int channel = 0; channel < channels; ++channel) {
                    const auto output = audio.getSample(channel, i);
                    if (!std::isfinite(output) || (channel == 1 && !rightDriven && std::abs(output) > .001f))
                        return 1;
                    if (block == warmup + static_cast<int>(times.size()) - 1) outEnergy[channel] += output * output;
                }
            }
        }
        const double ratio = std::sqrt(outEnergy / inEnergy);
        std::sort(times.begin(), times.end());
        std::cout << channels << " channel(s), output/input RMS " << ratio
                  << ", median " << (times[4] + times[5]) / 2
                  << " ms/512, max " << times.back() << " ms/512" << '\n';
        if (!(ratio > .001 && ratio < .95)) return 1;
        if (channels == 1) instance->releaseResources();
    }
    instance->releaseResources();
    if (!checkHighFrequencyResponse(*instance)) return 1;
    {
        // A platform host wraps the plug-in editor in a native child window.
        // Capture the same editor directly so it can paint headlessly.
        const std::unique_ptr<juce::AudioProcessor> panelOwner(createPluginFilter());
        if (!checkBypass(*panelOwner, false, 44100)
                || !checkBypass(*panelOwner, true, 44100)
                || !checkBypass(*panelOwner, false, 88200)) return 1;
        std::cout << "Bypass impulse alignment verified at 44.1 and 88.2 kHz\n";
        panelOwner->prepareToPlay(44100, 512);
        juce::AudioBuffer<float> stereo(2, 512);
        juce::MidiBuffer midi;
        for (int block = 0; block < 7; ++block) {
            for (int i = 0; i < 512; ++i) {
                const auto sine = static_cast<float>(std::sin(2 * std::numbers::pi * 1000
                                                           * (block * 512 + i) / 44100.));
                stereo.setSample(0, i, .045f * sine);
                stereo.setSample(1, i, .023f * sine);
            }
            panelOwner->processBlock(stereo, midi);
        }
        const std::unique_ptr<juce::AudioProcessorEditor> editor(panelOwner->createEditorIfNeeded());
        if (!editor || editor->getWidth() != 1200 || editor->getHeight() != 247) return 1;
        const auto image = editor->createComponentSnapshot(editor->getLocalBounds());
        auto output = juce::File(juce::String::fromUTF8(argv[2])).createOutputStream();
        if (!output || !output->setPosition(0) || output->truncate().failed()
                || !juce::PNGImageFormat().writeImageToStream(image, *output)) return 1;
    }
    instance->releaseResources();
    return 0;
}
