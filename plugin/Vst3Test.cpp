#include <juce_audio_processors/juce_audio_processors.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <numbers>

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
        double outEnergy = 0, inEnergy = 0;
        for (int block = 0; block < 10; ++block) {
            for (int i = 0; i < 512; ++i) {
                const auto input = static_cast<float>(.10 * std::sin(2 * std::numbers::pi * 1000
                                                                    * (block * 512 + i) / 44100.));
                audio.setSample(0, i, input);
                if (channels == 2) audio.setSample(1, i, 0);
                if (block == 9) inEnergy += input * input;
            }
            instance->processBlock(audio, midi);
            for (int i = 0; i < 512; ++i) {
                const auto output = audio.getSample(0, i);
                if (!std::isfinite(output)) return 1;
                if (channels == 2 && std::abs(audio.getSample(1, i)) > .001f) return 1;
                if (block == 9) outEnergy += output * output;
            }
        }
        const double ratio = std::sqrt(outEnergy / inEnergy);
        std::cout << channels << " channel(s), output/input RMS " << ratio << '\n';
        if (!(ratio > .001 && ratio < .95)) return 1;
        if (channels == 1) instance->releaseResources();
    }
    {
        const std::unique_ptr<juce::AudioProcessorEditor> editor(instance->createEditorIfNeeded());
        if (!editor || editor->getWidth() != 1200 || editor->getHeight() != 247) return 1;
        const auto image = editor->createComponentSnapshot(editor->getLocalBounds());
        auto output = juce::File(juce::String::fromUTF8(argv[2])).createOutputStream();
        if (!output || !juce::PNGImageFormat().writeImageToStream(image, *output)) return 1;
    }
    instance->releaseResources();
    return 0;
}
