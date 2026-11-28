#include <juce_audio_processors/juce_audio_processors.h>
#include "RateAdaptedEngine.hpp"
#include "BinaryData.h"

#include <array>
#include <atomic>
#include <cmath>
#include <memory>

namespace {

using juce::Colour;
using juce::Graphics;
using juce::Rectangle;
using juce::String;
using State = juce::AudioProcessorValueTreeState;

constexpr int panelWidth = 1200;
constexpr int panelHeight = 247;
constexpr std::array ratios{2, 3, 4, 8, 12, 20};
const Colour white{0xffffffff};
const Colour ink{0xff191322};
const Colour glossInk{0xff0d0816};
const Colour black{0xff101016};
const Colour menuSurface{0xffead8fa};
const Colour needle{0xff801da0};
constexpr float rotaryStart = 3.9269908f;
constexpr float rotaryEnd = 8.6393798f;

enum Control { drive, sensitivity, ratio, attack, release, makeup, hpf, topology, knee, bypass, controlCount };
constexpr std::array ids{"input", "threshold", "ratio", "attack", "release", "output",
                         "sidechain_hpf", "topology", "knee", "bypass"};

class Processor final : public juce::AudioProcessor {
public:
    Processor() : AudioProcessor(BusesProperties().withInput("Input", juce::AudioChannelSet::stereo(), true)
                                                 .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
                  state(*this, nullptr, "OpticalLimiter", makeParameters()) {
        for (std::size_t i = 0; i < ids.size(); ++i)
            raw[i] = state.getRawParameterValue(ids[i]);
    }

    const String getName() const override { return "Optical Limiter"; }
    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }
    bool isMidiEffect() const override { return false; }
    double getTailLengthSeconds() const override { return 3.0; }
    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram(int) override {}
    const String getProgramName(int) override { return {}; }
    void changeProgramName(int, const String&) override {}
    bool hasEditor() const override { return true; }
    juce::AudioProcessorEditor* createEditor() override;
    juce::AudioProcessorParameter* getBypassParameter() const override { return state.getParameter(ids[bypass]); }

    bool isBusesLayoutSupported(const BusesLayout& layout) const override {
        const auto channels = layout.getMainInputChannelSet();
        return channels == layout.getMainOutputChannelSet()
            && (channels == juce::AudioChannelSet::mono() || channels == juce::AudioChannelSet::stereo());
    }

    void prepareToPlay(double rate, int) override {
        ready = true;
        const auto settings = readSettings();
        for (int ch = 0; ch < getTotalNumInputChannels(); ++ch) {
            ready = engines[static_cast<std::size_t>(ch)].prepare(rate, settings) && ready;
        }
        setLatencySamples(optical::RateAdaptedEngine::latencyForRate(rate));
        status.store(ready, std::memory_order_relaxed);
        for (auto& meter : meters) meter.store(0, std::memory_order_relaxed);
    }

    void releaseResources() override {
        ready = false;
        status.store(false, std::memory_order_relaxed);
    }

    void processBlock(juce::AudioBuffer<float>& audio, juce::MidiBuffer& midi) override {
        if (read(bypass) >= .5f) {
            processBlockBypassed(audio, midi);
            return;
        }
        juce::ScopedNoDenormals noDenormals;
        if (!ready || audio.getNumSamples() == 0)
            return;
        const auto settings = readSettings();
        bool converged = true;
        for (int ch = 0; ch < audio.getNumChannels(); ++ch) {
            const auto index = static_cast<std::size_t>(ch);
            auto& engine = engines[index];
            engine.configure(settings);
            float* samples = audio.getWritePointer(ch);
            const auto before = engine.failedSteps();
            for (int i = 0; i < audio.getNumSamples(); ++i) {
                const double input = std::isfinite(samples[i]) ? static_cast<double>(samples[i]) : 0.;
                samples[i] = static_cast<float>(engine.step(input * optical::fullScalePeakVolts)
                                                 / optical::fullScalePeakVolts);
            }
            converged = engine.failedSteps() == before && converged;
            meters[index].store(static_cast<float>(engine.reductionDb()), std::memory_order_relaxed);
        }
        if (audio.getNumChannels() < 2) meters[1].store(0, std::memory_order_relaxed);
        status.store(converged, std::memory_order_relaxed);
    }

    void processBlockBypassed(juce::AudioBuffer<float>& audio, juce::MidiBuffer&) override {
        zeroMeters();
        if (!ready) return;
        for (int ch = 0; ch < audio.getNumChannels(); ++ch) {
            auto& engine = engines[static_cast<std::size_t>(ch)];
            float* samples = audio.getWritePointer(ch);
            for (int i = 0; i < audio.getNumSamples(); ++i) {
                const double input = std::isfinite(samples[i]) ? static_cast<double>(samples[i]) : 0.;
                samples[i] = static_cast<float>(engine.step(input * optical::fullScalePeakVolts, true)
                                                 / optical::fullScalePeakVolts);
            }
        }
    }

    void getStateInformation(juce::MemoryBlock& data) override {
        if (const auto xml = state.copyState().createXml()) copyXmlToBinary(*xml, data);
    }
    void setStateInformation(const void* data, int size) override {
        if (const auto xml = getXmlFromBinary(data, size); xml && xml->hasTagName(state.state.getType()))
            state.replaceState(juce::ValueTree::fromXml(*xml));
    }

    State state;
    std::array<std::atomic<float>, 2> meters{};
    std::atomic<bool> status{false};

private:
    static State::ParameterLayout makeParameters() {
        State::ParameterLayout p;
        const auto range = [&](Control control, const char* name, juce::NormalisableRange<float> values,
                               float initial, const char* unit) {
            p.add(std::make_unique<juce::AudioParameterFloat>(
                juce::ParameterID{ids[control], 1}, name, values, initial,
                juce::AudioParameterFloatAttributes().withLabel(unit)));
        };
        range(drive, "Input", {-18.f, 18.f, .1f}, 0, "dB");
        range(sensitivity, "Threshold Sensitivity", {0.f, 30.f, .1f}, 15, "dB");
        p.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID{ids[ratio], 1}, "Ratio",
              juce::StringArray{"2:1", "3:1", "4:1", "8:1", "12:1", "20:1"}, 2));
        range(attack, "Attack", {1.f, 40.f, .1f, .63f}, 10, "ms");
        range(release, "Release", {10.f, 3000.f, 1.f, .30f}, 100, "ms");
        range(makeup, "Output Makeup", {0.f, 30.f, 1.f}, 0, "dB");
        juce::StringArray highPass{"Off"};
        for (int hz = 10; hz <= 200; hz += 10)
            highPass.add(String(hz) + " Hz");
        p.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID{ids[hpf], 1}, "Sidechain HPF", highPass, 0));
        p.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID{ids[topology], 1}, "Detector Topology",
              juce::StringArray{"Feedback", "Feed-forward"}, 1));
        p.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID{ids[knee], 1}, "Feed-forward Knee",
              juce::StringArray{"Hard", "Soft"}, 1));
        p.add(std::make_unique<juce::AudioParameterBool>(juce::ParameterID{ids[bypass], 1}, "Bypass", false));
        return p;
    }

    float read(Control control) const { return raw[control]->load(std::memory_order_relaxed); }

    optical::Settings readSettings() const {
        optical::Settings s;
        s.gdb = read(drive);
        s.sens = read(sensitivity);
        s.ratio = ratios[static_cast<std::size_t>(read(ratio))];
        s.att = .001 * read(attack);
        s.rel = .001 * read(release);
        s.makeup = read(makeup);
        s.hpf = 10 * read(hpf);
        s.ff = read(topology);
        s.soft = s.ff > .5 ? read(knee) : 0.;
        return s;
    }

    void zeroMeters() {
        for (auto& meter : meters) meter.store(0, std::memory_order_relaxed);
    }

    std::array<optical::RateAdaptedEngine, 2> engines;
    std::array<const std::atomic<float>*, controlCount> raw{};
    bool ready = false;
};

void label(Graphics& g, const String& text, Rectangle<int> bounds, float size, Colour colour = ink,
           int style = juce::Font::plain) {
    g.setColour(colour);
    g.setFont(juce::Font(juce::FontOptions(size, style)));
    g.drawFittedText(text, bounds, juce::Justification::centred, 1);
}

class Hardware final : public juce::LookAndFeel_V4 {
public:
    Hardware() {
        setColour(juce::ComboBox::backgroundColourId, menuSurface);
        setColour(juce::ComboBox::textColourId, glossInk);
        setColour(juce::ComboBox::outlineColourId, Colour{0xff672480});
        setColour(juce::ComboBox::arrowColourId, ink);
        setColour(juce::PopupMenu::backgroundColourId, menuSurface);
        setColour(juce::PopupMenu::textColourId, ink);
        setColour(juce::PopupMenu::highlightedBackgroundColourId, Colour{0xffcf9eed});
        setColour(juce::PopupMenu::highlightedTextColourId, glossInk);
        setColour(juce::TextButton::buttonColourId, menuSurface);
        setColour(juce::TextButton::buttonOnColourId, Colour{0xffd8acf4});
        setColour(juce::TextButton::textColourOffId, glossInk);
        setColour(juce::TextButton::textColourOnId, glossInk);
    }

    void drawRotarySlider(Graphics& g, int x, int y, int w, int h, float proportion,
                          float start, float end, juce::Slider&) override {
        const float d = static_cast<float>(std::min(w, h));
        const auto area = Rectangle<float>{static_cast<float>(x), static_cast<float>(y),
                                           static_cast<float>(w), static_cast<float>(h)}
                              .withSizeKeepingCentre(d, d);
        const auto centre = area.getCentre();
        const bool large = d > 65;
        const float radius = d * .5f;
        g.setColour(Colour{0xff130b1d});
        g.fillEllipse(area.translated(2, 4));
        g.setGradientFill({Colour{0xff8b8997}, centre.x, area.getY(),
                           Colour{0xff100e16}, centre.x, area.getBottom(), false});
        g.fillEllipse(area);
        g.setColour(black);
        g.fillEllipse(area.reduced(1.5f));
        const auto bevel = area.reduced(3.f);
        g.setGradientFill({Colour{0xff575461}, centre.x, bevel.getY(),
                           Colour{0xff191720}, centre.x, bevel.getBottom(), false});
        g.fillEllipse(bevel);
        for (int i = 0; i < (large ? 36 : 24); ++i) {
            const float angle = juce::MathConstants<float>::twoPi * i / (large ? 36.f : 24.f);
            const auto a = centre.getPointOnCircumference(radius - 8.f, angle);
            const auto b = centre.getPointOnCircumference(radius - 3.5f, angle);
            g.setColour(Colour{0xff0b0a12}.withAlpha(.7f));
            g.drawLine({a, b}, large ? 2.1f : 1.3f);
        }
        const auto cap = area.reduced(large ? 10.f : 7.f);
        g.setGradientFill({Colour{0xff292733}, centre.x - radius * .35f, cap.getY(),
                           Colour{0xff14131c}, centre.x + radius * .45f, cap.getBottom(), false});
        g.fillEllipse(cap);
        g.setColour(Colour{0xff858291}.withAlpha(.45f));
        g.drawEllipse(cap, .75f);
        const float angle = start + proportion * (end - start);
        g.setColour(white);
        g.drawLine({centre.getPointOnCircumference(radius * .47f, angle),
                    centre.getPointOnCircumference(radius - (large ? 12.f : 8.f), angle)}, large ? 3.f : 2.f);
    }
};

struct Dial {
    Control parameter;
    const char* title;
    int x;
    bool large;
};
constexpr std::array dials{
    Dial{drive, "INPUT", 125, true}, Dial{sensitivity, "THRESHOLD", 247, false},
    Dial{ratio, "RATIO", 337, false}, Dial{attack, "ATTACK", 863, false},
    Dial{release, "RELEASE", 953, false}, Dial{makeup, "OUTPUT", 1075, true}
};

class Panel final : public juce::AudioProcessorEditor, private juce::Timer {
public:
    explicit Panel(Processor& owner) : AudioProcessorEditor(owner), processor(owner) {
        setLookAndFeel(&hardware);
        for (std::size_t i = 0; i < dials.size(); ++i) {
            auto& slider = knobs[i];
            const auto& dial = dials[i];
            slider.setName(dial.title);
            slider.setSliderStyle(juce::Slider::RotaryHorizontalVerticalDrag);
            slider.setRotaryParameters(rotaryStart, rotaryEnd, true);
            slider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
            slider.setPopupDisplayEnabled(true, true, this);
            slider.onValueChange = [this] { repaint(); };
            addAndMakeVisible(slider);
            bindings[i] = std::make_unique<State::SliderAttachment>(processor.state, ids[dial.parameter], slider);
        }

        hpfSelector.setName("Sidechain high-pass");
        configureSelector(hpfSelector);
        const auto& cutoff = static_cast<juce::AudioParameterChoice&>(*processor.state.getParameter(ids[hpf]));
        for (int i = 0; i < cutoff.choices.size(); ++i)
            hpfSelector.addItem(cutoff.choices[i], i + 1);
        hpfBinding = std::make_unique<State::ComboBoxAttachment>(processor.state, ids[hpf], hpfSelector);

        modeSelector.setName("Detector mode");
        configureSelector(modeSelector);
        modeSelector.addItemList({"Feedback", "Feed-forward Hard", "Feed-forward Soft"}, 1);
        auto& topologyParam = *processor.state.getParameter(ids[topology]);
        auto& kneeParam = *processor.state.getParameter(ids[knee]);
        const auto syncMode = [this](float) {
            const auto ff = processor.state.getRawParameterValue(ids[topology])->load();
            const auto soft = processor.state.getRawParameterValue(ids[knee])->load();
            modeSelector.setSelectedId(ff < .5f ? 1 : soft < .5f ? 2 : 3, juce::dontSendNotification);
        };
        topologyBinding = std::make_unique<juce::ParameterAttachment>(topologyParam, syncMode);
        kneeBinding = std::make_unique<juce::ParameterAttachment>(kneeParam, syncMode);
        modeSelector.onChange = [this] {
            const int mode = modeSelector.getSelectedId();
            if (mode < 1) return;
            if (mode > 1) kneeBinding->setValueAsCompleteGesture(mode == 3 ? 1.f : 0.f);
            topologyBinding->setValueAsCompleteGesture(mode == 1 ? 0.f : 1.f);
        };
        topologyBinding->sendInitialUpdate();
        bypassButton.setName("Bypass");
        addAndMakeVisible(bypassButton);
        bypassBinding = std::make_unique<State::ButtonAttachment>(processor.state, ids[bypass], bypassButton);
        bypassButton.onStateChange = [this] {
            bypassButton.setButtonText(bypassButton.getToggleState() ? "BYPASSED" : "IN");
        };
        bypassButton.setButtonText(bypassButton.getToggleState() ? "BYPASSED" : "IN");

        setSize(panelWidth, panelHeight);
        startTimerHz(30);
    }

    ~Panel() override { setLookAndFeel(nullptr); }

    void resized() override {
        for (std::size_t i = 0; i < dials.size(); ++i) {
            const auto& dial = dials[i];
            const int size = dial.large ? 93 : 47;
            knobs[i].setBounds(dial.x - size / 2, (dial.large ? 126 : 113) - size / 2, size, size);
        }
        hpfSelector.setBounds(222, 180, 114, 22);
        modeSelector.setBounds(399, 199, 197, 22);
        bypassButton.setBounds(605, 199, 197, 22);
    }

    void paint(Graphics& g) override {
        paintChassis(g);
        label(g, "OPTICAL LIMITER", {400, 14, 400, 31}, 26, glossInk, juce::Font::bold);
        for (std::size_t i = 0; i < dials.size(); ++i)
            paintDial(g, dials[i], knobs[i]);

        paintMeter(g, 399, processor.meters[0].load(std::memory_order_relaxed), "LEFT");
        paintMeter(g, 605, processor.meters[1].load(std::memory_order_relaxed), "RIGHT");
        label(g, "SIDECHAIN HPF", {222, 166, 114, 12}, 9.f, ink, juce::Font::bold);
        label(g, "DETECTOR", {399, 185, 197, 12}, 9.f, ink, juce::Font::bold);
        label(g, "SIGNAL PATH", {605, 185, 197, 12}, 9.f, ink, juce::Font::bold);
        g.setColour(processor.status.load(std::memory_order_relaxed) ? Colour{0xffacde91} : Colour{0xff94524b});
        g.fillEllipse(785.f, 188.f, 6.f, 6.f);
    }

private:
    void configureSelector(juce::ComboBox& box) {
        box.setJustificationType(juce::Justification::centred);
        addAndMakeVisible(box);
    }

    void paintChassis(Graphics& g) const {
        const Rectangle<float> panel{42.f, 7.f, 1116.f, 233.f};
        g.fillAll(black);
        g.setGradientFill({Colour{0xff8d42cc}, 0.f, 0.f,
                           Colour{0xff5a2685}, 0.f, 247.f, false});
        g.fillRoundedRectangle({1.f, 1.f, 1198.f, 245.f}, 6.f);

        g.setGradientFill({Colour{0xffcb7df8}, 0.f, panel.getY(),
                           Colour{0xffad62ec}, 0.f, panel.getBottom(), false});
        g.fillRoundedRectangle(panel, 4.f);
        if (faceplate.isValid()) {
            Graphics::ScopedSaveState state(g);
            juce::Path trim;
            trim.addRoundedRectangle(panel, 4.f);
            g.reduceClipRegion(trim);
            g.setOpacity(.30f);
            g.drawImage(faceplate, panel, juce::RectanglePlacement::stretchToFit);
        }
        g.setColour(Colour{0xfff0c7ff}.withAlpha(.74f));
        g.drawRoundedRectangle(panel.reduced(.5f), 4.f, 1.f);
        for (int x : {13, 1163}) {
            for (int y : {28, 199}) {
                const Rectangle<float> hole{static_cast<float>(x), static_cast<float>(y), 24.f, 18.f};
                g.setColour(Colour{0xffd49cec}.withAlpha(.62f));
                g.fillRoundedRectangle(hole.translated(0.f, 1.f), 8.f);
                g.setColour(black);
                g.fillRoundedRectangle(hole.reduced(1.1f), 7.f);
            }
        }
        g.setColour(ink.withAlpha(.47f));
        g.drawHorizontalLine(43, 61, 378);
        g.drawHorizontalLine(43, 822, 1139);
        g.drawVerticalLine(382, 51, 223);
        g.drawVerticalLine(818, 51, 223);
        for (int x : {56, 1144}) {
            for (int y : {19, 226}) {
                const auto screw = Rectangle<float>{x - 4.f, y - 4.f, 8.f, 8.f};
                g.setGradientFill({Colour{0xff57505e}, screw.getCentreX(), screw.getY(),
                                   black, screw.getCentreX(), screw.getBottom(), false});
                g.fillEllipse(screw);
                g.setColour(black);
                g.drawLine(x - 2.7f, y + 1.5f, x + 2.7f, y - 1.5f, 1.f);
            }
        }
    }

    static String valueText(Control c, double value) {
        switch (c) {
            case drive: return String(value, 1) + " dB";
            case makeup: return "+" + String(static_cast<int>(value)) + " dB";
            case sensitivity: return value == 0 ? "OFF" : String(-value, 1) + " dB";
            case ratio: return String(ratios[static_cast<std::size_t>(value)]) + ":1";
            case attack: return String(value, 1) + " ms";
            case release: return String(static_cast<int>(value)) + " ms";
            default: return {};
        }
    }

    static void paintDial(Graphics& g, const Dial& dial, const juce::Slider& slider) {
        const auto center = juce::Point<float>{static_cast<float>(dial.x), dial.large ? 126.f : 113.f};
        const float radius = dial.large ? 54.f : 29.f;
        const int tickCount = dial.parameter == ratio ? 6 : dial.large ? 13 : 9;
        for (int tick = 0; tick < tickCount; ++tick) {
            const float angle = rotaryStart + (rotaryEnd - rotaryStart) * tick / (tickCount - 1.f);
            g.setColour(tick % 2 == 0 || dial.parameter == ratio ? glossInk : ink.withAlpha(.65f));
            g.drawLine({center.getPointOnCircumference(radius, angle),
                        center.getPointOnCircumference(radius + (tick % 2 == 0 ? 4.5f : 3.f), angle)}, 1.f);
        }
        const int labelY = dial.large ? 50 : 64;
        label(g, dial.title, {dial.x - 52, labelY, 104, 15}, dial.large ? 12.f : 10.f, glossInk, juce::Font::bold);

        const int valueY = dial.large ? 194 : 146;
        const Rectangle<float> readout(dial.x - 34.f, static_cast<float>(valueY), 68.f, 17.f);
        g.setColour(menuSurface);
        g.fillRoundedRectangle(readout, 2.f);
        g.setColour(ink.withAlpha(.43f));
        g.drawRoundedRectangle(readout, 2.f, .75f);
        label(g, valueText(dial.parameter, slider.getValue()), readout.toNearestInt(), 10.f);
    }

    static void paintMeter(Graphics& g, int x, float value, const String& channel) {
        constexpr int y = 54;
        const Rectangle<float> bezel(static_cast<float>(x), y, 197.f, 119.f);
        g.setColour(Colour{0xff180b25}.withAlpha(.78f));
        g.fillRoundedRectangle(bezel.translated(1.f, 3.f), 7.f);
        g.setGradientFill({Colour{0xff6d687b}, bezel.getX(), bezel.getY(),
                           Colour{0xff17141e}, bezel.getX(), bezel.getBottom(), false});
        g.fillRoundedRectangle(bezel, 7.f);
        g.setColour(Colour{0xff17151f});
        g.fillRoundedRectangle(bezel.reduced(2.f), 5.f);
        const auto glass = bezel.reduced(6.f);
        g.setGradientFill({Colour{0xffffffff}, glass.getCentreX(), glass.getY(),
                           Colour{0xffedebf4}, glass.getCentreX(), glass.getBottom(), false});
        g.fillRoundedRectangle(glass, 3.f);
        g.setColour(Colour{0xffaaa1b6}.withAlpha(.8f));
        g.drawRoundedRectangle(glass, 3.f, 1.f);

        const juce::Point<float> pivot{bezel.getCentreX(), bezel.getBottom() - 17.f};
        constexpr std::array marks{40.f, 30.f, 20.f, 15.f, 10.f, 7.f, 5.f, 3.f, 1.f, 0.f};
        const auto angleAt = [](float db) { return .94f - 1.88f * std::log1p(db * .22f) / std::log1p(40.f * .22f); };
        const auto point = [&](float db, float radius) { return pivot.getPointOnCircumference(radius, angleAt(db)); };

        juce::Path arc;
        for (int i = 0; i <= 50; ++i) {
            const auto p = pivot.getPointOnCircumference(71.f, -.94f + 1.88f * i / 50.f);
            if (i == 0) arc.startNewSubPath(p);
            else arc.lineTo(p);
        }
        g.setColour(ink);
        g.strokePath(arc, juce::PathStrokeType(.8f));
        for (float mark : marks) {
            const auto p = point(mark, 82.f);
            g.setColour(mark >= 20.f ? glossInk : ink);
            g.drawLine({point(mark, 66.f), point(mark, 72.f)}, 1.f);
            label(g, String(static_cast<int>(mark)), {static_cast<int>(p.x) - 9, static_cast<int>(p.y) - 6, 18, 12},
                  9.f, mark >= 20.f ? glossInk : ink, juce::Font::bold);
        }
        label(g, channel + "  GR", {x + 35, y + 58, 127, 13}, 11.f, ink, juce::Font::bold);
        label(g, String(value, 1) + " dB", {x + 48, y + 74, 101, 13}, 10.f, ink);
        g.setColour(needle);
        g.drawLine({pivot, point(juce::jlimit(0.f, 40.f, value), 74.f)}, 1.5f);
        g.setColour(ink);
        g.fillEllipse(Rectangle<float>{6.f, 6.f}.withCentre(pivot));
        g.setColour(white.withAlpha(.23f));
        g.drawHorizontalLine(static_cast<int>(glass.getY()) + 2, glass.getX() + 4, glass.getRight() - 4);
    }

    void timerCallback() override {
        repaint(394, 50, 412, 124);
        repaint(779, 183, 19, 15);
    }

    Processor& processor;
    Hardware hardware;
    std::array<juce::Slider, dials.size()> knobs;
    juce::ComboBox hpfSelector, modeSelector;
    juce::TextButton bypassButton;
    const juce::Image faceplate = juce::ImageCache::getFromMemory(BinaryData::PurpleFaceplate_jpg,
                                                               BinaryData::PurpleFaceplate_jpgSize);
    std::array<std::unique_ptr<State::SliderAttachment>, dials.size()> bindings;
    std::unique_ptr<State::ComboBoxAttachment> hpfBinding;
    std::unique_ptr<juce::ParameterAttachment> topologyBinding, kneeBinding;
    std::unique_ptr<State::ButtonAttachment> bypassBinding;
};

juce::AudioProcessorEditor* Processor::createEditor() { return new Panel(*this); }

} // namespace

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() { return new Processor(); }
