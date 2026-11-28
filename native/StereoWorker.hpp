#pragma once

#include "RateAdaptedEngine.hpp"

#include <atomic>
#include <system_error>
#include <thread>
#include <vector>

#if defined(__SSE2__) || defined(_M_X64)
#include <xmmintrin.h>
#endif

namespace optical {

// An independent right-channel solver can publish completed sections while
// the callback solves the left channel. The callback can resume the right
// channel from its newest published state, even if this thread is preempted.
class StereoWorker {
public:
    ~StereoWorker() { stop(); }

    void prepare(double rate, const Settings& settings, int preferredFrames) {
        stop();
        if (preferredFrames < minimumFrames || !engine.prepare(rate, settings)) return;
        const int capacity = std::min(preferredFrames, maximumFrames);
        samples.resize(capacity);
        checkpoints.resize((capacity + sectionSize - 1) / sectionSize);
        phase.store(Phase::idle, std::memory_order_relaxed);
        try {
            worker = std::thread([this] { run(); });
        } catch (const std::system_error&) {
            phase.store(Phase::stopped, std::memory_order_relaxed);
        }
    }

    void stop() {
        phase.store(Phase::stopped, std::memory_order_release);
        phase.notify_one();
        if (worker.joinable()) worker.join();
    }

    bool submit(const RateAdaptedEngine& source, const Settings& next, const float* input, int frames) {
        if (frames < minimumFrames || static_cast<std::size_t>(frames) > samples.size()
            || phase.load(std::memory_order_acquire) != Phase::idle) return false;

        source.save(initial);
        settings = next;
        std::copy_n(input, frames, samples.begin());
        length = frames;
        published.store(0, std::memory_order_relaxed);
        phase.store(Phase::queued, std::memory_order_release);
        phase.notify_one();
        return true;
    }

    // Called exactly once for each successful submission. Applying a
    // checkpoint never reads the solver or any unpublished sample.
    int collect(RateAdaptedEngine& destination, float* output) {
        auto current = phase.load(std::memory_order_acquire);
        while (current == Phase::queued || current == Phase::running) {
            const auto replacement = current == Phase::queued ? Phase::idle : Phase::abandoned;
            if (phase.compare_exchange_strong(current, replacement, std::memory_order_acq_rel))
                break;
        }

        const int done = published.load(std::memory_order_acquire);
        if (done > 0) {
            destination.restore(checkpoints[(done - 1) / sectionSize]);
            std::copy_n(samples.begin(), done, output);
        }
        if (current == Phase::ready) phase.store(Phase::idle, std::memory_order_release);
        return done;
    }

private:
    enum class Phase { stopped, idle, queued, running, ready, abandoned };
    static constexpr int sectionSize = 48;
    static constexpr int minimumFrames = 128;
    static constexpr int maximumFrames = 2048;

    void run() {
#if defined(__SSE2__) || defined(_M_X64)
        _mm_setcsr(_mm_getcsr() | 0x8040); // The callback also flushes subnormals.
#endif
        for (;;) {
            auto current = phase.load(std::memory_order_acquire);
            if (current == Phase::stopped) return;
            if (current != Phase::queued) {
                phase.wait(current, std::memory_order_acquire);
                continue;
            }
            if (!phase.compare_exchange_strong(current, Phase::running, std::memory_order_acq_rel))
                continue;

            engine.restore(initial);
            engine.configure(settings);
            for (int start = 0; start < length; start += sectionSize) {
                if (phase.load(std::memory_order_relaxed) != Phase::running) break;
                const int count = std::min(sectionSize, length - start);
                engine.process(samples.data() + start, count);
                engine.save(checkpoints[start / sectionSize]);
                published.store(start + count, std::memory_order_release);
            }

            current = Phase::running;
            if (!phase.compare_exchange_strong(current, Phase::ready, std::memory_order_acq_rel)
                && current == Phase::abandoned) {
                phase.compare_exchange_strong(current, Phase::idle, std::memory_order_acq_rel);
            }
        }
    }

    RateAdaptedEngine engine;
    RateAdaptedEngine::State initial{};
    Settings settings;
    std::vector<RateAdaptedEngine::State> checkpoints;
    std::vector<float> samples;
    std::atomic<Phase> phase{Phase::stopped};
    std::atomic<int> published{0};
    int length = 0;
    std::thread worker;
};

} // namespace optical
