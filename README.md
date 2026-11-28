# Optical Limiter

Manual optical compressor and limiter based on the project’s accepted SPICE circuit. The VST3 runs a separate native nonlinear MNA circuit per audio channel for mono and dual-mono stereo. The detector menu offers Feedback, Feed-forward Hard, and Feed-forward Soft. Feedback uses the hard-knee command and develops its knee in the loop. The existing knee automation controls the feed-forward behavior.

The plugin aligns a full-scale digital sine with +24 dBu at unity input. Threshold is the circuit’s sensitivity in dB below its reference. The zero-sensitivity position switches gain reduction off. Output controls the circuit makeup stage. Both meters report optical gain reduction before makeup from the solved LDR conductance. At host rates below 88.2 kHz the engine runs at twice the host rate with a linear-phase half-band FIR in each direction. At 44.1 kHz the circuit solves at 88.2 kHz and the plugin reports 96 samples of latency. The Signal Path switch delays dry audio by that same amount for host alignment.

## Build and test

Requires CMake 3.22 or newer, Python 3.10 or newer, a C++20 compiler, and ngspice to export the operating point. A pre-exported operating point can be supplied with `-DOPTICAL_SPICE_SEED=/path/to/operating-point.raw`. CMake fetches JUCE 8.0.10 and generates native stamps from the reduced sidechain deck. The pre-exported seed must come from that deck. Generate a standalone ngspice input with `python manual/build_fast_deck.py build/Optical-Limiter-Fast.cir`. The generated file refers to the vendor libraries in this checkout.

```sh
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The Windows x64 GitHub Actions job builds with clang-cl and targets CPUs with AVX2 and FMA. It exports the operating point, builds the plugin, loads the VST3 in a test host, and renders audio in mono and stereo. It publishes the VST3 bundle and a panel capture.

Use `native/compare_spice.py` with the built `OpticalEngineTest` executable to compare native output, gain reduction, phase, and THD with a fresh ngspice run of the full vendor reference deck. `--short` runs both detector branches across the three modes. The 5% flag evaluates fundamental amplitude and gain reduction. Phase and THD are reported separately.

## Reduced sidechain

The BJT receiver, optical audio stage, makeup and output amplifiers retain their manufacturer models. The sidechain uses buffered RC stages, a precision rectifier and a current-limited peak-charge model with both edge and crest charging. The native solver eliminates these sidechain equations algebraically. It solves the remaining 129 audio and control unknowns with the detector's derivative in the Newton matrix. Trapezoidal detector history advances only when a circuit step is accepted.

The peak branch integrates four intervals per native step with causal interpolation to resolve short charge pulses. The audio rate and oversampling are shared with the average detector. Detector selection by ratio, HPF, timing and the photocell conductance curve follow the reduced SPICE deck. See `manual/manual.cir` for the full vendor reference.
