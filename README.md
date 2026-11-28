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

The BJT receiver, optical audio stage, makeup and output amplifiers retain their manufacturer models. The reduced sidechain uses buffered RC stages, a piecewise precision rectifier and a current-limited peak-charge model. It retains the detector selection by ratio, HPF, timing states and photocell conductance curve. SPICE and the VST3 use the same reduced deck. `manual/manual.cir` is the full vendor reference for comparison. The reduced model approximates the vendor stages' switching transients and low-level distortion.

The packaged VST3 test reports median and maximum block times for mono, driven stereo, and threshold-off response. It also measures the 1 to 20 kHz through response and verifies bypass delay at 44.1 and 88.2 kHz.
