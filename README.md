# Optical Limiter

Manual optical compressor and limiter based on the project’s accepted SPICE circuit. The VST3 runs a separate native nonlinear MNA circuit per audio channel for mono and dual-mono stereo. It uses the existing FF/FB detector selector, hard and soft knees, and six ratio positions.

The plugin aligns a full-scale digital sine with +24 dBu at unity input. Threshold is the circuit’s sensitivity in dB below its reference. The zero-sensitivity position switches gain reduction off. Output controls the circuit makeup stage. Both meters report optical gain reduction before makeup from the solved LDR conductance. At host rates below 88.2 kHz the engine runs at twice the host rate with a linear-phase half-band FIR in each direction. At 44.1 kHz the circuit solves at 88.2 kHz and the plugin reports 96 samples of latency. The Signal Path switch delays dry audio by that same amount for host alignment.

## Build and test

Requires CMake 3.22 or newer, Python 3.10 or newer, a C++20 compiler, and ngspice to export the operating point. A pre-exported operating point can be supplied with `-DOPTICAL_SPICE_SEED=/path/to/operating-point.raw`. CMake fetches JUCE 8.0.10 and generates the full native stamps from `manual/manual.cir` and its included vendor models.

```sh
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The Windows x64 GitHub Actions job builds with clang-cl and targets CPUs with AVX2 and FMA. It exports the operating point, builds the plugin, loads the VST3 in a test host, and renders audio in mono and stereo. It publishes the VST3 bundle and a panel capture.

Use `native/compare_spice.py` with the built `OpticalEngineTest` executable to compare native output, gain reduction, phase, and THD with a fresh ngspice run of the same deck. `--short` runs both detector branches across all four modes. The 5% flag evaluates fundamental amplitude and gain reduction. Phase and THD are reported separately.

## Playback status

The engine reuses constant conductance and source stamps and the invariant part of sparse elimination. It preallocates trapezoidal charge histories. A few detector equations screen Newton trials when their individual residual already exceeds the full convergence merit. Trials that pass the screen still run through the complete circuit, including the manufacturer device models.

The plugin is an offline-render preview. The full circuit’s driven solve times exceed the 11.61 ms per 512-sample block budget at 44.1 kHz in development tests. The packaged test host logs median and maximum block times with both stereo channels driven. It compares the threshold-off gain at 1 and 20 kHz at 44.1 kHz against the SPICE through response. Playback on the Ryzen 7 8745HS remains unverified. Auto and EPR are outside the manual circuit and are absent from the plugin.
