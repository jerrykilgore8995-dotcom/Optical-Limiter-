# Optical Limiter

Manual optical compressor and limiter based on the project’s accepted SPICE circuit. The VST3 runs a separate native nonlinear MNA circuit per audio channel for mono and dual-mono stereo. The detector menu offers Feedback, Feed-forward Hard, and Feed-forward Soft. Feedback uses the hard-knee command and develops its knee in the loop. The existing knee automation controls the feed-forward behavior.

The plugin aligns a full-scale digital sine with +24 dBu at unity input. Threshold spans +30 to −30 dB relative to the detector reference, with OFF at the fully counterclockwise position. The digital detector compensates for Input trim so drive and gain reduction can be set independently. Input trim tops out at +18 dB. Clipping headroom also depends on incoming peak level. The fitted photocell K2 and K4 correction runs at twice the catalog-fit amplitude. Output controls the circuit makeup stage. Both meters report optical gain reduction before makeup from the solved LDR conductance. At host rates below 88.2 kHz the engine runs at twice the host rate with a linear-phase half-band FIR in each direction. At 44.1 kHz the circuit solves at 88.2 kHz and the plugin reports 96 samples of latency. The Signal Path switch delays dry audio by that same amount for host alignment.

In stereo, an independent solver can assist the callback by processing the right channel while the callback processes the left. It publishes completed audio and circuit history in sections. The callback can finish from the latest published section without waiting and reports the same signal latency.

## Build and test

Requires CMake 3.22 or newer, Python 3.10 or newer, a C++20 compiler, and ngspice to export the operating point. A pre-exported operating point can be supplied with `-DOPTICAL_SPICE_SEED=/path/to/operating-point.raw`. CMake fetches JUCE 8.0.10 and generates native stamps from the reduced sidechain deck. The pre-exported seed must come from that deck. Generate a standalone ngspice input with `python manual/build_fast_deck.py build/Optical-Limiter-Fast.cir`. The generated file refers to the vendor libraries in this checkout.

```sh
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

The Windows x64 GitHub Actions job builds with clang-cl and targets CPUs with AVX2 and FMA. It exports the operating point, builds the plugin, loads the VST3 in a test host, and renders audio in mono and stereo. It publishes the VST3 bundle and a panel capture.

Use `native/compare_spice.py` with the built `OpticalEngineTest` executable to compare native output, gain reduction, phase, and THD with a fresh ngspice run of the full vendor reference deck. `--short` runs both detector branches across the three modes. The 5% flag evaluates fundamental amplitude and gain reduction. Phase and THD are reported separately.

## Photocell harmonics

The native photocell fit uses the terminal RMS voltage and conductance. Digitized Figure 3(a) traces provide the total THD target as a function of these two inputs. The fit allocates the cell correction to K2 and a smaller K4. The chart does not resolve the harmonic split. Outside its resistance range the nearest trace is used. Voltages above the printed range are held at the last point. The 10 kHz K2 reference applies to the deeply lit region outside the measured resistance range.

The main audio circuit solves the original odd current law. The even correction follows the photocell current through a linearized feedback and output network. This avoids moving the existing odd harmonics through interaction with the added even current at high drive.

## Digital sidechain

The average detector uses an ideal full-wave rectifier and the peak detector an ideal peak follower with its original discharge time. The threshold, ratio curves and attack/release equations remain in the control path. A small steady-state peak calibration preserves the measured response of the reference at high frequencies. The control uses a cycle mean while a tone is stable. The unfiltered audio signal continues to drive detection. Large level changes and loss of periodicity restore the timing calculation directly.

The original SPICE reference and the generated reference deck retain the component-level behavior for comparison.
