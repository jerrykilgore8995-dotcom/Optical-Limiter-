# Stage 1 — input receiver
Status: provisional SPICE topology only. No simulation results or native MNA implementation yet.

Digital gain -20 to +20 dB, default 0; split into +/- half amplitude. Four complementary emitter followers on +/-24 V, differential OPA2134 channel on +/-15 V. Initial 10k output load. Ideal bias offsets in receiver.cir are placeholders; replace with physical bias network and tune for 3mA per device. Verify that neither device cuts off throughout intended Class A operating range.

Model sources:
- https://www.ti.com/product/OPA2134 (OPAx134 PSpice model)
- https://www.onsemi.com/ (KSC1845 and KSA992 device model downloads)
Vendor model files are not yet imported. Check pin order and licensing before use. Do not substitute generic models silently.

Validation:
- DC nodes, currents, differential gain, common-mode behavior, output impedance and headroom.
- SPICE reference AC and transient; loads 600, 2k, 10k, 47k.
- Native MNA at 44100 Hz, 512 samples; warm-start DC solution; no corrective EQ.
- Compare magnitude in linear units (5% relative tolerance, approximately +/-0.45dB), phase via complex transfer error and absolute degrees (avoid relative division near zero), THD/harmonic amplitudes in linear units with explicit numerical floor.
- Frequencies 100, 1000, 5000, 10000 Hz; levels -18,-12,-8,-4,0,4,8,12,18 dBu. Only report harmonics below Nyquist for native sampled measurements; use high-rate reference for higher orders. Check aliasing independently.
- Report actual analog response, not an imposed flat target. Investigate additional native rolloff at 20kHz.
- Performance target whole stereo plugin 2-4ms/512 on Ryzen 8745HS; measure rather than claim.
- Internal dBFS/dBu calibration not yet finalized. Reference peak amplitude = sqrt(2)*0.775*10^(dBu/20).
No optical section, GUI or VST3 at this stage.
