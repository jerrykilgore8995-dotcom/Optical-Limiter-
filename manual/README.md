# Manual optical limiter
Prototype netlist: manual.cir. No Auto circuitry.

This is NOT a finished VST3 or validated MNA model. Vendor libraries are missing; verify original model pinouts and model names. The bias offsets, ideal switches, log/threshold/knee/ratio blocks, timing steering, optical dynamics and inverse LED shaping are provisional behavioral fixtures, not physical component implementations. Replace with agreed physical circuits before declaring circuit-level completion.

Locked design: +/-24V BJT receiver, +/-15V OPA2134 stages; WIMA MKP2 coupling caps matched within 1%; 10k optical input and feedback, 100pF; VTL5C1; 0–30dB stepped makeup; 47 ohm output. FF/FB, HPF off–200Hz/10Hz steps, 220pF pre-rectifier, 1N4148W precision full-wave rectification; average 10k/100n with 1.571 gain for 2/3/4:1, peak 47k/100n for 8/12/20:1; hard/soft 10dB knee, continuous 50k threshold pot, 100mV/dB; manual timing 1–40ms and 10ms–3s nominal; actual timings may vary. Current sink 2N2222A, 100 ohm sense, nominal 30mA cap. Check actual current limit and driver dissipation.

Validation order: import verified manufacturer models; DC bias/headroom; each audio stage FR/phase/load/harmonics; rectifier polarity/accuracy and recovery; optical static/dynamic curves; physical threshold/ratio/timing/driver; closed-loop stability; native MNA comparison at 44.1kHz; 512-sample performance; only then VST3.
Use 5% reference tolerance where meaningful; report phase and near-zero harmonic errors without unstable relative division. Measure ratio segments 0–10,10–20,20–30dB GR, both knees and topologies; actual attack/release at 3/6/10/20/30dB GR. No corrective EQ or fake saturation. Record unreachable operating points rather than fabricate results.
