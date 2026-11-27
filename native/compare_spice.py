"""Compare the locked deck and the native solve at identical input and settings."""
import argparse
import io
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RATE = 44100
SIGNALS = ('pre', 'gr', 'out', 'cell', 'det', 'timing')


def db(value):
    return 20 * np.log10(max(abs(value), 1e-30))


def steady(wave, frequency):
    cycles = max(2, min(20, int(wave[-1, 0] * frequency / 3)))
    period = cycles / frequency
    t = wave[-1, 0] - period + np.arange(2048) * period / 2048
    output = [np.interp(t, wave[:, 0], wave[:, i]) for i in range(1, 4)]
    basis = np.exp(-2j * np.pi * frequency * t)
    fundamental = [2 * np.mean(s * basis) for s in output]
    harmonic_count = min(20, int((RATE / 2 - 1) / frequency))
    harmonics = [2 * np.mean(output[2] * basis ** k) for k in range(2, harmonic_count + 1)]
    return {
        'pre_peak_V': abs(fundamental[0]),
        'output_peak_V': abs(fundamental[2]),
        'gain_reduction_dB': db(abs(fundamental[0]) / max(abs(fundamental[1]), 1e-30)),
        'phase_deg': float(np.angle(fundamental[2] / fundamental[0], deg=True)),
        'thd_percent': 100 * np.linalg.norm(harmonics) / max(abs(fundamental[2]), 1e-30),
    }


def spice_case(ff, soft, ratio, sensitivity, level, seconds, frequency, work):
    source = ROOT / 'manual/manual.cir'
    deck = re.sub(r'(?im)^\.include\s+(\S+)',
                  lambda m: f'.include "{(source.parent / m[1]).resolve().as_posix()}"',
                  source.read_text())
    settings = dict(FF=ff, SOFT=soft, RATIO=ratio, SENS=sensitivity,
                    AMP=2 ** .5 * .775 * 10 ** (level / 20), FREQ=frequency)
    for key, value in settings.items():
        deck = re.sub(rf'\b{key}=\S+', f'{key}={value}', deck, count=1)
    deck = re.sub(r'(?im)^\.(op|tran|end)\b[^\n]*', '', deck)
    signals = ' '.join(f'v({s})' for s in SIGNALS)
    deck += f'''
.control
set wr_singlescale
set wr_vecnames
tran 10u {seconds} 0 10u
wrdata wave.txt {signals}
quit
.endc
.end
'''
    (work / 'case.cir').write_text(deck)
    result = subprocess.run([os.environ.get('NGSPICE', 'ngspice'), '-b', 'case.cir'],
                            cwd=work, capture_output=True, text=True, timeout=300, check=True)
    wave = np.loadtxt(work / 'wave.txt', skiprows=1)
    if re.search(r'timestep too small|simulation aborted|error:', result.stderr, re.I) or wave[-1, 0] < .999 * seconds:
        raise RuntimeError(result.stdout + result.stderr)
    return steady(wave, frequency)


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('native_executable', type=Path)
    args.add_argument('output', type=Path)
    args.add_argument('--short', action='store_true', help='Limit the run to one drive level per mode and detector')
    options = args.parse_args()
    exe = options.native_executable.resolve()
    rows = []
    cases = [(ff, soft, ratio, 30, level, 1000)
             for ff in (1, 0) for soft in (0, 1)
             for ratio in ((4, 12) if options.short else (2, 3, 4, 8, 12, 20))
             for level in ((-16,) if options.short else (-36, -26, -16, -6, 4))]
    if not options.short:
        cases += [(1, 0, 4, 0, level, frequency)
                  for frequency in (100, 1000, 5000, 10000)
                  for level in (-18, -12, -8, -4, 0, 4, 8, 12, 18)]
    for ff, soft, ratio, sens, level, frequency in cases:
        seconds = .22
        with tempfile.TemporaryDirectory() as directory:
            reference = spice_case(ff, soft, ratio, sens, level, seconds, frequency, Path(directory))
        result = subprocess.run([str(exe), *map(str, (ff, soft, ratio, sens, level, seconds, RATE, frequency))],
                                capture_output=True, text=True, timeout=900)
        row = dict(ff=bool(ff), soft=bool(soft), ratio=ratio, sensitivity=sens,
                   input_dBu=level, frequency_Hz=frequency, spice=reference,
                   native_status=result.stderr.strip())
        if result.stdout.startswith('time,'):
            wave = np.loadtxt(io.StringIO(result.stdout), delimiter=',', skiprows=1)
            native = steady(wave, frequency)
            row['native'] = native
            row['output_gain_error_dB'] = db(native['output_peak_V'] / reference['output_peak_V'])
            row['gain_reduction_error_dB'] = native['gain_reduction_dB'] - reference['gain_reduction_dB']
            row['within_5_percent'] = result.returncode == 0 and abs(
                native['output_peak_V'] / reference['output_peak_V'] - 1) <= .05 and abs(
                row['gain_reduction_error_dB']) <= max(.10, abs(reference['gain_reduction_dB']) * .05)
        else:
            row['within_5_percent'] = False
        rows.append(row)
        options.output.write_text(json.dumps(rows, indent=2, allow_nan=False) + '\n')
        print(json.dumps(row, allow_nan=False), flush=True)
    return 0 if all(r['within_5_percent'] for r in rows) else 1


if __name__ == '__main__':
    raise SystemExit(main())
