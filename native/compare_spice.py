"""Compare the locked deck and the native solve at identical input and settings."""
import argparse
from concurrent.futures import ThreadPoolExecutor
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
    # Fit at the original sample times. Reinterpolating the native wave
    # onto a finer grid would attenuate high-frequency tones and harmonics.
    t = wave[:, 0]
    boundary = (t[1:] + t[:-1]) * .5
    left = np.maximum(np.r_[t[0], boundary], t[-1] - period)
    right = np.r_[boundary, t[-1]]
    mask = right > left
    weight = np.sqrt(right[mask] - left[mask])
    harmonic_count = min(20, int((RATE / 2 - 1) / frequency))
    angle = 2 * np.pi * frequency * t[mask, None] * np.arange(1, harmonic_count + 1)
    design = np.column_stack((np.ones(len(weight)), np.cos(angle), np.sin(angle)))
    fitted = np.linalg.lstsq(design * weight[:, None], wave[mask, 1:4] * weight[:, None], rcond=None)[0]
    phasors = fitted[1:1 + harmonic_count] - 1j * fitted[1 + harmonic_count:]
    fundamental = phasors[0]
    return {
        'pre_peak_V': float(abs(fundamental[0])),
        'output_peak_V': float(abs(fundamental[2])),
        'gain_reduction_dB': float(db(abs(fundamental[0]) / max(abs(fundamental[1]), 1e-30))),
        'phase_deg': float(np.angle(fundamental[2] / fundamental[0], deg=True)),
        'thd_percent': float(100 * np.linalg.norm(phasors[1:, 2]) / max(abs(fundamental[2]), 1e-30)),
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
    args.add_argument('--jobs', default=1, type=int, help='Independent cases to run concurrently')
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
    def compare(case):
        ff, soft, ratio, sens, level, frequency = case
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
            row['within_5_percent'] = bool(result.returncode == 0 and abs(
                native['output_peak_V'] / reference['output_peak_V'] - 1) <= .05 and abs(
                row['gain_reduction_error_dB']) <= max(.10, abs(reference['gain_reduction_dB']) * .05))
        else:
            row['within_5_percent'] = False
        return row

    with ThreadPoolExecutor(max_workers=options.jobs) as workers:
        for row in workers.map(compare, cases):
            rows.append(row)
            options.output.write_text(json.dumps(rows, indent=2, allow_nan=False) + '\n')
            print(json.dumps(row, allow_nan=False), flush=True)
    return 0 if all(r['within_5_percent'] for r in rows) else 1


if __name__ == '__main__':
    raise SystemExit(main())
