"""Assemble the reduced sidechain with the reference deck's audio path."""
import argparse
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
MARKER = '* Sidechain source selector.'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    reference = (ROOT / 'manual.cir').read_text(encoding='utf-8')
    audio, separator, _ = reference.partition(MARKER)
    if not separator:
        raise ValueError(f'Missing section marker {MARKER}')
    audio = re.sub(r'(?im)^\.include\s+(\S+)',
                   lambda m: f'.include "{(ROOT / m[1]).resolve().as_posix()}"', audio)
    sidechain = (ROOT / 'sidechain_fast.cir.inc').read_text(encoding='utf-8')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(audio + sidechain + '\n.op\n.tran 0.5u 200m\n.end\n',
                           encoding='utf-8')


if __name__ == '__main__':
    main()
