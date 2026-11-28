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
    linear_cell = 'Bldr gr sum I={V(gr,sum)*cg(V(cell))}'
    if audio.count(linear_cell) != 1:
        raise ValueError('Missing reference photocell stamp')
    # Figure 3(a) constrains the low-voltage even floor at 30k and 3k.
    # Add a small polarity asymmetry to the reduced deck's photocell only.
    even_cell = """.param LDR_EVEN_30K=35u LDR_EVEN_3K=45u LDR_EVEN_KNEE=0.10
.func ldrmix(g) {min(1,max(0,(g-1/30000)/(1/3000-1/30000)))}
.func ldrshape(v) {v*v/sqrt(v*v+LDR_EVEN_KNEE*LDR_EVEN_KNEE)}
Bldr gr sum I={cg(V(cell))*(V(gr,sum)+(LDR_EVEN_30K+(LDR_EVEN_3K-LDR_EVEN_30K)*ldrmix(cg(V(cell))))*ldrshape(V(gr,sum)))}"""
    audio = audio.replace(linear_cell, even_cell)
    sidechain = (ROOT / 'sidechain_fast.cir.inc').read_text(encoding='utf-8')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(audio + sidechain + '\n.op\n.tran 0.5u 200m\n.end\n',
                           encoding='utf-8')


if __name__ == '__main__':
    main()
