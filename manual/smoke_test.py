"""Run a real transient and verify completion, not merely process exit."""
from pathlib import Path
import subprocess, re, math
root=Path(__file__).resolve().parent
out=root/'results'
out.mkdir(exist_ok=True)
source=(root/'manual.cir').read_text()
# Use an explicit control block: batch ngspice otherwise has no requested
# print, plot or rawfile output in the original smoke-test invocation.
deck=re.sub(r'(?im)^\.end\s*$',""".control
set wr_singlescale
set wr_vecnames
tran 0.5u 200m
wrdata results/manual_transient.dat v(pre) v(gr) v(out) v(rect) v(timing) v(sense)
quit
.endc
.end
""",source)
deck=re.sub(r'(?im)^\.tran[^\n]*\n','',deck)
# Remove only the original .tran directive (control command has no dot).
p=root/'smoke_generated.cir'
p.write_text(deck)
data=out/'manual_transient.dat'
data.unlink(missing_ok=True)
try:
    proc=subprocess.run(['ngspice','-b',str(p)],cwd=root,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=240)
    (out/'manual.log').write_text(proc.stdout)
    print(proc.stdout)
    if proc.returncode or re.search(r'error:|timestep too small|fatal|simulation\s+aborted',proc.stdout,re.I):
        raise RuntimeError('ngspice reported an error')
    if not data.exists():
        raise RuntimeError('No transient output produced')
    last=None
    count=0
    with data.open() as stream:
        for line in stream:
            fields=line.split()
            if not fields: continue
            if fields[0].lower()=='time': continue
            row=list(map(float,fields))
            if not all(map(math.isfinite,row)): raise RuntimeError('Nonfinite transient output')
            last=row; count+=1
    if last is None or last[0]<0.19999:
        raise RuntimeError(f'Transient did not reach 200 ms; last row={last}')
    print(f'PASS: {count} finite transient samples, final time {last[0]:.9f} s')
finally:
    p.unlink(missing_ok=True)
