"""Tune one FF hard ratio, keeping the shared detector unmodified."""
import sys, subprocess, json, math
from pathlib import Path
ROOT=Path(__file__).resolve().parent
ratio=int(sys.argv[1])
assert ratio in (8,12,20)
out=ROOT/'measurements'/f'tune_{ratio}'
out.mkdir(parents=True,exist_ok=True)
nominal=10000*(1-1/ratio)
results={}
def evaluate(offset):
    offset=round(float(offset),3)
    if offset in results:return results[offset]
    suite=f'rf_{ratio}_{offset}'
    p=subprocess.run([sys.executable,str(ROOT/'ff_ratio_resistor.py'),suite],cwd=ROOT)
    path=ROOT/'measurements'/suite/'results.json'
    rows=json.loads(path.read_text()) if path.exists() else []
    local=[x for x in rows if x.get('case')=='local']
    valid=[x for x in local if x.get('local_ratio') is not None and x.get('local_ratio')>0]
    complete=p.returncode==0 and len(local)==5 and len(valid)==5
    errors=[100*(x['local_ratio']/ratio-1) for x in valid]
    result={'offset_ohm':offset,'feedback_ohm':nominal+offset,'complete':complete,'max_error_percent':max(map(abs,errors)) if complete else None,'rms_error_percent':math.sqrt(sum(e*e for e in errors)/5) if complete else None,'ratios':[{k:x.get(k) for k in ('target_GR_dB','actual_GR_dB','local_ratio','result')} for x in local]}
    results[offset]=result
    print(json.dumps({'case':'candidate',**result}),flush=True)
    (out/'candidates.json').write_text(json.dumps(list(results.values()),indent=2))
    return result
def best():
    valid=[x for x in results.values() if x['complete']]
    return min(valid,key=lambda x:(x['max_error_percent'],x['rms_error_percent'])) if valid else None
# Initial candidates, then progressively narrower local search.
for off in (0,-50,50):
    evaluate(off)
    if best() and best()['max_error_percent']<=1:break
for step in (50,25,10,5,2,1):
    b=best()
    if b and b['max_error_percent']<=1:break
    center=b['offset_ohm'] if b else 0
    for off in (center-step,center+step):
        if abs(off)<=250:evaluate(off)
        if best() and best()['max_error_percent']<=1:break
b=best()
summary={'ratio':ratio,'target_tolerance_percent':1,'best':b,'status':'within_target' if b and b['max_error_percent']<=1 else 'not_within_target','note':'No production circuit change. Fixed resistor may not correct detector curvature.'}
(out/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({'case':'tuning_summary',**summary}),flush=True)
if not b:sys.exit(2)
