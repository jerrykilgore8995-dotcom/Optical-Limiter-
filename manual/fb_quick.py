"""Compact FB hard/soft static survey; locked production circuit, no tuning."""
from pathlib import Path
p=Path(__file__).with_name('measure.py')
exec(compile(p.read_text().split("if SUITE=='audio' or")[0],str(p),'exec'))
_,mode,rs=SUITE.split('_')
ratio=int(rs); soft=int(mode=='soft')
assert mode in ('hard','soft') and ratio in (2,3,4,8,12,20)
cache={}
def sample(level):
    level=round(float(level),5)
    if level in cache:return cache[level]
    params={'FF':0,'SOFT':soft,'RATIO':ratio,'SENS':30,'AMP':peak(level),'FREQ':1000}
    for duration in (.22,.45,.9):
        x=run(f'level_{level}_{duration}',params,seconds=duration,step=10e-6)
        st=steady(x)
        if st['settled']:break
    a=tone(x,1000,1)['peak_V'];b=tone(x,1000,2)['peak_V']
    row={'case':'point','status':'OK' if st['settled'] else 'UNSETTLED','mode':mode,'ratio':ratio,'input_dBu':level,'pre_dBu':db(a/math.sqrt(2)/.775),'output_dBu':db(b/math.sqrt(2)/.775),'GR_dB':db(a/b),**st}
    cache[level]=row;save(row);return row
base=sample(-45)['GR_dB']
# Actual input sweep, not assuming feedback knee width or ratio.
for level in [-38,-36,-34,-32,-31,-30,-29,-28,-26,-24,-21,-17,-12,-7,-2,4,10,18]:
    sample(level)
points=sorted((p for p in cache.values() if p['settled']),key=lambda p:p['pre_dBu'])
def cross(target):
    for a,b in zip(points,points[1:]):
        ga,gb=a['GR_dB']-base,b['GR_dB']-base
        if ga<=target<=gb and gb>ga:
            w=(target-ga)/(gb-ga)
            return a['pre_dBu']+w*(b['pre_dBu']-a['pre_dBu']),a['output_dBu']+w*(b['output_dBu']-a['output_dBu']),a['input_dBu']+w*(b['input_dBu']-a['input_dBu'])
    return None
# True hard threshold reference is provisional from earlier detector tests;
# use 0.05 dB measured onset proxy for first interval, not a local 9-10 slope.
for lo,hi in [(0,10),(10,20),(20,30)]:
    a,b=cross(.05 if lo==0 else lo),cross(hi)
    save({'case':'segment','status':'OK','mode':mode,'ratio':ratio,'GR_interval_dB':[lo,hi],'onset_proxy_dB':.05 if lo==0 else None,'effective_ratio':(b[0]-a[0])/(b[1]-a[1]) if a and b and abs(b[1]-a[1])>1e-9 else None,'method':'interpolated measured crossings'})
for target in [1,3,5,10,15,20,25,30]:
    c=cross(target)
    if c is None or c[2]+.5>18:continue
    a,b=sample(c[2]-.5),sample(c[2]+.5)
    if not a['settled'] or not b['settled']:continue
    dx=b['pre_dBu']-a['pre_dBu'];dy=b['output_dBu']-a['output_dBu']
    save({'case':'local','status':'OK','mode':mode,'ratio':ratio,'target_GR_dB':target,'actual_midpoint_GR_dB':(a['GR_dB']+b['GR_dB'])/2-base,'local_ratio':dx/dy if abs(dy)>1e-9 else None,'input_bracket_dBu':[a['input_dBu'],b['input_dBu']]})
save({'case':'summary','status':'OK','mode':mode,'ratio':ratio,'max_measured_GR_dB':max(p['GR_dB']-base for p in points),'baseline_GR_dB':base,'unsettled_count':sum(not p['settled'] for p in cache.values())})
