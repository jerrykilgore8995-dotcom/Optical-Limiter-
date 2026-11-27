"""FF soft knee: measure actual transfer relative to measured hard threshold."""
from pathlib import Path
p=Path(__file__).with_name('measure.py')
exec(compile(p.read_text().split("if SUITE=='audio' or")[0],str(p),'exec'))
ratio=int(SUITE.split('_')[-1])
assert ratio in (2,3,4,8,12,20)
# Fit hard knee threshold using established 5-10dB region, no circuit modifications.
a=1-1/ratio
def sample(level,soft):
    tag=f"{'soft' if soft else 'hard'}_{level:.5f}"
    params={'FF':1,'SOFT':soft,'RATIO':ratio,'SENS':30,'AMP':peak(level)}
    for duration in (.22,.45,.9):
        x=run(tag+'_'+str(duration),params,seconds=duration,step=10e-6,signals=('pre','gr','out','rect','timing','sense','cell','det','logv','pbuf','cmd','over','knee'))
        st=steady(x)
        if st['settled']:break
    if not st['settled']:raise RuntimeError('Unsettled '+tag)
    inp=tone(x,1000,1)['peak_V'];out=tone(x,1000,2)['peak_V']
    t=x[-1,0]-.02+np.arange(2000)*.02/2000
    over=np.interp(t,x[:,0],x[:,12]);knee=np.interp(t,x[:,0],x[:,13])
    row={'case':'point','status':'OK','ratio':ratio,'soft':soft,'input_dBu':level,'pre_dBu':db(inp/math.sqrt(2)/.775),'out_dBu':db(out/math.sqrt(2)/.775),'GR_dB':db(inp/out),'detector_over_mean_dB':float(np.mean(over)),'knee_mean_dB':float(np.mean(knee)),**st}
    save(row);return row
base=sample(-45,0)['GR_dB']
guess=-31.04 if ratio<8 else -30.6
h=[sample(guess+g/a,0) for g in (5,10)]
slope=(h[1]['GR_dB']-h[0]['GR_dB'])/(h[1]['pre_dBu']-h[0]['pre_dBu'])
threshold=h[0]['pre_dBu']-(h[0]['GR_dB']-base)/slope
source_offset=np.mean([p['pre_dBu']-p['input_dBu'] for p in h])
save({'case':'threshold','status':'OK','ratio':ratio,'threshold_pre_dBu':threshold,'hard_ratio_fit':1/(1-slope),'baseline_GR_dB':base})
def ideal(x):
    return a*(0 if x<=-5 else x if x>=5 else (x+5)**2/20)
# One dB input increments through knee; add endpoints and established region.
points=[]
for d in [-8,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,8,10,15,20,25,30,35,40]:
    level=float(threshold-source_offset+d)
    if level>18:continue
    row=sample(level,1)
    x=row['pre_dBu']-threshold
    g=row['GR_dB']-base
    point={'case':'curve','status':'OK','ratio':ratio,'input_relative_threshold_dB':x,'GR_dB':g,'output_relative_threshold_dB':x-g,'ideal_GR_dB':ideal(x),'GR_error_dB':g-ideal(x),'ideal_local_ratio':1/(1-a*min(1,max(0,(x+5)/10)))}
    points.append(point);save(point)
for l,r in zip(points,points[1:]):
    dx=r['input_relative_threshold_dB']-l['input_relative_threshold_dB']
    dy=r['output_relative_threshold_dB']-l['output_relative_threshold_dB']
    save({'case':'local_interval','status':'OK','ratio':ratio,'input_interval_dB':[l['input_relative_threshold_dB'],r['input_relative_threshold_dB']],'GR_interval_dB':[l['GR_dB'],r['GR_dB']],'local_ratio':dx/dy if abs(dy)>1e-9 else None,'ideal_ratio':dx/(dx-ideal(r['input_relative_threshold_dB'])+ideal(l['input_relative_threshold_dB']))})
# Actual gain reduction targets, use interpolated crossing then measured symmetric local slope.
for target in range(1,16):
    pair=next(((l,r) for l,r in zip(points,points[1:]) if l['GR_dB']<=target<=r['GR_dB'] and r['GR_dB']>l['GR_dB']),None)
    if pair is None:continue
    l,r=pair
    d=l['input_relative_threshold_dB']+(r['input_relative_threshold_dB']-l['input_relative_threshold_dB'])*(target-l['GR_dB'])/(r['GR_dB']-l['GR_dB'])
    level=float(threshold-source_offset+d)
    if level+.25>18:continue
    lo=sample(level-.25,1);hi=sample(level+.25,1)
    dx=hi['pre_dBu']-lo['pre_dBu'];dy=hi['out_dBu']-lo['out_dBu']
    save({'case':'GR_local','status':'OK','ratio':ratio,'target_GR_dB':target,'actual_midpoint_GR_dB':(lo['GR_dB']+hi['GR_dB'])/2-base,'local_ratio':dx/dy if abs(dy)>1e-9 else None})
