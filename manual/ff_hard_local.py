"""FF hard knee only. Full circuit, no calibration edits."""
from pathlib import Path
p=Path(__file__).with_name('measure.py')
exec(compile(p.read_text().split("if SUITE=='audio' or")[0],str(p),'exec'))
ratio=int(SUITE.split('_')[-1])
assert ratio in (2,3,4,8,12,20)
data={}
def sample(level):
    level=round(float(level),5)
    if level in data:return data[level]
    tag='level_'+str(level)
    params={'FF':1,'SOFT':0,'RATIO':ratio,'SENS':30,'AMP':peak(level)}
    for duration in [.22,.45,.9]:
        x=run(tag+'_'+str(duration),params,seconds=duration,step=10e-6)
        st=steady(x)
        if st['settled']:break
    if not st['settled']:raise RuntimeError('Unsettled '+tag)
    a=tone(x,1000,1)['peak_V'];b=tone(x,1000,2)['peak_V']
    t=x[-1,0]-np.arange(2000)/100000
    cv=np.interp(t,x[:,0],x[:,5]);cell=np.interp(t,x[:,0],x[:,7])
    row={'case':tag,'status':'OK','input_dBu':level,'pre_dBu':db(a/math.sqrt(2)/.775),'output_dBu':db(b/math.sqrt(2)/.775),'GR_dB':db(a/b),'mean_cell_GR_dB':float(np.mean(20*np.log10(1+10000*np.clip(cell,2e-8,.005)))),'timing_CV_mean_V':float(np.mean(cv)),'timing_CV_pp_V':float(np.ptp(cv)),**st}
    data[level]=row;save(row);return row
def find(target):
    for _ in range(16):
        pts=sorted(data.values(),key=lambda a:a['input_dBu'])
        closest=min(pts,key=lambda a:abs(a['GR_dB']-target))
        if abs(closest['GR_dB']-target)<.003:return closest
        pair=next(((a,b) for a,b in zip(pts,pts[1:]) if a['GR_dB']<=target<=b['GR_dB'] and b['GR_dB']>a['GR_dB']),None)
        if pair is None:return None
        a,b=pair
        lo,hi=a['input_dBu'],b['input_dBu']
        frac=(target-a['GR_dB'])/(b['GR_dB']-a['GR_dB'])
        sample(lo+(hi-lo)*min(.9,max(.1,frac)))
    return None
for level in [-42,-30,-18,-6,6,18]:sample(level)
baseline=sample(-42)['GR_dB']
cross={}
for target in [.05,.5,1,2,3,5,10,20,30]:
    center=find(target+baseline)
    if center is None:
        save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'result':'unreachable_or_not_converged'});continue
    cross[target]=center
    if target==.05:continue
    level=center['input_dBu']
    if level+.1>18:
        save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'result':'outside_input_range'});continue
    a,b=sample(level-.1),sample(level+.1)
    dx=b['pre_dBu']-a['pre_dBu'];dy=b['output_dBu']-a['output_dBu']
    slope=dy/dx
    save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'actual_GR_dB':center['GR_dB']-baseline,'center_input_dBu':level,'delta_input_dB':dx,'delta_output_dB':dy,'local_ratio':1/slope if abs(slope)>1e-9 else None,'slope':slope,'ratio_error_percent':100*(1/slope/ratio-1) if abs(slope)>1e-9 else None,'bracket':[a,b]})
# Extrapolate the established post-knee line from measured 3/5 dB crossings,
# independently of the near-zero onset proxy. Report both, do not hide curvature.
if 3 in cross and 5 in cross:
    a,b=cross[3],cross[5]
    slope=(b['GR_dB']-a['GR_dB'])/(b['pre_dBu']-a['pre_dBu'])
    intercept=a['pre_dBu']-(a['GR_dB']-baseline)/slope
    save({'case':'threshold_fit','status':'OK','threshold_pre_dBu':intercept,'fit_GR_slope':slope,'baseline_GR_dB':baseline})
for lo,hi in [(0,10),(10,20),(20,30)]:
    a,b=cross.get(.05 if lo==0 else lo),cross.get(hi)
    rr=None
    if a and b:
        dx=b['pre_dBu']-a['pre_dBu'];dy=b['output_dBu']-a['output_dBu']
        if abs(dy)>1e-9:rr=dx/dy
    save({'case':'segment','status':'OK','GR_interval_dB':[lo,hi],'onset_proxy_dB':.05 if lo==0 else None,'effective_ratio':rr})
