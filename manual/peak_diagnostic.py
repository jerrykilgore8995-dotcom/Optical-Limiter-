"""FF hard knee only. Full circuit, no calibration edits."""
from pathlib import Path
p=Path(__file__).with_name('measure.py')
exec(compile(p.read_text().split("if SUITE=='audio' or")[0],str(p),'exec'))
variant=SUITE.split('_')[1]
ratio=int(SUITE.split('_')[-1])
assert variant in ('stock','r22','r10','reference')
def fixture(s):
    s=s.replace('.option reltol=1e-4 abstol=1e-11 method=gear','.option reltol=1e-6 abstol=1e-13 vntol=1e-8 method=gear')
    if variant in ('r22','r10'):s=s.replace('RPC charge hold 100','RPC charge hold '+variant[1:])
    if variant=='reference':
        # Diagnostic only: scale-invariant active peak hold with same nominal RC.
        start=s.index('RPKFB pkout');end=s.index('* Provisional ideal log',start)
        s=s[:start]+'''CP hold 0 100n\nRP hold 0 47k\nBCH 0 hold I={max(0,V(rect,hold))/10}\nXPKBUF hold pkout op on pkout OPA134\n'''+s[end:]
    return s
assert ratio in (2,3,4,8,12,20)
data={}
def sample(level):
    level=round(float(level),5)
    if level in data:return data[level]
    tag='level_'+str(level)
    params={'FF':1,'SOFT':0,'RATIO':ratio,'SENS':30,'AMP':peak(level)}
    for duration in [.18,.36,.72]:
        x=run(tag+'_'+str(duration),params,seconds=duration,step=2e-6,edit=fixture,signals=('pre','gr','out','rect','timing','sense','cell','det','logv','pbuf','cmd','hold','hfout','rrout'))
        st=steady(x)
        if st['settled']:break
    if not st['settled']:raise RuntimeError('Unsettled '+tag)
    a=tone(x,1000,1)['peak_V'];b=tone(x,1000,2)['peak_V']
    t=x[-1,0]-np.arange(2000)/100000
    cv=np.interp(t,x[:,0],x[:,5]);cell=np.interp(t,x[:,0],x[:,7])
    row={'case':tag,'status':'OK','input_dBu':level,'pre_dBu':db(a/math.sqrt(2)/.775),'output_dBu':db(b/math.sqrt(2)/.775),'GR_dB':db(a/b),'mean_cell_GR_dB':float(np.mean(20*np.log10(1+10000*np.clip(cell,2e-8,.005)))),'timing_CV_mean_V':float(np.mean(cv)),'timing_CV_pp_V':float(np.ptp(cv)),**st}
    for name,col in [('rect',4),('det',8),('hold',12),('hfout',13),('rrout',14),('logv',9),('cmd',11)]:
        v=np.interp(t,x[:,0],x[:,col])
        row[name+'_mean_V']=float(np.mean(v))
        row[name+'_peak_V']=float(np.max(v))
        row[name+'_pp_V']=float(np.ptp(v))
    row['detector_gain_dB']=db(row['det_mean_V']/a)
    row['variant']=variant
    data[level]=row;save(row);return row
def find(target):
    for _ in range(24):
        pts=sorted(data.values(),key=lambda a:a['input_dBu'])
        closest=min(pts,key=lambda a:abs(a['GR_dB']-target))
        if abs(closest['GR_dB']-target)<.003:return closest
        pair=next(((a,b) for a,b in zip(pts,pts[1:]) if a['GR_dB']<=target<=b['GR_dB'] and b['GR_dB']>a['GR_dB']),None)
        if pair is None:return None
        a,b=pair
        lo,hi=a['input_dBu'],b['input_dBu']
        frac=(target-a['GR_dB'])/(b['GR_dB']-a['GR_dB'])
        sample(lo+(hi-lo)*min(.8,max(.2,frac)))
    return None
for level in [-36,-30,-18,-6,6,18]:sample(level)
baseline=sample(-36)['GR_dB']
cross={}
for target in [.05,.5,1,2,5,10,20,30]:
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
    save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'actual_GR_dB':center['GR_dB']-baseline,'center_input_dBu':level,'delta_input_dB':dx,'delta_output_dB':dy,'local_ratio':1/slope if abs(slope)>1e-9 else None,'slope':slope,'ratio_error_percent':100*(1/slope/ratio-1) if abs(slope)>1e-9 else None,'detector_slope_dB_per_dB':(db(b['det_mean_V'])-db(a['det_mean_V']))/dx,'command_slope_dB_per_dB':10*(b['timing_CV_mean_V']-a['timing_CV_mean_V'])/dx,'optical_slope_dB_per_dB':(b['mean_cell_GR_dB']-a['mean_cell_GR_dB'])/dx,'bracket':[a,b]})
# Extrapolate the established post-knee line from measured 3/5 dB crossings,
# independently of the near-zero onset proxy. Report both, do not hide curvature.
if 2 in cross and 5 in cross:
    a,b=cross[2],cross[5]
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
