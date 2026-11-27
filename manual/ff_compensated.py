"""FF hard knee only. Full circuit, no calibration edits."""
from pathlib import Path
p=Path(__file__).with_name('measure.py')
exec(compile(p.read_text().split("if SUITE=='audio' or")[0],str(p),'exec'))
_,ratio_s,cap_s,audit=SUITE.split('_')
cap=int(cap_s)
assert audit in ('peakonly','all','baseline')
offset_s='0'
ratio=int(ratio_s)
offset=float(offset_s)
rf=10000*(1-1/ratio)+offset
variant='stock_peak'
assert ratio in (2,3,4,8,12,20)
def fixture(s):
    old='RRF rrout rrneg {10000*(FF ? 1-1/RATIO : RATIO-1)}'
    assert old in s
    # Retain original ratio and timing values.
    if cap:
        s=s.replace('RPKFB pkout pkneg 10k',f'RPKFB pkout pkneg 10k\nCPKFB pkout pkneg {cap}p')
    if audit=='all':
        for resistor,capname,nodes in [
            ('RSF','CSF','rect sa'),
            ('RAF','CAF','avgo an'),
            ('RSUB','CSUB','subout subneg'),
            ('RIF','CIF','ov100 invneg'),
            ('RRF','CRRF','rrout rrneg')]:
            s=re.sub(r'(?m)^'+resistor+r' .*
assert ratio in (2,3,4,8,12,20)
data={}
def sample(level):
    level=round(float(level),5)
    if level in data:return data[level]
    tag='level_'+str(level)
    params={'FF':1,'SOFT':0,'RATIO':ratio,'SENS':30,'AMP':peak(level)}
    for duration in [.18,.36,.72]:
        x=run(tag+'_'+str(duration),params,seconds=duration,step=10e-6,edit=fixture,signals=('pre','gr','out','rect','timing','sense','cell','det','logv','pbuf','cmd','hold','hfout','rrout'))
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
    row['variant']=audit
    row['peak_feedback_cap_pF']=cap
    row['ratio_feedback_ohm']=rf
    data[level]=row;save(row);return row
def find(target):
    for _ in range(12):
        pts=sorted(data.values(),key=lambda a:a['input_dBu'])
        closest=min(pts,key=lambda a:abs(a['GR_dB']-target))
        if abs(closest['GR_dB']-target)<.015:return closest
        pair=next(((a,b) for a,b in zip(pts,pts[1:]) if a['GR_dB']<=target<=b['GR_dB'] and b['GR_dB']>a['GR_dB']),None)
        if pair is None:return None
        a,b=pair
        lo,hi=a['input_dBu'],b['input_dBu']
        frac=(target-a['GR_dB'])/(b['GR_dB']-a['GR_dB'])
        sample(lo+(hi-lo)*min(.95,max(.05,frac)))
    return None
for level in [-36,-30,-12,6,18]:sample(level)
baseline=sample(-36)['GR_dB']
cross={}
ceiling=sample(18)['GR_dB']-baseline
save({'case':'maximum','status':'OK','ratio_selected':ratio,'feedback_ohm':rf,'peak_feedback_cap_pF':cap,'compensation':audit,'ratio_feedback_ohm':rf,'max_measured_GR_dB':ceiling,'input_dBu':18,'note':'Saturation ceiling is not a compression ratio'})
for target in [1,5,10,20,round(ceiling-0.65,3)]:
    center=find(target+baseline)
    if center is None:
        save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'result':'unreachable_or_not_converged'});continue
    cross[target]=center
    if target==.05:continue
    level=center['input_dBu']
    if level+.5>18:
        save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'result':'outside_input_range'});continue
    a,b=sample(level-.5),sample(level+.5)
    dx=b['pre_dBu']-a['pre_dBu'];dy=b['output_dBu']-a['output_dBu']
    slope=dy/dx
    save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'actual_GR_dB':center['GR_dB']-baseline,'center_input_dBu':level,'delta_input_dB':dx,'delta_output_dB':dy,'local_ratio':1/slope if abs(slope)>1e-9 else None,'slope':slope,'ratio_error_percent':100*(1/slope/ratio-1) if abs(slope)>1e-9 else None,'detector_slope_dB_per_dB':(db(b['det_mean_V'])-db(a['det_mean_V']))/dx,'command_slope_dB_per_dB':10*(b['timing_CV_mean_V']-a['timing_CV_mean_V'])/dx,'optical_slope_dB_per_dB':(b['mean_cell_GR_dB']-a['mean_cell_GR_dB'])/dx,'bracket':[a,b]})
# Extrapolate the established post-knee line from measured 3/5 dB crossings,
# independently of the near-zero onset proxy. Report both, do not hide curvature.
if 5 in cross and 10 in cross:
    a,b=cross[5],cross[10]
    slope=(b['GR_dB']-a['GR_dB'])/(b['pre_dBu']-a['pre_dBu'])
    intercept=a['pre_dBu']-(a['GR_dB']-baseline)/slope
    save({'case':'threshold_fit','status':'OK','threshold_pre_dBu':intercept,'fit_GR_slope':slope,'baseline_GR_dB':baseline})
,lambda m:m.group(0)+'\n'+capname+' '+nodes+' 22p',s)
    return s
assert ratio in (2,3,4,8,12,20)
data={}
def sample(level):
    level=round(float(level),5)
    if level in data:return data[level]
    tag='level_'+str(level)
    params={'FF':1,'SOFT':0,'RATIO':ratio,'SENS':30,'AMP':peak(level)}
    for duration in [.18,.36,.72]:
        x=run(tag+'_'+str(duration),params,seconds=duration,step=10e-6,edit=fixture,signals=('pre','gr','out','rect','timing','sense','cell','det','logv','pbuf','cmd','hold','hfout','rrout'))
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
    row['ratio_feedback_ohm']=rf
    data[level]=row;save(row);return row
def find(target):
    for _ in range(12):
        pts=sorted(data.values(),key=lambda a:a['input_dBu'])
        closest=min(pts,key=lambda a:abs(a['GR_dB']-target))
        if abs(closest['GR_dB']-target)<.015:return closest
        pair=next(((a,b) for a,b in zip(pts,pts[1:]) if a['GR_dB']<=target<=b['GR_dB'] and b['GR_dB']>a['GR_dB']),None)
        if pair is None:return None
        a,b=pair
        lo,hi=a['input_dBu'],b['input_dBu']
        frac=(target-a['GR_dB'])/(b['GR_dB']-a['GR_dB'])
        sample(lo+(hi-lo)*min(.95,max(.05,frac)))
    return None
for level in [-36,-30,-12,6,18]:sample(level)
baseline=sample(-36)['GR_dB']
cross={}
ceiling=sample(18)['GR_dB']-baseline
save({'case':'maximum','status':'OK','ratio_selected':ratio,'feedback_ohm':rf,'ratio_feedback_ohm':rf,'max_measured_GR_dB':ceiling,'input_dBu':18,'note':'Saturation ceiling is not a compression ratio'})
for target in [1,5,10,20,round(ceiling-0.65,3)]:
    center=find(target+baseline)
    if center is None:
        save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'result':'unreachable_or_not_converged'});continue
    cross[target]=center
    if target==.05:continue
    level=center['input_dBu']
    if level+.5>18:
        save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'result':'outside_input_range'});continue
    a,b=sample(level-.5),sample(level+.5)
    dx=b['pre_dBu']-a['pre_dBu'];dy=b['output_dBu']-a['output_dBu']
    slope=dy/dx
    save({'case':'local','status':'OK','ratio_selected':ratio,'target_GR_dB':target,'actual_GR_dB':center['GR_dB']-baseline,'center_input_dBu':level,'delta_input_dB':dx,'delta_output_dB':dy,'local_ratio':1/slope if abs(slope)>1e-9 else None,'slope':slope,'ratio_error_percent':100*(1/slope/ratio-1) if abs(slope)>1e-9 else None,'detector_slope_dB_per_dB':(db(b['det_mean_V'])-db(a['det_mean_V']))/dx,'command_slope_dB_per_dB':10*(b['timing_CV_mean_V']-a['timing_CV_mean_V'])/dx,'optical_slope_dB_per_dB':(b['mean_cell_GR_dB']-a['mean_cell_GR_dB'])/dx,'bracket':[a,b]})
# Extrapolate the established post-knee line from measured 3/5 dB crossings,
# independently of the near-zero onset proxy. Report both, do not hide curvature.
if 5 in cross and 10 in cross:
    a,b=cross[5],cross[10]
    slope=(b['GR_dB']-a['GR_dB'])/(b['pre_dBu']-a['pre_dBu'])
    intercept=a['pre_dBu']-(a['GR_dB']-baseline)/slope
    save({'case':'threshold_fit','status':'OK','threshold_pre_dBu':intercept,'fit_GR_slope':slope,'baseline_GR_dB':baseline})
