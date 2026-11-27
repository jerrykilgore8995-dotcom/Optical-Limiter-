"""Focused full-circuit tests, preserving vendor models and solver tolerances."""
from pathlib import Path
import sys, re, math, json
import numpy as np
# Reuse laboratory helpers without invoking the original CLI dispatcher.
p=Path(__file__).with_name('measure.py')
exec(compile(p.read_text().split("if SUITE=='audio' or")[0],str(p),'exec'))
def coupling():
    cap=float(SUITE.split('_')[1])
    def edit(s):
        s=s.replace('CVAL=1.005u',f'CVAL={cap*1.005}u').replace('CVAL=.995u',f'CVAL={cap*.995}u')
        s=s.replace('SIN(0 {AMP} {FREQ})','SIN(0 {AMP} {FREQ}) AC 1')
        # Zero volt series probes measure the actual current into each capacitor.
        s=s.replace('CIN in b {CVAL}','VSIN in capin 0\nCIN capin b {CVAL}')
        return s
    wd=OUT/'ac';wd.mkdir(exist_ok=True)
    shutil.copytree(ROOT.parent/'stage1'/'models',wd/'models',dirs_exist_ok=True)
    s=edit(BASE.replace('../stage1/models/','models/').replace('SENS=15','SENS=0'))
    s=re.sub(r'(?im)^\.(op|tran|end)\b[^\n]*\n?','',s)
    nodes=['ip','im','xp.b','xm.b','p','m','pre','gr','out']
    vec=['v('+n+')' for n in nodes]+['i(v.xp.vsin)','i(v.xm.vsin)']
    s+='\n.control\nset wr_singlescale\nset wr_vecnames\nac dec 60 2 40000\nwrdata wave.txt '+' '.join(vec)+'\nquit\n.endc\n.end\n'
    (wd/'case.cir').write_text(s)
    q=subprocess.run(['ngspice','-b','case.cir'],cwd=wd,capture_output=True,text=True,timeout=180)
    (wd/'ngspice.log').write_text(q.stdout+q.stderr)
    if q.returncode or re.search(r'error:|fatal',q.stdout+q.stderr,re.I):raise RuntimeError(q.stdout[-2000:]+q.stderr[-2000:])
    a=np.loadtxt(wd/'wave.txt',skiprows=1,ndmin=2)
    if a.shape[1]!=23 or not np.isfinite(a).all():raise RuntimeError(f'Unexpected AC data shape {a.shape}')
    z=a[:,1::2]+1j*a[:,2::2]
    f=a[:,0];mid=np.argmin(abs(f-1000))
    def polar(v):return {'mag':float(abs(v)),'phase_deg':float(np.angle(v,deg=True))}
    rows=[]
    for hz in [2,5,10,15,20,30,40,50,80,100,200,500,1000,5000,10000,20000,40000]:
        i=np.argmin(abs(f-hz));v=z[i];zin=v[0]/v[9];zload=v[2]/v[9]
        rows.append({'Hz':float(f[i]),'input_Z_ohm':polar(zin),'post_cap_Z_ohm':polar(zload),'cap_transfer_dB':db(v[2]/v[0]),'receiver_transfer_dB':db((v[4]-v[5])/(v[2]-v[3])),'pre_dB':db(v[6]),'gr_dB':db(v[7]),'out_dB':db(v[8]),'relative_out_dB':db(v[8]/z[mid,8]),'phase_deg':float(np.angle(v[8],deg=True))})
    save({'case':'ac','status':'OK','cap_uF':cap,'rows':rows})
    # Independent sine check (AC and transient should agree at small signal).
    for hz in [20,30,100,1000]:
        def case():
            x=run(f'tone_{hz}',{'SENS':0,'FREQ':hz,'AMP':peak(-18)},seconds=max(.12,12/hz),step=20e-6,edit=lambda s:s.replace('CVAL=1.005u',f'CVAL={cap*1.005}u').replace('CVAL=.995u',f'CVAL={cap*.995}u'))
            return {'Hz':hz,'cap_uF':cap,'transient_gain_dB':db(tone(x,hz,3)['peak_V']/peak(-18))}
        attempt(f'tone_{hz}',case)
def ratios():
    _,ff,soft,ratio=SUITE.split('_');ff,soft,ratio=map(int,(ff,soft,ratio))
    data={}
    def sample(level):
        level=round(level,4)
        if level in data:return data[level]
        tag='ratio_'+str(level)
        def case():
            x=run(tag,{'FF':ff,'SOFT':soft,'RATIO':ratio,'SENS':30,'AMP':peak(level)},seconds=.22,step=20e-6)
            st=steady(x)
            if not st['settled']:
                x=run(tag+'_long',{'FF':ff,'SOFT':soft,'RATIO':ratio,'SENS':30,'AMP':peak(level)},seconds=.45,step=20e-6)
                st=steady(x)
            a=tone(x,1000,1)['peak_V'];b=tone(x,1000,2)['peak_V']
            row={'input_dBu':level,'pre_dBu':db(a/math.sqrt(2)/.775),'output_dBu':db(b/math.sqrt(2)/.775),'GR_dB':db(a/max(b,1e-30)),**st}
            if st['settled']:data[level]=row
            return row
        attempt(tag,case)
        return data.get(level)
    for level in [-48,-36,-30,-24,-18,-12,-6,0,6,12,18]:sample(level)
    def bracket(g):
        a=sorted(data.values(),key=lambda x:x['input_dBu'])
        for p,q in zip(a,a[1:]):
            if p['GR_dB']<=g<=q['GR_dB'] and q['GR_dB']>p['GR_dB']:return p,q
        return None
    def crossing(g):
        for _ in range(4):
            pair=bracket(g)
            if pair is None:return None
            p,q=pair
            w=(g-p['GR_dB'])/(q['GR_dB']-p['GR_dB'])
            lv=p['input_dBu']+w*(q['input_dBu']-p['input_dBu'])
            if min(abs(p['GR_dB']-g),abs(q['GR_dB']-g))<.025:break
            if sample(lv) is None:break
        pair=bracket(g)
        if pair is None:return None
        p,q=pair;w=(g-p['GR_dB'])/(q['GR_dB']-p['GR_dB'])
        return {'pre_dBu':p['pre_dBu']+w*(q['pre_dBu']-p['pre_dBu']),'out_dBu':p['output_dBu']+w*(q['output_dBu']-p['output_dBu']),'bracket_GR_dB':[p['GR_dB'],q['GR_dB']]}
    c=[crossing(g) for g in [.05,10,20,30]]
    result=[]
    for i,(lo,hi) in enumerate([(0,10),(10,20),(20,30)]):
        a,b=c[i:i+2];ratio_eff=None
        if a and b and abs(b['out_dBu']-a['out_dBu'])>1e-8:ratio_eff=(b['pre_dBu']-a['pre_dBu'])/(b['out_dBu']-a['out_dBu'])
        result.append({'GR_interval_dB':[lo,hi],'onset_proxy_dB':.05 if lo==0 else None,'ratio':ratio_eff,'crossings':[a,b]})
    save({'case':'segments','status':'OK','selected_ratio':ratio,'FF':ff,'SOFT':soft,'segments':result})
if SUITE.startswith('cap_'):coupling()
elif SUITE.startswith('target_'):ratios()
else:raise SystemExit('Unknown diagnostic')
if any(r['status']!='OK' for r in ROWS):sys.exit(2)

