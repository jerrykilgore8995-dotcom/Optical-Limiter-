"""Measurement-only laboratory. No MNA, no Auto. All waveforms from ngspice."""
from pathlib import Path
import sys, subprocess, shutil, re, json, math
import numpy as np
ROOT=Path(__file__).resolve().parent
SUITE=sys.argv[1]
OUT=ROOT/'measurements'/SUITE
OUT.mkdir(parents=True,exist_ok=True)
BASE=(ROOT/'manual.cir').read_text()
ROWS=[]
def save(row):
    ROWS.append(row)
    (OUT/'results.json').write_text(json.dumps(ROWS,indent=2,allow_nan=False))
    print(json.dumps(row,allow_nan=False),flush=True)
def peak(db): return math.sqrt(2)*.775*10**(db/20)
def db(x): return 20*math.log10(max(abs(float(x)),1e-30))
def run(tag,params=None,seconds=.25,step=5e-6,edit=lambda x:x,signals=('pre','gr','out','rect','timing','sense','cell','det','logv','pbuf','cmd'),ac=False):
    wd=OUT/tag; wd.mkdir(exist_ok=True)
    shutil.copytree(ROOT.parent/'stage1'/'models',wd/'models',dirs_exist_ok=True)
    deck=BASE.replace('../stage1/models/','models/')
    for k,v in (params or {}).items():
        deck=re.sub(r'(?<![A-Za-z0-9_])'+k+r'=[^\s]+',f'{k}={v}',deck,count=1)
    deck=edit(deck)
    deck=re.sub(r'(?im)^\.(op|tran|end)\b[^\n]*\n?','',deck)
    cmd='ac dec 48 10 40000\nwrdata wave.txt v(out)' if ac else f'tran {step} {seconds} 0 {step}\nwrdata wave.txt '+' '.join('v('+s+')' for s in signals)
    deck+='\n.control\nset wr_singlescale\nset wr_vecnames\n'+cmd+'\nquit\n.endc\n.end\n'
    (wd/'case.cir').write_text(deck)
    p=subprocess.run(['ngspice','-b','case.cir'],cwd=wd,capture_output=True,text=True,timeout=900)
    log=p.stdout+p.stderr;(wd/'ngspice.log').write_text(log)
    if p.returncode or re.search(r'error:|timestep too small|simulation aborted|fatal',log,re.I): raise RuntimeError(log[-2500:])
    x=np.loadtxt(wd/'wave.txt',skiprows=1,ndmin=2)
    if not np.isfinite(x).all() or len(x)<3 or (not ac and x[-1,0]<seconds*.99999): raise RuntimeError('Incomplete or nonfinite waveform')
    return x
def tone(x,f,col):
    cycles=max(3,min(20,int(x[-1,0]*f/3)))
    rate=max(96000,64*f); n=int(math.ceil(rate*cycles/f))
    t=x[-1,0]-cycles/f+np.arange(n)*cycles/f/n
    y=np.interp(t,x[:,0],x[:,col]);y-=y.mean()
    a=[abs(2/n*np.dot(y,np.exp(-2j*np.pi*k*f*t))) for k in range(1,21)]
    return {'peak_V':float(a[0]),'THD_percent':float(100*np.linalg.norm(a[1:])/max(a[0],1e-30)),'harmonics_normalized_dBFS':{f'K{k}':db(a[k-1]/max(a[0],1e-30)) if k*f<22050 else None for k in range(1,21)}}
def attempt(tag,fn):
    try: save({'case':tag,'status':'OK',**fn()})
    except Exception as e: save({'case':tag,'status':'FAILED','error':str(e)})
def audio():
    def response(load):
        x=run(f'ac_{load}',{'SENS':0},edit=lambda s:s.replace('RL out 0 10k',f'RL out 0 {load}').replace('SIN(0 {AMP} {FREQ})','SIN(0 {AMP} {FREQ}) AC 1'),ac=True)
        z=x[:,1]+1j*x[:,2]; mid=np.argmin(abs(x[:,0]-1000))
        return {'load_ohm':load,'frequency_Hz':x[:,0].tolist(),'gain_dB':[db(a) for a in z],'relative_dB':[db(a/z[mid]) for a in z],'phase_deg':np.degrees(np.unwrap(np.angle(z))).tolist()}
    for l in [600,2000,10000,47000]: attempt(f'ac_{l}',lambda l=l:response(l))
    for f in [100,1000,5000,10000]:
        for level in [-18,-12,-8,-4,0,4,8,12,18]:
            tag=f'tone_{f}_{level}'
            def case():
                x=run(tag,{'SENS':0,'FREQ':f,'AMP':peak(level)},seconds=max(.15,25/f),step=min(5e-6,1/(64*f)))
                return {'Hz':f,'input_dBu':level,'receiver':tone(x,f,1),'optical_amp':tone(x,f,2),'output':tone(x,f,3)}
            attempt(tag,case)
def optical():
    # Isolated conductance law driven by specified sense voltage; tests model
    # implementation, not independent cell data. Remove actual driver only in fixture.
    def forced(s,ma):
        return re.sub(r'(?m)^XDRV .*\nRB drv base .*\nQD .*\nRS .*\nDLED .*\nDREV .*\n(?:\*[^\n]*\n)?DBE .*\nRBE .*\n',f'VSENSE sense 0 {ma/10}\n',s)
    for ma in [0,.1,.5,1,2,5,10,20,30,40]:
        tag=f'cell_{ma}'
        def case():
            x=run(tag,{'SENS':0,'AMP':.1},edit=lambda s:forced(s,ma),seconds=.08,signals=('sense','cell','gr','pre'))
            g=max(x[-1,2],2e-8)
            return {'forced_equivalent_mA':ma,'resistance_ohm':float(1/g),'static_GR_dB':db(1+10000*g),'model_only':True}
        attempt(tag,case)
    def recovery():
        x=run('recovery',{'SENS':0,'AMP':.1},edit=lambda s:forced(s,40).replace('VSENSE sense 0 4.0','VSENSE sense 0 PULSE(4 0 20m 1u 1u 100m 1)'),seconds=.09,signals=('sense','cell'))
        ix=np.where((x[:,0]>=.020001)&(x[:,2]<=1e-5))[0]
        return {'recovery_to_100k_ms':float((x[ix[0],0]-.020001)*1000) if len(ix) else None,'model_only':True,'voltage_dependent_LDR_distortion':'NOT IMPLEMENTED'}
    attempt('recovery',recovery)
def threshold():
    for sens in range(0,31,5):
        for level in [-30,-12,0,12]:
            tag=f'threshold_{sens}_{level}'
            def case():
                x=run(tag,{'SENS':sens,'AMP':peak(level)},seconds=.3)
                t=x[-1]; return {'sensitivity_dB':sens,'input_dBu':level,'pot_V':float(t[10]),'detector_V':float(t[8]),'log_V':float(t[9]),'expected_log_V':2*math.log10(max(1e-9,min(15,t[8]))),'command_V':float(t[11]),'GR_dB':db(tone(x,1000,1)['peak_V']/max(tone(x,1000,2)['peak_V'],1e-30))}
            attempt(tag,case)
def ratio_suite():
    _,_,ff,soft,ratio=SUITE.split('_');ff,soft,ratio=map(int,(ff,soft,ratio))
    good=[]
    for level in range(-48,19,2):
        tag=f'ratio_{level}'
        def case():
            x=run(tag,{'FF':ff,'SOFT':soft,'RATIO':ratio,'SENS':30,'AMP':peak(level)},seconds=.45)
            a=tone(x,1000,1)['peak_V'];b=tone(x,1000,2)['peak_V']
            row={'input_dBu':level,'pre_dBu':db(a/math.sqrt(2)/.775),'output_dBu':db(b/math.sqrt(2)/.775),'GR_dB':db(a/max(b,1e-30))}
            good.append(row);return row
        attempt(tag,case)
    # Interpolate crossing along increasing input, do not extrapolate.
    def cross(target):
        for a,b in zip(good,good[1:]):
            if b['input_dBu']-a['input_dBu']!=2:continue
            ga,gb=a['GR_dB'],b['GR_dB']
            if ga<=target<=gb and gb>ga:
                w=(target-ga)/(gb-ga)
                return tuple(a[k]+w*(b[k]-a[k]) for k in ('pre_dBu','output_dBu'))
        return None
    # 0 GR is asymptotic; use low-level baseline plus 0.05 dB onset proxy.
    segments=[]
    for lo,hi in [(0,10),(10,20),(20,30)]:
        a,b=cross(.05 if lo==0 else lo),cross(hi)
        segments.append({'GR_interval_dB':[lo,hi],'onset_proxy_dB':.05 if lo==0 else None,'effective_ratio':(b[0]-a[0])/(b[1]-a[1]) if a and b and abs(b[1]-a[1])>1e-9 else None,'status':'measured_crossings' if a and b else 'unreachable_or_missing'})
    (OUT/'segments.json').write_text(json.dumps(segments,indent=2))
def timing():
    # Separate independently dispatched timing suite per knob pair.
    for att,rel in ([(int(SUITE.split('_')[1]),int(SUITE.split('_')[2]))] if SUITE.startswith('timing_') else [(1,10),(10,100),(40,3000),(1,3000),(40,10)]):
        for level in [-12,-4,4,12,18]:
            tag=f'burst_{att}_{rel}_{level}'
            def case():
                rise=.15;fall=rise+max(.3,att/1000*8);end=fall+max(.2,rel/1000*5)
                def edit(s):
                    return re.sub(r'(?m)^VIN src.*$',f'BINPUT src 0 V={{sin(2*pi*1000*time)*({peak(level)})*ternary_fcn(time<{rise},0,ternary_fcn(time<{fall},1,0))}}',s)
                x=run(tag,{'ATT':att/1000,'REL':rel/1000,'SENS':30},seconds=end,step=20e-6,edit=edit,signals=('pre','gr','cell','timing'))
                t=x[:,0];g=20*np.log10(1+10000*np.clip(x[:,3],2e-8,.005))
                plateau=float(np.median(g[(t>fall-.05)&(t<fall)]));base=float(np.median(g[t<rise-.01]))
                ai=np.where((t>=rise)&(t<fall)&(g>=base+.632*(plateau-base)))[0];ri=np.where((t>=fall)&(g<=base+.368*(plateau-base)))[0]
                return {'attack_nominal_ms':att,'release_nominal_ms':rel,'input_dBu':level,'settled_GR_dB':plateau,'attack_63pct_ms':float((t[ai[0]]-rise)*1000) if len(ai) and plateau-base>.1 else None,'release_63pct_ms':float((t[ri[0]]-fall)*1000) if len(ri) and plateau-base>.1 else None,'GR_source':'optical conductance, not audio envelope'}
            attempt(tag,case)
def ripple():
    for f in [20,50,100]:
        for ff in [0,1]:
            for ratio in [4,12]:
                for level in [-4,4,12]:
                    tag=f'ripple_{f}_{ff}_{ratio}_{level}'
                    def case():
                        x=run(tag,{'FREQ':f,'FF':ff,'RATIO':ratio,'SENS':30,'AMP':peak(level)},seconds=1,step=10e-6)
                        a=x[x[:,0]>.8];g=20*np.log10(1+10000*np.clip(a[:,7],2e-8,.005))
                        return {'Hz':f,'FF':ff,'ratio':ratio,'input_dBu':level,'mean_GR_dB':float(g.mean()),'GR_ripple_pp_dB':float(np.ptp(g)),'timing_CV_pp_V':float(np.ptp(a[:,5])),'LED_sense_pp_V':float(np.ptp(a[:,6])),'output':tone(x,f,3)}
                    attempt(tag,case)
if SUITE=='audio':audio()
elif SUITE=='optical':optical()
elif SUITE=='threshold':threshold()
elif SUITE.startswith('ratio_case_'):ratio_suite()
elif SUITE=='timing' or SUITE.startswith('timing_'):timing()
elif SUITE=='ripple':ripple()
else:raise SystemExit('Unknown suite')
if any(r['status']!='OK' for r in ROWS):sys.exit(2)
