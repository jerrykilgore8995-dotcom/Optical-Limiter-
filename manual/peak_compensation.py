"""Fast detector-only tests; no production edits or timing changes."""
from pathlib import Path
p=Path(__file__).with_name('measure.py')
exec(compile(p.read_text().split("if SUITE=='audio' or")[0],str(p),'exec'))
_,rf_s,cf_s,diode,precision=SUITE.split('_')
rf=int(rf_s);cf=int(cf_s)
variant='direct'
assert diode in ('stock','lowvf')

# Remove everything except supply, rectifier (optional) and physical peak stage.
def fixture(s):
    head=s[:s.index('VCC vp')]
    rect=s[s.index('* Precision rectifier:'):s.index('* Average branch.')]
    pk=s[s.index('* Peak branch:'):s.index('* Provisional ideal log')]
    if variant.startswith('direct'):
        rect='BINPUT rect 0 V={abs(AMP*sin(2*pi*FREQ*time))}\n'
    else:
        rect='VIN hfout 0 SIN(0 {AMP} {FREQ})\n'+rect
    pk=pk.replace('RPKFB pkout pkneg 10k',f'RPKFB pkout pkneg {rf}')
    if cf:pk=pk.replace('XPK rect',f'CPKFB pkout pkneg {cf}p\nXPK rect')
    if diode=='lowvf':
        pk=pk.replace('DPC pkamp charge D4148W','DPC pkamp charge DTEST\n.model DTEST D(IS=1e-7 N=1.1 RS=1 CJO=2p)')
    if 'noclamp' in variant:
        pk=pk.replace('DCL pkneg pkamp D4148W','')
    if 'fb1k' in variant:
        pk=pk.replace('RPKFB pkout pkneg 10k','RPKFB pkout pkneg 1k')
    if 'fb100k' in variant:
        pk=pk.replace('RPKFB pkout pkneg 10k','RPKFB pkout pkneg 100k')
    if 'schottky' in variant:
        # Diagnostic generic diode only; not a vendor model or final choice.
        pk=pk.replace('DPC pkamp charge D4148W','DPC pkamp charge DTEST\n.model DTEST D(IS=1e-7 N=1.1 RS=1 CJO=2p)')
    if precision=='fine':
        head=head.replace('reltol=1e-4 abstol=1e-11','reltol=1e-6 abstol=1e-13 vntol=1e-8')
    return head+'VOP op 0 15\nVON on 0 -15\n'+rect+pk
prev=None
levels=[-36,-33,-30,-27,-24,-21,-18,-15,-12,-9,-6,-3,0,3,6]
for level in levels:
    try:
        x=run('level_'+str(level),{'AMP':peak(level),'FREQ':1000},seconds=.035,step=2e-6 if precision=='fine' else 10e-6,edit=fixture,signals=('rect','hold','pkout','pkamp','pkneg','charge'))
        t=x[-1,0]-.01+np.arange(4000)*.01/4000
        z=np.column_stack([np.interp(t,x[:,0],x[:,i]) for i in range(1,7)])
        v=float(z[:,2].mean())
        row={'case':'transfer','status':'OK','input_dBu':level,'peak_input_V':peak(level),'rect_peak_V':float(z[:,0].max()),'rect_mean_V':float(z[:,0].mean()),'hold_mean_V':v,'hold_pp_V':float(np.ptp(z[:,1])),'amp_min_V':float(z[:,3].min()),'amp_max_V':float(z[:,3].max()),'clamp_forward_fraction':float(np.mean(z[:,4]-z[:,3]>.35)),'detector_gain_dB':db(v/peak(level))}
        if prev:
            row['detector_slope']= (db(v)-db(prev['hold_mean_V']))/(level-prev['input_dBu'])
        row.update({'feedback_ohm':rf,'feedback_cap_pF':cf,'diode':diode})
        prev=row
        save(row)
        if level in (-30,-12,0):
            # Resolve high-frequency ringing separately over the last 1ms.
            # Start with a settled waveform and use an independent fine-step run.
            q=run('ring_'+str(level),{'AMP':peak(level),'FREQ':1000},seconds=.008,step=50e-9,edit=fixture,signals=('rect','hold','pkout','pkamp','pkneg','charge'))
            tt=q[-1,0]-.001+np.arange(20000)*50e-9
            yy=np.interp(tt,q[:,0],q[:,4])
            spec=np.abs(np.fft.rfft((yy-yy.mean())*np.hanning(len(yy))))
            hz=np.fft.rfftfreq(len(yy),50e-9)
            mask=hz>=20000
            idx=np.where(mask)[0][np.argmax(spec[mask])]
            save({'case':'ringing','status':'OK','input_dBu':level,'feedback_ohm':rf,'feedback_cap_pF':cf,'diode':diode,'dominant_above_20k_Hz':float(hz[idx]),'above_20k_rms_V':float(np.sqrt(2*np.sum(spec[mask]**2)/(len(yy)*np.sum(np.hanning(len(yy))**2)))),'amp_pp_V':float(np.ptp(yy)),'note':'HF energy includes switching edges; inspect waveform before claiming oscillation'})

    except Exception as e:save({'case':'transfer','status':'FAILED','input_dBu':level,'error':str(e)})
if any(x['status']!='OK' for x in ROWS):sys.exit(2)
