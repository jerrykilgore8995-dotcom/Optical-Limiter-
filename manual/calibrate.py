from pathlib import Path
import subprocess, re, json, math
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'; OUT.mkdir(exist_ok=True)
def run(name, deck):
    p=OUT/(name+'.cir'); p.write_text(deck)
    # includes resolve relative to manual working directory
    proc=subprocess.run(['ngspice','-b',str(p)],cwd=ROOT,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=180)
    (OUT/(name+'.log')).write_text(proc.stdout)
    if proc.returncode or re.search(r'(?im)^\s*(error:|fatal|.*timestep too small)',proc.stdout):
        raise RuntimeError(name+' failed; see log')
    return proc.stdout
def bias():
    deck=(ROOT/'receiver_bias.cir').read_text()
    results=[]
    for reval in [22,33,47,68,100]:
        name=f'bias_{reval}'
        log=run(name,deck.replace('REB=47',f'REB={reval}'))
        values=[float(x) for x in re.findall(r'i\(v\.x[pm]\.vq[np]\)\s*=\s*([-+\d.eE]+)',log,re.I)]
        if len(values)!=4: raise RuntimeError('missing bias currents: '+name)
        results.append({'REB':reval,'collector_mA':[abs(x)*1000 for x in values],'error':sum((abs(x)*1000-3)**2 for x in values)})
    (OUT/'bias.json').write_text(json.dumps(sorted(results,key=lambda x:x['error']),indent=2))
def ratios():
    source=(ROOT/'manual.cir').read_text()
    rows=[]
    for ff in [0,1]:
        for soft in [0,1]:
            for ratio in [2,3,4,8,12,20]:
                for db in range(-30,19,3):
                    amp=math.sqrt(2)*.775*10**(db/20)
                    deck=re.sub(r'(?m)^\.param GDB=.*$',f'.param GDB=0 FREQ=1000 AMP={amp} FF={ff} RATIO={ratio} SOFT={soft} SENS=30 HPF=0 ATT=10m REL=100m MAKEUP=0',source)
                    deck=deck.replace('.tran 0.5u 200m','.tran 2u 400m\n.measure tran pre_rms RMS v(pre) FROM=350m TO=400m\n.measure tran gr_rms RMS v(gr) FROM=350m TO=400m')
                    name=f'r_{ff}_{soft}_{ratio}_{db}'
                    log=run(name,deck)
                    vals=[re.search(k+r'\s*=\s*([-+\d.eE]+)',log) for k in ['pre_rms','gr_rms']]
                    if not all(vals): raise RuntimeError('missing measurement '+name)
                    a,b=[float(x.group(1)) for x in vals]
                    rows.append(dict(ff=ff,soft=soft,ratio=ratio,input_dbu=db,pre_rms=a,gr_rms=b,gr_db=20*math.log10(a/b)))
                    (OUT/'ratio_raw.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':
    import sys
    bias()
    if '--ratios' in sys.argv: ratios()
