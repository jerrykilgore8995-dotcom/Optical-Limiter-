"""Flatten the locked SPICE topology into a device-level native-engine IR.
Preserve all model parameters and behavioral expressions; never silently drop a device.
"""
import json,re,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sub={}; models={}; top=[]; functions=[]; options=[]
def read(path, title=False):
    lines=[]
    for raw in path.read_text(encoding='latin1').splitlines()[int(title):]:
        line=raw.strip()
        if not line or line.startswith('*'):continue
        if line.startswith('+'):lines[-1]+=' '+line[1:].strip()
        else:lines.append(line)
    current=top
    for line in lines:
        t=line.split(); cmd=t[0].lower()
        if cmd=='.include':read((path.parent/t[1]).resolve())
        elif cmd=='.subckt':
            pos=next((i for i,v in enumerate(t) if v.lower()=='params:'),len(t))
            body=[];sub[t[1].lower()]={'pins':t[2:pos],'defaults':t[pos+1:],'body':body};current=body
        elif cmd=='.ends':current=top
        elif cmd=='.func':functions.append(line)
        elif cmd in ('.op','.tran','.end'):continue
        elif cmd in ('.option','.options'):options.append(line)
        else:current.append(line)
read(root/'manual/manual.cir',True)
devices=[];parameters={}
def expand(lines,prefix='',mapping=None,scope=None):
    mapping=mapping or {};scope=dict(scope or {})
    def node(n):
        n=n.lower()
        return '0' if n=='0' else mapping.get(n,prefix+n)
    def branch(n):return prefix+n.lower()
    localmodels={}
    for line in lines:
        t=line.split();cmd=t[0].lower()
        if cmd=='.param':
            for p in t[1:]:
                k,v=p.split('=',1);scope[k.lower()]=v
                if not prefix:parameters[k.lower()]=v
        elif cmd=='.model':
            name=prefix+t[1].lower();localmodels[t[1].lower()]=name
            models[name]=' '.join(t[2:])
    for line in lines:
        t=line.split();cmd=t[0].lower();kind=cmd[0]
        if cmd in ('.param','.model'):continue
        if kind=='x':
            j=next((i for i,v in enumerate(t) if '=' in v),len(t))
            model=sub[t[j-1].lower()];args=t[1:j-1]
            assert len(args)==len(model['pins']),line
            env=dict(scope)
            for p in model['defaults']+t[j:]:
                k,v=p.split('=',1);env[k.lower()]=v
            expand(model['body'],branch(cmd)+'.',dict(zip(map(str.lower,model['pins']),map(node,args))),env)
            continue
        assert kind in 'rclvidqjegfhbs',line
        count=3 if kind in 'qj' else 4 if kind=='s' else 2
        pins=list(map(node,t[1:1+count]))
        tail=' '.join(t[1+count:])
        rec={'name':branch(cmd),'kind':kind,'nodes':pins,'raw_value':tail,'scope':scope,'source':line}
        if kind in 'dqjs':
            key=t[1+count].lower()
            rec['model']=localmodels.get(key,key)
        if kind in 'egfh':
            # preserve voltage/current control references explicitly, including TI POLY sources
            a=t[1+count:]
            if a[0].lower().startswith('poly('):
                n=int(a[0][5:-1])
                if kind in 'eg':
                    pairs=re.findall(r'\(\s*([^,]+),\s*([^)]+)\)', ' '.join(a[1:]))
                    assert len(pairs)==n,line
                    rec['control_nodes']=[[node(u),node(v)] for u,v in pairs]
                    rec['coefficients']=re.sub(r'\([^)]*\)','', ' '.join(a[1:])).split()
                else:
                    rec['control_branches']=list(map(branch,a[1:1+n]));rec['coefficients']=a[1+n:]
            elif kind in 'eg':
                rec['control_nodes']=[[node(a[0]),node(a[1])]];rec['coefficients']=['0',a[2]]
            else:
                rec['control_branches']=[branch(a[0])];rec['coefficients']=['0',a[1]]
        if kind=='b':
            rec['expression']=re.sub(r'(?i)\bv\(([^)]+)\)',lambda m:'v('+','.join(node(n.strip()) for n in m[1].split(','))+')',tail)
        devices.append(rec)
expand(top)
for d in devices:
    if 'model' in d:assert d['model'] in models,d
from collections import Counter
out={'source':'manual/manual.cir','parameters':parameters,'functions':functions,'options':options,'models':models,'devices':devices,'counts':dict(Counter(d['kind'] for d in devices))}
dest=Path(sys.argv[1]) if len(sys.argv)>1 else root/'native/deck.json'
dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2))
print('Flattened device counts:',out['counts'],'models:',len(models))
