"""Generate full native DC/transient stamping code from the locked SPICE IR."""
import json,re,math
from pathlib import Path
root=Path(__file__).resolve().parent
d=json.loads((root/'deck.json').read_text())
dev=d['devices']
nodes=sorted({n for x in dev for n in x['nodes'] if n!='0'})
# Explicit internal terminals for semiconductor series resistances.
scale={'':1,'t':1e12,'g':1e9,'meg':1e6,'k':1e3,'m':1e-3,'u':1e-6,'n':1e-9,'p':1e-12,'f':1e-15}
def literal(m):return repr(float(m[1])*scale[m[2].lower()])
def expr(s):
    s=s.strip().lower().strip('{}')
    s=re.sub(r'(?<![\w.])([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?)(meg|[tgkmunpf])(?![\w.])',lambda m:literal(m),s)
    # SPICE conditional expressions, including nested ternary inside parentheses.
    def tern(t):
        depth=0;q=-1
        for i,c in enumerate(t):
            if c=='(':depth+=1
            elif c==')':depth-=1
            elif c=='?' and depth==0:q=i;break
        if q>=0:
            depth=0;nest=0
            for i in range(q+1,len(t)):
                c=t[i]
                if c=='(':depth+=1
                elif c==')':depth-=1
                elif depth==0:
                    if c=='?':nest+=1
                    elif c==':':
                        if nest==0:return '('+tern(t[:q])+'?'+tern(t[q+1:i])+':'+tern(t[i+1:])+')'
                        nest-=1
            raise ValueError(t)
        out='';i=0
        while i<len(t):
            if t[i]=='(':
                j=i+1;k=1
                while k:
                    if t[j]=='(':k+=1
                    elif t[j]==')':k-=1
                    j+=1
                out+='('+tern(t[i+1:j-1])+')';i=j
            else:out+=t[i];i+=1
        return out
    return tern(s)
def model(name):
    s=d['models'][name].lower()
    typ=s.split('(')[0].strip()
    vals={k:float(expr(v)) for k,v in re.findall(r'(\w+)\s*=\s*([^\s,)]+)',s)}
    return typ,vals
for x in dev:
    if x['kind'] in 'dqj':
        t,p=model(x['model'])
        terminals={'d':[(0,'rs')],'q':[(0,'rc'),(1,'rb'),(2,'re')],'j':[(0,'rd'),(2,'rs')]}[x['kind']]
        for i,k in terminals:
            if p.get(k,0)>0:
                internal=x['name']+'#'+k
                nodes.append(internal)
                x.setdefault('series',[]).append((x['nodes'][i],internal,p[k]))
                x['nodes'][i]=internal
nid={n:i for i,n in enumerate(nodes)}
branches={x['name']:len(nodes)+i for i,x in enumerate(x for x in dev if x['kind'] in 'veh' or x['kind']=='b' and x['raw_value'].lower().startswith('v='))}
def ix(n):return -1 if n=='0' else nid[n]
def val(x,s):
    s=s.strip()
    if s.startswith('{') and s.endswith('}'):
        inner=s[1:-1].lower()
        if inner in x['scope'] and inner not in d['parameters']:s=x['scope'][inner]
    return expr(s)
def num(x):return format(x,'.17g')
def params(p,keys,defaults):return ','.join(num(p.get(k,v)) for k,v in zip(keys.split(),defaults))
out=['#pragma once','#include "MNA.hpp"','#include <cmath>','#include <algorithm>','#include <vector>','namespace optical {',
'struct Settings { double gdb=0,freq=1000,amp=1,ff=1,ratio=4,soft=1,sens=15,hpf=0,att=.01,rel=.1,makeup=0; };',
'inline double min(double a,double b){return std::min(a,b);} inline double max(double a,double b){return std::max(a,b);}',
'inline double cg(double x){return min(.005,max(2e-8,x));} inline double cs(double x){return min(4,max(0,x));}',
'inline double ternary_fcn(bool a,double b,double c){return a?b:c;}',
f'inline constexpr int nodeCount={len(nodes)}, unknownCount={len(nodes)+len(branches)};',
f'inline constexpr int outputNode={ix("out")}, preNode={ix("pre")}, grNode={ix("gr")}, cellNode={ix("cell")};',
'inline void stampDeck(MNA& m,const Settings& s,double input,double dt,const std::vector<double>& previous){',
'using std::exp;using std::log;using std::log10;using std::pow;using std::sqrt;',
'const auto [gdb,freq,amp,ff,ratio,soft,sens,hpf,att,rel,makeup]=s;',
'auto v=[&](int n){return m.voltage(n);};auto old=[&](int n){return n<0?0.:previous[n];};']
for x in dev:
    k=x['kind'];a,b,*rest=map(ix,x['nodes']);name=x['name'];raw=x['raw_value'];line=[]
    for p,n,r in x.get('series',[]):line.append(f'm.conductance({ix(p)},{ix(n)},{num(1/r)});')
    if k=='r':line.append(f'm.conductance({a},{b},1./({val(x,raw)}));')
    elif k=='c':line.append(f'if(dt>0)m.capacitor({a},{b},{val(x,raw.split()[0])},old({a})-old({b}),dt);')
    elif k=='v':
        v='input' if raw.lower().startswith('sin(') else val(x,re.sub(r'(?i)^dc\s+','',raw))
        line.append(f'm.source({a},{b},{branches[name]},{v});')
    elif k=='i':line.append(f'm.residual({a},{b},{val(x,re.sub(r"(?i)^dc\s+","",raw))});')
    elif k in 'egfh':
        c=x['coefficients'];assert len(c)>=2,x
        controls=[f'v({ix(p)})-v({ix(n)})' for p,n in x.get('control_nodes',[])] if k in 'eg' else [f'v({branches[n]})' for n in x['control_branches']]
        assert len(c)==len(controls)+1,x
        if k in 'eh':line.append(f'm.source({a},{b},{branches[name]},{val(x,c[0])});')
        else:line.append(f'm.residual({a},{b},{val(x,c[0])});')
        for i,g in enumerate(c[1:]):
            if k in 'eg':
                p,n=map(ix,x['control_nodes'][i])
                if k=='e':
                    line.append(f'm.residual({branches[name]},-1,-({val(x,g)})*({controls[i]}));m.derivative({branches[name]},-1,{p},{n},-({val(x,g)}));')
                else:line.append(f'm.controlled({a},{b},{p},{n},{val(x,g)});')
            else:
                p=branches[x['control_branches'][i]]
                if k=='h':line.append(f'm.residual({branches[name]},-1,-({val(x,g)})*v({p}));m.derivative({branches[name]},-1,{p},-1,-({val(x,g)}));')
                else:line.append(f'm.cccs({a},{b},{p},{val(x,g)});')
    elif k=='d':
        t,p=model(x['model']);line.append(f'm.diode({a},{b},{params(p,"is n",[1e-14,1])});')
    elif k=='q':
        t,p=model(x['model']);line.append(f'm.bjt({a},{b},{rest[0]},{-1 if t=="pnp" else 1},{{{params(p,"is bf nf br nr vaf var ikf ikr ise ne isc nc",[1e-16,100,1,1,1,0,0,0,0,0,1.5,0,2])}}});')
    elif k=='j':
        t,p=model(x['model']);line.append(f'm.jfet({a},{b},{rest[0]},{-1 if t=="pjf" else 1},{params(p,"vto beta is lambda",[-2,1e-4,1e-14,0])});')
    elif k=='s':
        t,p=model(x['model']);line.append(f'm.conductance({a},{b},v({rest[0]})-v({rest[1]})>{num(p.get("vt",0))}?{num(1/p.get("ron",1))}:{num(1/p.get("roff",1e12))});')
    elif k=='b':
        mode,body=raw.split('=',1)
        refs=[]
        def replace(m):
            ns=m[1].lower().split(',');ids=[ix(n.strip()) for n in ns]
            for n in ids:
                if n>=0 and n not in refs:refs.append(n)
            return '('+'-'.join(f'w({n})' for n in ids)+')'
        body=re.sub(r'(?i)\bv\(([^)]+)\)',replace,body)
        body=expr(body)
        key=branches[name] if mode.lower()=='v' else -1
        line.append(f'm.behavior({a},{b},{key},{{{",".join(map(str,refs))}}},[&](auto w){{return {body};}});')
    else:raise ValueError(x)
    out.append(' '.join(line)+' // '+name)
out+=['}','}']
(root/'FullDeck.hpp').write_text('\n'.join(out)+'\n')
print('Generated full circuit stamps:',len(dev),'unknowns:',len(nodes)+len(branches))
