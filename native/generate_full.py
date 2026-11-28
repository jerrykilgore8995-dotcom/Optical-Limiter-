"""Generate full native DC/transient stamping code from the locked SPICE IR."""
import json,re,math,sys
from pathlib import Path
root=Path(__file__).resolve().parent
deck_path=Path(sys.argv[1]) if len(sys.argv)>1 else root/'deck.json'
seed_path=Path(sys.argv[2]) if len(sys.argv)>2 else None
output_path=Path(sys.argv[3]) if len(sys.argv)>3 else root/'FullDeck.hpp'
d=json.loads(deck_path.read_text())
dev=d['devices']
condense = '--condense' in sys.argv[4:]
sidechain = {}
if condense:
    start = next(i for i, element in enumerate(dev) if element['name'] == 'bselect')
    sidechain = {element['name']: element for element in dev[start:]}
    dev = dev[:start]
nodes=sorted({n for x in dev for n in x['nodes'] if n!='0'})
# Explicit internal terminals for semiconductor series resistances.
scale={'':1,'t':1e12,'g':1e9,'meg':1e6,'k':1e3,'m':1e-3,'u':1e-6,'n':1e-9,'p':1e-12,'f':1e-15}
def literal(m):return repr(float(m[1])*scale[m[2].lower()])
def expr(s):
    s=s.strip().lower().strip('{}')
    s=re.sub(r'(?<![\w.])([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?)(meg|[tgkmunpf])(?![\w.])',lambda m:literal(m),s)
    s=re.sub(r'(?<![\w.])(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?(?![\w.])',
             lambda m:repr(float(m[0])),s)
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
    def lazy_conditionals(text):
        prefix = 'ternary_fcn('
        start = text.find(prefix)
        if start < 0:
            return text
        args = []
        depth = 0
        begin = start + len(prefix)
        for i in range(begin, len(text)):
            c = text[i]
            if c == '(':
                depth += 1
            elif c == ')' and depth:
                depth -= 1
            elif depth == 0 and c in ',)':
                args.append(lazy_conditionals(text[begin:i]))
                begin = i + 1
                if c == ')':
                    test, yes, no = args
                    return text[:start] + f'(({test}) ? ({yes}) : ({no}))' + lazy_conditionals(text[i+1:])
        raise ValueError(text)
    return lazy_conditionals(tern(s))
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
def val(x,s):
    s=s.strip()
    if s.startswith('{') and s.endswith('}'):
        inner=s[1:-1].lower()
        if inner in x['scope'] and inner not in d['parameters']:s=x['scope'][inner]
    return expr(s)
# Ideal source terminals with no branch-current consumer can be removed
# exactly. Their voltages and every load attached to them are still stamped.
observed = {branch for x in dev for branch in x.get('control_branches', [])}
known = {}
for x in dev:
    name = x['name']
    if x['kind']=='v' and name not in observed:
        a,b=x['nodes']
        raw=x['raw_value']
        source='input' if raw.lower().startswith('sin(') else val(x,re.sub(r'(?i)^dc\s+','',raw))
        if b=='0':
            known[a]=source
        elif a=='0':
            known[b]='-('+source+')'
    elif x['kind']=='e' and name not in observed and x['nodes'][1]=='0':
        if all(n=='0' or n in known for pair in x['control_nodes'] for n in pair):
            def voltage(n):
                return '0.' if n=='0' else '('+known[n]+')'
            coeff=x['coefficients']
            known[x['nodes'][0]]=val(x,coeff[0]) + ''.join(
                f'+({val(x,gain)})*({voltage(a)}-{voltage(b)})'
                for gain,(a,b) in zip(coeff[1:],x['control_nodes'],strict=True))
nodes=[n for n in nodes if n not in known]
stamps=[x for x in dev if not (x['kind'] in 've' and x['name'] not in observed
        and ((x['nodes'][0] in known and x['nodes'][1]=='0')
             or (x['nodes'][0]=='0' and x['nodes'][1] in known)))]
branches={x['name']:len(nodes)+i for i,x in enumerate(x for x in stamps
          if x['kind'] in 'veh' or x['kind']=='b' and x['raw_value'].lower().startswith('v='))}
nid={n:i for i,n in enumerate([*nodes,*branches,*known]) if n not in branches}
def ix(n):return -1 if n=='0' else nid[n]

def bias():
    if seed_path is None:
        return [0.] * (len(nodes)+len(branches))
    header, data = seed_path.read_text().split('Values:', 1)
    labels = [line.split()[1] for line in header.split('Variables:\n', 1)[1].splitlines() if line.strip()]
    values = dict(zip(labels, map(float, data.split()[1:]), strict=True))
    internal = {}
    for x in dev:
        for outer, inner, _ in x.get('series', []):
            terminal = inner.split('#')[-1]
            device = x['name']
            prefix, _, leaf = device.rpartition('.')
            spice_device = f'{leaf[0]}.{prefix}.{leaf}' if prefix else leaf
            suffix = {'rs':'internal' if x['kind']=='d' else 'source',
                      'rd':'drain','rb':'base','rc':'collector','re':'emitter'}[terminal]
            internal[inner] = values.get(f'v({spice_device}#{suffix})', values.get(f'v({outer})', 0.))
    initial = [values.get(f'v({n})', internal.get(n, 0.)) for n in nodes]
    for n in branches:
        prefix, _, leaf = n.rpartition('.')
        name = f'{leaf[0]}.{prefix}.{leaf}' if prefix else leaf
        candidates = [f'i({name})', f'i(a$poly$e.{n})', f'i({n})']
        initial.append(next((values[k] for k in candidates if k in values), 0.))
    return initial

def num(x):return format(x,'.17g')
def params(p,keys,defaults):return ','.join(num(p.get(k,v)) for k,v in zip(keys.split(),defaults))
out=['#pragma once','#include "MNA.hpp"','#include <cmath>','#include <algorithm>','#include <array>','namespace optical {',
'struct Settings { double gdb=0,freq=1000,amp=1,ff=1,ratio=4,soft=1,sens=15,hpf=0,att=.01,rel=.1,makeup=0; };',
'inline double min(double a,double b){return std::min(a,b);} inline double max(double a,double b){return std::max(a,b);}',
'inline double cg(double x){return min(.005,max(2e-8,x));} inline double cs(double x){return min(4,max(0,x));}',
f'inline constexpr int nodeCount={len(nodes)}, unknownCount={len(nodes)+len(branches)}, fixedCount={len(known)};',
f'inline constexpr int outputNode={ix("out")}, preNode={ix("pre")}, grNode={ix("gr")}, cellNode={ix("cell")};',
f'inline constexpr int timingNode={-2 if condense else ix("timing")}, detectorNode={-3 if condense else ix("det")};',
'inline constexpr std::array<double,unknownCount> dcBias{'+','.join(num(x) for x in bias())+'};',
 'inline constexpr std::array<const char*,unknownCount> unknownNames{'+','.join('"'+x+'"' for x in [*nodes,*branches])+'};']
if condense:
    out.append('namespace detector_model {')
    for name in ('ch', 'cfb', 'rhf', 'rhalf', 'chb', 'rrect', 'crect', 'ra', 'ca',
                 'raf', 'caf', 'cp', 'rp', 'rsub', 'csub', 'rif', 'cif', 'crrf', 'cpkfb', 'cpkvel'):
        element = sidechain[name]
        out.append(f'inline constexpr double {name} = {val(element,element["raw_value"])};')
    for name, i in (('eavg',1), ('gpkint',1), ('gpkdamp',1),
                    ('epkedge',1), ('epkedge',2), ('epkedge',3)):
        element = sidechain[name]
        out.append(f'inline constexpr double {name}{i} = {val(element, element["coefficients"][i])};')
    for name, args in (('bcrest',('rect','hfout','pkcomp')), ('bpkshift',('pkedge','hfout')), ('bpkdrive',('pkshift',))):
        body = expr(sidechain[name]['raw_value'].split('=',1)[1])
        body = re.sub(r'\bv\(([^)]+)\)', lambda m:'('+'-'.join(m[1].split(','))+')', body)
        out.append(f'template<class V> inline V {name}(' +
                   ', '.join(f'V {arg}' for arg in args) + ') { using std::abs; return ' + body + '; }')
    out.append('} // namespace detector_model')
known_lines = [f'm.v[{ix(n)}]={value};' for n,value in known.items()]
projection=[]
projected_nodes = set()
linear=[]
dynamic=[]
capacitor_jacobians=[]
# These high-curvature detector equations often reject a Newton trial by
# themselves. Only emit screens for rows whose complete residual is known.
trial_rows = {ix(n): [] for n in ('charge', 'dpc#rs', 'dcl#rs', 'pkamp',
                                  'oa', 'dha#rs', 'dhb#rs') if n in nid}
companion_slot = 0
trial_dependencies = {row: set() for row in trial_rows}
control_names = re.compile(r'\b(gdb|freq|amp|ff|ratio|soft|sens|hpf|att|rel|makeup|input)\b')
for x in stamps:
    k=x['kind'];a,b,*rest=map(ix,x['nodes']);name=x['name'];raw=x['raw_value'];line=[]; series=[]
    outside_base=b
    for p,n,r in x.get('series',[]):
        if k=='q' and n.endswith('#rb'):
            outside_base=ix(p)
        else:
            series.append(f'm.conductance({ix(p)},{ix(n)},{num(1/r)});')
            for row in (ix(p), ix(n)):
                if row in trial_dependencies:
                    trial_dependencies[row].update((ix(p), ix(n)))
    if k=='r':line.append(f'm.conductance({a},{b},1./({val(x,raw)}));')
    elif k=='c':
        capacitance = val(x, raw.split()[0])
        call = 'm.capacitive<FixedCaps>' if condense else 'm.capacitive'
        line.append(f'{call}({a},{b},{capacitance},{companion_slot});')
        if condense:
            if control_names.search(capacitance):
                raise ValueError(f'Parameter-dependent capacitance needs a dynamic stamp: {name}')
            capacitor_jacobians.append(f'm.derivative({a},{b},{a},{b},inverse*({capacitance}));')
        companion_slot += 1
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
        t,p=model(x['model'])
        spec = params(p,"is n bv ibv cjo vj m fc tt",[1e-14,1,0,.001,0,1,.5,.5,0])
        stores_charge = p.get('cjo',0) != 0 or p.get('tt',0) != 0
        slot = companion_slot if stores_charge else -1
        if stores_charge:
            companion_slot += 1
        line.append(f'm.diodeDevice({a},{b},{{{spec}}},{slot});')
        for row, sign in ((a, '-'), (b, '+')):
            if row in trial_rows:
                trial_dependencies[row].update((a, b))
                trial_rows[row].append(f'r {sign}= m.trialDiodeCurrent({a},{b},{{{spec}}},{slot},inverse);')
    elif k=='q':
        t,p=model(x['model'])
        stores_charge = any(p.get(key,0) != 0 for key in ('cje','cjc','tf','tr'))
        slot = companion_slot if stores_charge else -1
        if stores_charge:
            companion_slot += 2
        line.append(f'm.bjtDevice({a},{b},{rest[0]},{outside_base},{-1 if t=="pnp" else 1},{{{params(p,"is bf nf br nr vaf var ikf ikr ise ne isc nc nk cje vje mje cjc vjc mjc fc tf tr xtf itf vtf rb rbm irb",[1e-16,100,1,1,1,0,0,0,0,0,1.5,0,2,.5,0,.75,.33,0,.75,.33,.5,0,0,0,0,0,0,p.get('rb',0),0])}}},{slot});')
    elif k=='j':
        t,p=model(x['model']);line.append(f'm.jfet({a},{b},{rest[0]},{-1 if t=="pjf" else 1},{params(p,"vto beta is lambda",[-2,1e-4,1e-14,0])});')
    elif k=='s':
        t,p=model(x['model']);line.append(f'm.conductance({a},{b},v({rest[0]})-v({rest[1]})>{num(p.get("vt",0))}?{num(1/p.get("ron",1))}:{num(1/p.get("roff",1e12))});')
    elif condense and name == 'btarget':
        key = branches[name]
        projection.append(f'm.v[{a}]=target;')
        line.append(f'm.source({a},{b},{key},target);m.derivative({key},-1,preNode,-1,-targetPre,true);'
                    f'm.derivative({key},-1,grNode,-1,-targetGr,true);')
    elif name == 'bldr' and 'ldr_even_30k' in d['parameters']:
        even = ','.join(expr(d['parameters'][key])
                        for key in ('ldr_even_30k', 'ldr_even_3k', 'ldr_even_knee'))
        line.append(f'm.photocell({a},{b},{ix("cell")},{even});')
    elif k=='b':
        mode,body=raw.split('=',1)
        refs=[]
        def replace(m):
            ns=m[1].lower().split(',');ids=[ix(n.strip()) for n in ns]
            for n in ids:
                if 0<=n<len(nodes) and n not in refs:refs.append(n)
            return '('+'-'.join(f'w({n})' for n in ids)+')'
        body=expr(body)
        body=re.sub(r'(?i)\bv\(([^)]+)\)',replace,body)
        key=branches[name] if mode.lower()=='v' else -1
        if key>=0:
            assert b==-1
            projected_nodes.add(a)
            projection.append(f'm.v[{a}]={body}; // {name}')
        line.append(f'm.behavior({a},{b},{key},{{{",".join(map(str,refs))}}},[&](auto w){{return {body};}});')
    else:raise ValueError(x)
    # Series resistance and parameter-independent linear elements have
    # a constant Jacobian. Known terminal voltages still vary per sample.
    coefficient = ' '.join(line)
    linear.extend(series)
    if k in 'rievfgh' and not control_names.search(coefficient):
        linear.extend(line)
        dependencies = [a, b, *(ix(n) for pair in x.get('control_nodes', []) for n in pair)]
        for row in (a, b):
            if row in trial_dependencies:
                trial_dependencies[row].update(dependencies)
    else:
        dynamic.extend(line)
        if k != 'd':
            terminals = [a,b,*rest] if k in 'qj' else [a,b]
            for terminal in terminals:
                trial_rows.pop(terminal, None)
            if k == 'q':
                trial_rows.pop(outside_base, None)
out+=['inline void stampFixed(LinearStamp& m) {', 'auto v=[&](int n){return m.voltage(n);};']
out+=linear
out+=['}', 'inline void setKnownTerminals(MNA& m,const Settings& s,double input) {',
      'using std::exp;using std::log; const auto [gdb,freq,amp,ff,ratio,soft,sens,hpf,att,rel,makeup]=s;']
out+=known_lines
if condense:
    out += ['}', 'inline void stampCapacitorJacobian(MNA& m,double inverse) {', *capacitor_jacobians]
control_args = ',double target,double targetPre,double targetGr' if condense else ''
project_args = ',double target' if condense else ''
template = 'template<bool FixedCaps> ' if condense else ''
out+=['}', f'{template}inline void stampDynamic(MNA& m,const Settings& s,double inverse{control_args}) {{',
      'using std::abs;using std::exp;using std::log;using std::log10;using std::pow;using std::sqrt;',
      'const auto [gdb,freq,amp,ff,ratio,soft,sens,hpf,att,rel,makeup]=s;',
      'auto v=[&](int n){return m.voltage(n);};m.beginCompanions(inverse);']
out+=dynamic
out+=['m.endCompanions();','}', f'inline constexpr int companionCount = {companion_slot};']
out+=[f'inline void projectBehavior(MNA& m,const Settings& s{project_args}) {{',
      'using std::abs;using std::exp;using std::log10;using std::pow;using std::sqrt;',
      'const auto& [gdb,freq,amp,ff,ratio,soft,sens,hpf,att,rel,makeup]=s;',
      'auto w=[&](int n){return m.voltage(n);};']
out+=projection
out+=['}', 'inline bool rejectTrial(const MNA& m, const LinearStamp& linear, double inverse, double limit) {']
for row, contributions in trial_rows.items():
    if trial_dependencies[row] & projected_nodes:
        continue
    out += [f'    {{ double r = linear.residualForRow({row}, m.v.data());',
            *('      '+line for line in contributions),
            '      if (r*r > limit) return true; }']
out += ['    return false;', '}', '}']
output_path.write_text('\n'.join(out)+'\n')
print('Generated full circuit stamps:',len(dev),'unknowns:',len(nodes)+len(branches))
