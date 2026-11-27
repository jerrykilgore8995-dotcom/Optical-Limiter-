"""Generate a typed, compile-time C++ passive-network deck from flattened SPICE IR.
Dynamic/behavioral devices are enumerated explicitly as pending, never discarded.
"""
import json,re
from pathlib import Path
root=Path(__file__).resolve().parent
d=json.loads((root/'deck.json').read_text())
def number(v,scope):
    v=v.strip()
    if v.startswith('{') and v.endswith('}'):
        v=scope.get(v[1:-1].lower(),v)
    if '=' in v or '{' in v or '}' in v or '(' in v:return None
    m=re.fullmatch(r'([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?)(meg|[tgkmunpf])?',v,re.I)
    if not m:return None
    x=float(m.group(1));u=(m.group(2) or '').lower()
    return x*{'':1,'t':1e12,'g':1e9,'meg':1e6,'k':1e3,'m':1e-3,'u':1e-6,'n':1e-9,'p':1e-12,'f':1e-15}[u]
nodes=sorted({n for x in d['devices'] for n in x['nodes'] if n!='0'})
nid={n:i for i,n in enumerate(nodes)}
rows=[];caps=[];unsupported=[]
for x in d['devices']:
    k=x['kind'];v=number(x['raw_value'],x['scope'])
    if k=='r':
        if v is not None and v>0: rows.append((nid[x['nodes'][0]] if x['nodes'][0]!='0' else -1,nid[x['nodes'][1]] if x['nodes'][1]!='0' else -1,v))
        else:unsupported.append((x['name'],'resistor value expression',x['raw_value']))
    elif k=='c':
        if v is not None and v>=0: caps.append((nid[x['nodes'][0]] if x['nodes'][0]!='0' else -1,nid[x['nodes'][1]] if x['nodes'][1]!='0' else -1,v))
        else:unsupported.append((x['name'],'capacitor value expression',x['raw_value']))
    else:unsupported.append((x['name'],k,x.get('raw_value','')))
out=['#pragma once','#include <array>','#include <cstdint>','namespace optical::generated {',
     'struct Passive { int a,b; double value; };',
     f'inline constexpr int nodeCount={len(nodes)};',
     f'inline constexpr std::array<Passive,{len(rows)}> resistors={{{{']
out += [f'{{{a},{b},{v:.17g}}},' for a,b,v in rows]
out += ['}};']
out += [f'inline constexpr std::array<Passive,{len(caps)}> capacitors={{{{']
out += [f'{{{a},{b},{v:.17g}}},' for a,b,v in caps]
out += ['}};']
out += [f'inline constexpr int unsupportedCount={len(unsupported)};','}']
(root/'DeckGenerated.hpp').write_text('\n'.join(out)+'\n')
(root/'unsupported.json').write_text(json.dumps(unsupported,indent=2))
print(f'Generated {len(nodes)} nodes, {len(rows)} resistors, {len(caps)} capacitors; {len(unsupported)} dynamic/non-passive elements explicitly pending')
