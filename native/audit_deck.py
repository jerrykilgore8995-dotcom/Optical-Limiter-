"""Audit flattened optical-limiter circuit before native MNA code generation."""
import json,sys
from collections import Counter
from pathlib import Path
p=Path(__file__).with_name('deck.json')
d=json.loads(p.read_text())
required={'r','c','l','v','i','d','q','j','e','g','f','h','b','s'}
counts=Counter(x['kind'] for x in d['devices'])
assert set(counts)<=required, set(counts)-required
assert counts['q']>=8 and counts['d']>=90 and counts['j']>=30
assert counts['b']==14 and counts['s']==2
assert len(d['models'])>=40
assert any('OPA2134' in k.upper() for k in d['models'])
assert any('KSC1845' in k.upper() for k in d['models'])
assert any('KSA992' in k.upper() for k in d['models'])
assert any('D4148W' in k.upper() for k in d['models'])
nodes={n for x in d['devices'] for n in x['nodes']}
assert '0' in nodes
# Explicitly identify nonlinear/behavioral constructs that require custom stamps.
beh=[]
for x in d['devices']:
    if x['kind']=='b':
        beh.append({'name':x['name'],'expression':x.get('expression',x['raw_value'])})
    elif x['kind'] in 'efgh':
        beh.append({'name':x['name'],'kind':x['kind'],'controls':x.get('control_nodes',x.get('control_branches',[])),'coefficients':x.get('coefficients',[])})
out={'device_counts':dict(counts),'total_devices':len(d['devices']),'model_count':len(d['models']),'node_count':len(nodes),'behavioral_stamps':beh,'model_names':sorted(d['models'])}
Path('native/audit.json').write_text(json.dumps(out,indent=2))
print('PASS: all device families retained; models:',len(d['models']),'devices:',len(d['devices']),'nodes:',len(nodes))
print('Behavioral and controlled-source stamps requiring native implementation:',len(beh))
for x in beh:print(x['name'],x.get('kind','B'),x.get('expression',x.get('coefficients')))
