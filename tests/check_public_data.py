#!/usr/bin/env python3
"""No-dependency checks of the already generated public atlas.js and its audit."""
import json
from pathlib import Path
root=Path(__file__).resolve().parent.parent
text=(root/'data/atlas.js').read_text(encoding='utf8')
assert text.startswith('/* pAML aggregate')
assert '\nwindow.PAML_METHYLATION_ATLAS=' in text
D=json.loads(text.split('window.PAML_METHYLATION_ATLAS=',1)[1].removesuffix(';\n'))
audit=json.loads((root/'data/build-audit.json').read_text())
assert D['schema']=='pAML-methylation-v2'
assert len(D['rows'])==D['meta']['associations']==79920
assert len(D['probes'])==D['meta']['uniqueCpGs']==50013
assert len(D['directions'])==119
assert len(D['groups'])==23
assert len(D['clinical'])==1254
assert len(D['flags'])==22
assert len(D['geneRows'])==D['meta']['geneLabels']==20277
assert D['meta']['uniqueCpGsByRegion']=={'Promoter':15048,'Gene body':26012,'Intergenic':8953}
assert sum(D['meta']['uniqueCpGsByRegion'].values())==50013
assert audit['input_11d_12c_exact_order_match']
assert audit['group_probe_original_mean_records']==134073
assert audit['group_mean_consistency_cross_checks']>25000
assert sum(1 for c in D['clinical'] if c[3]!='eligible' and c[4:7]==[None,None,None])==962
assert all((0<=r[2]<=1) if r[2] is not None else True for r in D['rows'])
assert all(r[5] in ('Promoter','Gene body','Intergenic') for r in D['probes'])
assert all(len(r)==28 for r in D['rows'])
assert all(len(r)==25 for r in D['clinical'])
assert all(len(r)==11 for r in D['probes'])
assert not ({'USI','Patient_ID','Individual_Beta','Patient_H'}&set(D))
print('PASS: complete 11D/12C keys, data schema, 23 groups, all contexts, full suite flags, 13C redaction.')
