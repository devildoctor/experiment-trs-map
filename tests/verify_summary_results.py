"""Cross-check final artifacts independently of summary presentation."""
import json
from pathlib import Path
import sys

root=Path(sys.argv[1]);d=json.loads((root/(root.name+'.structures.json')).read_text(encoding='utf-8'))
assert [m['id'] for m in d['motifs']]==[f'M{i+1}' for i in range(len(d['motifs']))]
assert all(m['base']==m['id'] for m in d['motifs'])
fragments={f['id']:f for r in d['rows'] for f in r['fragments']}
fragment_rows={f['id']:r for r in d['rows'] for f in r['fragments']}


def molecule(assignment):
    row=fragment_rows[assignment['fragment_id']]
    assert assignment['sample']==row['sample'] and assignment['read_id']==row['read_id']
    key=(row['sample'],row.get('namespace',''),row['read_id'])
    stored=json.loads(assignment['molecule'])
    # 兼容旧存档的双段键；新结果必须保留 namespace，不能忽略非空值。
    if len(stored)==2:
        assert 'namespace' not in row and 'namespace' not in assignment
        stored.insert(1,'')
    assert tuple(stored)==key
    assert assignment.get('namespace','')==key[1]
    return key


summary=d['summary'];assignments=summary['assignments']
assert len(assignments)==len(fragments)
assert len({a['fragment_id'] for a in assignments})==len(fragments)
for f in fragments.values():assert ''.join(t['exact_bases'] for t in f['tokens'])==f['sequence']
for assignment in assignments:molecule(assignment)
for group in summary['classes']:
    members=[a for a in assignments if a['class_id']==group['id']]
    assert group['support']==len({molecule(a) for a in members})
    assert all(all(t['status']=='unit' for t in fragments[a['fragment_id']]['tokens']) for a in members)
    for subtype in group['subtypes']:
        example=subtype['representative'];f=fragments[example['fragment_id']]
        assert example['sequence']==f['sequence']
        assert subtype['counts']==[b['count'] for b in f['blocks']]
for module in summary['modules']:
    assert module['support']==len({molecule(a) for a in module['instances']})
    assert module['occurrences']==len({(molecule(a),fragments[a['fragment_id']]['qstart']+a['start'],a['sequence']) for a in module['instances']})
    for a in module['instances']:
        assert fragments[a['fragment_id']]['sequence'][a['start']:a['end']]==a['sequence']
reference=d['reference_structure']
assert ''.join(t['exact_bases'] for t in reference['tokens'])==reference['sequence']
assert len(reference['sequence'])==d['locus']['end']-d['locus']['start']
assert len(d['rows'])==d['audit']['read_rows']
for name in ('reference_structure.json','structure_classes.tsv','local_modules.tsv','read_assignments.tsv'):
    assert (root/name).is_file()
print('STRUCTURE_SUMMARY_VERIFIED',root.name,len(d['rows']),len(fragments),summary['audit'],flush=True)
