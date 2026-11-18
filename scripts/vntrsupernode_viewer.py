"""Embed a pinned Cytoscape runtime and path data for direct file:// use."""
from collections import defaultdict
from pathlib import Path
import json


def write_viewer(out,display,visits,variants,records,events,policy,rows):
    root=Path(__file__).resolve().parents[1]
    from vntrsupernode_export import path_budget
    samples={}
    path_info={row['path']:{k:row[k] for k in ('role','length','edited_bp','fraction',
               'represented_is_observed_whole_sequence')} for row in rows}
    for r in records:
        active=set(r['molecules']) if r['role']!='reference' else set()
        path_info[r['id']].update(budget_bp=path_budget(r,policy),molecule_keys=sorted(active))
        by_sample=defaultdict(set)
        for o in r['observations']:
            if o['molecule'] in active: by_sample[o['sample']].add(o['molecule'])
        samples[r['id']]={s:sorted(ms) for s,ms in sorted(by_sample.items())}
    bundle=dict(graph=display,paths=visits,variants=variants,samples=samples,path_info=path_info,
                losses={e['id']:e for es in events.values() for e in es},policy=policy)
    data=json.dumps(bundle,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    library=root/'third_party/cytoscape-3.30.4'
    html=(root/'scripts/supernodes.html').read_text(encoding='utf8')
    html=html.replace('/* COMPARE */',(root/'scripts/supernode_compare.js').read_text(encoding='utf8'))
    html=html.replace('/* CYTOSCAPE */',(library/'cytoscape.min.js').read_text(encoding='utf8')).replace('/* DATA */',data)
    (out/(display['locus']+'.supernodes.html')).write_text(html,encoding='utf8')

