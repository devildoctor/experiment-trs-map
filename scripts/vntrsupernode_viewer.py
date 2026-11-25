"""Embed a pinned Cytoscape runtime and path data for direct file:// use."""
from collections import defaultdict
from html import escape
from pathlib import Path
import json


BASE_COLORS={'A':'#00FF00','C':'#0000FF','G':'#D17105','T':'#FF0000','N':'#B6B6B6'}


def motif_aliases(model):
    if not model:
        return {}
    names={m['id']:f'U{i}' for i,m in enumerate(model['motifs'],start=1)}
    mapping=model.get('efficient_motif_selection',{}).get('mapping',{})
    return {**names,**{raw:names[target] for raw,target in mapping.items()}}


def dictionary_table(model):
    if not model or not model['motifs']:
        return ''
    aliases=motif_aliases(model)
    motifs=model['motifs']
    aligned=all('aligned_sequence' in m and m['aligned_sequence'].replace('-','')==m['sequence']
                for m in motifs) and len({len(m['aligned_sequence']) for m in motifs})==1
    sequences=[m['aligned_sequence'] if aligned else m['sequence'] for m in motifs]
    columns=range(1,max(map(len,sequences))+1)
    caption='重复单元字典 · 对齐列' if aligned else '重复单元字典'
    header=''.join(f'<th scope="col">{i}</th>' for i in columns)
    rows=[]
    for motif,sequence in zip(motifs,sequences):
        name=escape(aliases[motif['id']])
        template_bases=escape(motif['sequence'])
        cells=[]
        for base in sequence:
            color=BASE_COLORS.get(base,'transparent')
            ink='white' if base=='C' else 'black'
            cells.append(f'<td style="background-color:{color};color:{ink}">{escape(base)}</td>')
        rows.append(f'<tr data-motif="{name}" data-template-bases="{template_bases}"><th scope="row">{name}</th>{"".join(cells)}</tr>')
    return (f'<div id="fiveDictionary"><table><caption>{caption}</caption>'
            f'<thead><tr><th>单元</th>{header}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>')


def write_viewer(out,display,visits,variants,records,events,policy,rows,*,display_model=None):
    root=Path(__file__).resolve().parents[1]
    from vntrsupernode_export import path_budget
    samples={}
    path_info={row['path']:{k:row[k] for k in ('role','length','edited_bp','fraction',
               'represented_is_observed_whole_sequence')} for row in rows}
    for r in records:
        active=set(r['molecules']) if r['role']!='reference' else set()
        path_info[r['id']].update(budget_bp=path_budget(r,policy),molecule_keys=sorted(active))
        known={o['molecule'] for o in r['observations'] if isinstance(o.get('core_ambiguous'),bool)} & active
        tied={o['molecule'] for o in r['observations'] if o.get('core_ambiguous') is True} & active
        path_info[r['id']].update(boundary_assessed=len(known),boundary_tied=len(tied))
        by_sample=defaultdict(set)
        for o in r['observations']:
            if o['molecule'] in active: by_sample[o['sample']].add(o['molecule'])
        samples[r['id']]={s:sorted(ms) for s,ms in sorted(by_sample.items())}
    bundle=dict(graph=display,paths=visits,variants=variants,samples=samples,path_info=path_info,
                losses={e['id']:e for es in events.values() for e in es},policy=policy)
    data=json.dumps(bundle,ensure_ascii=False,separators=(',',':')).replace('<','\\u003c')
    library=root/'third_party/cytoscape-3.30.4'
    html=(root/'scripts/supernodes.html').read_text(encoding='utf8')
    if display_model is None:
        for filename in ('locus_model.efficient.json','locus_model.original.json'):
            if (out/filename).is_file():
                display_model=json.loads((out/filename).read_bytes())
                break
    html=html.replace('/* MOTIF ALIASES */',json.dumps(motif_aliases(display_model)).replace('<','\\u003c'))
    html=html.replace('<!-- MOTIF DICTIONARY -->',dictionary_table(display_model))
    html=html.replace('/* PROJECTION */',(root/'scripts/supernode_projection.js').read_text(encoding='utf8'))
    html=html.replace('/* COMPARE */',(root/'scripts/supernode_compare.js').read_text(encoding='utf8'))
    html=html.replace('/* CYTOSCAPE */',(library/'cytoscape.min.js').read_text(encoding='utf8')).replace('/* DATA */',data)
    (out/(display['locus']+'.supernodes.html')).write_text(html,encoding='utf8')

