"""Independent real-file audit, including stored SAM slices and reversible loss."""
from collections import Counter
from hashlib import sha256
from pathlib import Path
import json
import sys
import xml.etree.ElementTree as ET

from verify_pathgraph_results import fasta,verify_gfa


def verify(root,baseline):
    summary=json.loads((root/'summary.json').read_text())
    old=json.loads((baseline/'summary.json').read_text())
    originals=fasta(root/'paths.original.fa');represented=fasta(root/'paths.represented.fa')
    previous=fasta(baseline/'paths.original.fa')
    records=json.loads((root/'source_members.json').read_text())
    members={r['id']:set(r['molecules']) for r in records}
    old_records=json.loads((baseline/'source_members.json').read_text())
    assert originals==previous
    assert members=={r['id']:set(r['molecules']) for r in old_records}
    assert summary['paths']==len(originals) and summary['molecules']==sum(map(len,members.values()))
    raw={}
    for line in (root/'records.jsonl').read_text().splitlines():
        r=json.loads(line);raw[(r['row_id'],r['record_id'])]=r
    for r in records:
        assert sha256(r['sequence'].encode()).hexdigest()==r['sha256']
        for o in r['observations']:
            sam=raw[(o['row_id'],o['record_id'])]
            assert sam['query_sequence'].upper()[o['qstart']:o['qend']]==originals[r['id']]
            assert sam['molecule_key']==o['molecule']
    losses={r['path']:r for r in map(json.loads,(root/'losses.jsonl').read_text().splitlines())}
    features={r['path']:r for r in map(json.loads,(root/'path_features.jsonl').read_text().splitlines())}
    for p,seq in represented.items():
        loss=losses.get(p,{});events=loss.get('events',[]);cursor=0;restored=[]
        for e in events:
            assert cursor<=e['rstart']<=e['rend']<=len(seq)
            assert seq[e['rstart']:e['rend']]==e['represented']
            restored += [seq[cursor:e['rstart']],e['original']];cursor=e['rend']
        restored.append(seq[cursor:]);assert ''.join(restored)==originals[p]
        if loss:
            assert 0<loss['loss_bp']==len(events)<=summary['max_loss_bp']
            counts=Counter(e['op'] for e in events)
            assert len(seq)==len(originals[p])+counts['I']-counts['D']
            assert features[p]['signature']==features[loss['target_path']]['signature']
            assert features[p]['merge_eligible'] and features[loss['target_path']]['merge_eligible']
    locus=summary['locus'];graphs={}
    for name,seqs in [('exact',originals),('structure',represented),
                     ('preview',{p:represented[p] for p in summary['preview']['paths']})]:
        graphs[name]=verify_gfa(root/f'{locus}.{name}.gfa',seqs,members)
        assert graphs[name]==summary['graph_stats'][name]
    # The drawn nodes and edges must match the selected GFA, including panel ports.
    ns={'s':'http://www.w3.org/2000/svg'};svg=ET.parse(root/f'{locus}.preview.svg')
    nodes=set();edges=set()
    for line in (root/f'{locus}.preview.gfa').read_text().splitlines():
        f=line.split('\t')
        if f[0]=='S':nodes.add(f[1])
        if f[0]=='L':edges.add((f[1],f[3]))
    drawn_nodes=Counter();drawn_edges=Counter()
    for element in svg.findall('.//s:g',ns):
        title=element.find('s:title',ns)
        if title is None or title.text is None:continue
        label=title.text
        if element.get('class')=='node' and label.startswith('s') and label[1:].isdigit():drawn_nodes[label[1:]]+=1
        if element.get('class')=='edge':
            a,b=label.split('->')
            if a.startswith('port_in_'):edge=tuple(a[len('port_in_'):].split('_'))
            elif b.startswith('port_out_'):edge=tuple(b[len('port_out_'):].split('_'))
            else:edge=(a[1:],b[1:])
            drawn_edges[edge]+=1
    assert set(drawn_nodes)==nodes and all(n==1 for n in drawn_nodes.values())
    assert set(drawn_edges)==edges
    report=dict(locus=locus,status='passed',paths=len(originals),molecules=summary['molecules'],
        observation_set_unchanged=True,raw_sam_slices=True,represented_and_original_reconstruction=True,
        protected_structure=True,static_svg_topology=True,lossy_paths=len(losses),loss_bp=summary['loss_bp'],
        previous_sequence_graph=old['graph_stats']['structure'],structure_graph=graphs['structure'],
        exact_structure_graph=graphs['exact'],annotation_status=summary['annotation_status'],
        scope='Sequence, source, display structure and loss consistency; not a biological truth or genotype validation.')
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':verify(Path(sys.argv[1]),Path(sys.argv[2]))
