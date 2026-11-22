"""Audit persisted graph, TXT, reference, SVG and original BAM slices."""
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrflankalign import core_eligible
from vntrpathgraph import validate_gfa,digest
from vntrgraph import IndexedFastaReference,read_catalog


def verify(root,baseline,reference,catalog):
    summary=json.loads((root/'summary.json').read_text());locus=summary['locus']
    paths=json.loads((root/'source_members.json').read_text())
    raw={(r['row_id'],r['record_id']):r for r in map(json.loads,(root/'records.jsonl').read_text().splitlines())}
    passed={key:r for key,r in raw.items() if core_eligible(r)}
    observed=[p for p in paths if p['role']=='observation'];refs=[p for p in paths if p['role']=='reference']
    assert [p['id'] for p in paths]==[f'P{i:04d}' for i in range(1,len(paths)+1)]
    assert len(refs)==1 and refs[0]['id']=='P0001' and not refs[0]['molecules']
    interval=next(l for l in read_catalog(catalog) if l.name==locus)
    with IndexedFastaReference(reference) as fasta:ref=fasta.fetch(interval.chrom,interval.start,interval.end)
    assert refs[0]['sequence']==ref
    assert {p['sequence'] for p in observed}=={r['sequence'] for r in passed.values()}
    assert {(o['row_id'],o['record_id']) for p in observed for o in p['observations']}==set(passed)
    for p in observed:
        assert p['sha256']==digest(p['sequence'])
        for o in p['observations']:
            r=raw[(o['row_id'],o['record_id'])]
            assert r['query_sequence'].upper()[o['qstart']:o['qend']]==p['sequence']
            assert r['molecule_key']==o['molecule']
            assert (r.get('mapq'),r.get('flag'),r.get('alignment_type'))==(o['mapq'],o['flag'],o['alignment_type'])
    expected={p['id']:p['sequence'] for p in paths}
    for kind in ['exact','structure','preview']:validate_gfa(root/f'{locus}.{kind}.gfa',expected)
    with (root/f'{locus}.path_sequences.txt').open() as handle:
        txt=list(csv.DictReader((s for s in handle if not s.startswith('#')),delimiter='\t'))
    assert {p['name']:p['sequence'] for p in txt}==expected
    assert all(int(p['length_bp'])==len(p['sequence']) for p in txt)
    for filename in ['paths.original.fa','paths.represented.fa']:
        seqs={}
        for line in (root/filename).read_text().splitlines():
            if line.startswith('>'):name=line[1:];seqs[name]=''
            else:seqs[name]+=line
        assert seqs==expected
    old=json.loads((baseline/'source_members.json').read_text())
    assert {(r['sequence'],m) for r in old for m in r['molecules']}<={(r['sequence'],m) for r in observed for m in r['molecules']}
    previous={(r['row_id'],r['record_id']):r for r in map(json.loads,(baseline/'records.jsonl').read_text().splitlines())}
    assert raw==previous,'Raw extraction changed'
    nodes=set();edges=set();walks={};node_support={};edge_support={}
    for line in (root/f'{locus}.structure.gfa').read_text().splitlines():
        f=line.split('\t')
        if f[0]=='S':nodes.add(f[1]);node_support[f[1]]=int(next(x[5:] for x in f if x.startswith('mc:i:')))
        elif f[0]=='L':
            edge=(f[1],f[3]);edges.add(edge);edge_support[edge]=int(next(x[5:] for x in f if x.startswith('mc:i:')))
        elif f[0]=='P':walks[f[1]]=[n[:-1] for n in f[2].split(',')]
    members={p['id']:set(p['molecules']) for p in paths}
    for node in nodes:assert node_support[node]==len(set().union(*(members[p] for p,w in walks.items() if node in w)))
    for edge in edges:assert edge_support[edge]==len(set().union(*(members[p] for p,w in walks.items() if edge in set(zip(w,w[1:])))))
    ns={'s':'http://www.w3.org/2000/svg'};svg=ET.parse(root/f'{locus}.preview.svg');drawn_nodes=[];drawn_edges=[]
    for element in svg.findall('.//s:g',ns):
        title=element.find('s:title',ns)
        if title is None or not title.text:continue
        if element.get('class')=='node':drawn_nodes.append(title.text[1:])
        if element.get('class')=='edge':
            a,b=title.text.split('->');drawn_edges.append((a[1:],b[1:]))
    assert set(drawn_nodes)==nodes and len(drawn_nodes)==len(nodes)
    assert set(drawn_edges)==edges and len(drawn_edges)==len(edges)
    assert summary['preview']['panels']==1 and summary['preview']['omitted_paths']==0
    report=dict(locus=locus,passed=True,input_records=len(raw),B2_passing_records=len(passed),
                observed_sequences=len(observed),reference_paths=1,total_paths=len(paths),
                molecules=len(set().union(*(members[p['id']] for p in observed))),nodes=len(nodes),edges=len(edges),
                reference_length=len(ref),all_B2_records_included=True,raw_extraction_unchanged=True,
                old_members_preserved=True,reference_exact=True,TXT_FASTA_GFA_identical=True,
                SVG_full_topology=True,support_independently_deduplicated=True)
    (root/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
    return report


if __name__=='__main__':verify(*(Path(s) for s in sys.argv[1:5]))
