"""Verify an actual run against input BAMs and exact GFA/SVG path sequences."""
import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

# Stored ONTUL SEQ can exceed the csv module's 128 KiB default field limit.
csv.field_size_limit(64 * 1024 * 1024)

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrgraph import IndexedFastaReference, read_catalog


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--inventory',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--catalog',type=Path,required=True)
    p.add_argument('--samtools',required=True)
    args=p.parse_args()
    out=args.run/'results'/'D1S80'
    status=json.loads((out/'status.json').read_text())
    assert not status['errors'],status
    decisions=[json.loads(line) for line in (out/'anchor_decisions.jsonl').read_text().splitlines()]
    decision_keys={r['molecule_key']:r for r in decisions}
    assert len(decision_keys)==len(decisions),'duplicate decision molecule'
    expected_sources={s for d in decisions for s in d['sources']}
    records,hashes={},[]
    inventory=list(csv.DictReader(args.inventory.open(),delimiter='\t'))
    for item in inventory:
        bam=Path(item['bam'])
        digest=sha256()
        with bam.open('rb') as f:
            for block in iter(lambda:f.read(4*1024*1024),b''):digest.update(block)
        hashes.append(dict(bam=str(bam),sha256=digest.hexdigest(),bytes=bam.stat().st_size))
        if item['locus']!='D1S80':continue
        sam=subprocess.check_output([args.samtools,'view','-h',str(bam)],text=True)
        for i,line in enumerate(sam.splitlines(),1):
            if line.startswith('@'):continue
            fields=line.split('\t')
            tags={x[0]:x[2] for t in fields[11:] if len(x:=t.split(':',2))==3}
            source=f'{bam}:record{i}:flag{int(fields[1])}:RG={tags.get("RG", "")}'
            records[source]=fields
    assert set(records)==expected_sources,'input records lost from audit'
    assert len(records)==status['input_records']
    molecules=list(csv.DictReader((out/'molecules.tsv').open(),delimiter='\t'))
    assert len({r['molecule_key'] for r in molecules})==len(molecules)==status['molecules']
    for m in molecules:
        assert decision_keys[m['molecule_key']]['status']=='complete_support'
        start,end=int(m['qL']),int(m['qR'])
        assert m['query_sequence'][start:end]==m['sequence']
        sources=json.loads(m['records'])
        assert any(records[s][9].upper()==m['query_sequence'] for s in sources),'stored SEQ not found in BAM'
        assert set(m['sequence'])<=set('ACGT')
    alleles={r['allele']:r for r in csv.DictReader((out/'audit.tsv').open(),delimiter='\t')}
    segments,paths={},{}
    gfa=out/'D1S80.gfa'
    if gfa.exists():
        for line in gfa.read_text().splitlines():
            f=line.split('\t')
            if f[0]=='S':segments[f[1]]=f[2]
            if f[0]=='P':paths[f[1]]=[x[:-1] for x in f[2].split(',')]
    assert set(paths)==set(alleles)
    locus=next(x for x in read_catalog(args.catalog) if x.name=='D1S80')
    with IndexedFastaReference(args.reference) as ref:
        left=ref.fetch(locus.chrom,locus.start-25,locus.start)
        right=ref.fetch(locus.chrom,locus.end,locus.end+25)
    molecule_ids={r['molecule_key'] for r in molecules}
    for name,ids in paths.items():
        row=alleles[name]
        assert ''.join(segments[i] for i in ids)==left+row['sequence']+right
        support=json.loads(row['molecule_ids'])
        assert len(support)==len(set(support))==int(row['read_support'])
        assert set(support)<=molecule_ids
    svg=ET.parse(out/'D1S80.svg')
    rendered={name:'' for name in alleles}
    for e in svg.iter():
        if 'data-path' in e.attrib:
            rendered[e.attrib['data-path']]+=e.attrib['data-exact']*int(e.attrib['data-count'])
    assert all(rendered[name]==row['sequence'] for name,row in alleles.items())
    summary=dict(status='verified',input_bams=sum(r['locus']=='D1S80' for r in inventory),
                 input_records=len(records),audited_molecules=len(decisions),complete_molecules=len(molecules),
                 paths=len(paths),molecule_decisions=dict(Counter(d['status'] for d in decisions)),
                 per_sample=dict(Counter(m['sample'] for m in molecules)),
                 checks=['all_BAM_records_audited','one_vote_per_molecule','exact_stored_SEQ_slice']+
                        (['GFA_exact_path_reconstruction','SVG_exact_path_reconstruction'] if paths else
                         ['no_complete_GFA_paths_without_complete_support']))
    (args.run/'input_sha256.json').write_text(json.dumps(hashes,indent=2),encoding='utf-8')
    (out/'verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False))


if __name__=='__main__':main()
