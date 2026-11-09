"""Run the B2 extraction, frozen decomposition and complete structure graph."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


def read(path):return json.loads(path.read_text())
def lines(path):return [json.loads(s) for s in path.read_text().splitlines()]
def digest(path):return sha256(path.read_bytes()).hexdigest()

def source_signature(folder):
    return sorted((r['sha256'],sorted(r['molecules'])) for r in read(folder/'source_members.json'))

def cut_signature(folder):
    return {(r['row_id'],r['record_id']):(r.get('qstart'),r.get('qend'),r.get('sequence'),
             r.get('core_status'),r.get('core_ambiguous'),r.get('core_anchor_evidence',{}).get('pairs'))
            for r in lines(folder/'records.jsonl')}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--previous-b2',type=Path,required=True)
    p.add_argument('--samtools',required=True)
    args=p.parse_args();source=Path(__file__).resolve().parents[1]
    args.out.mkdir(parents=True,exist_ok=False)
    reports=[]
    for locus in ('ATXN2','D1S80'):
        common=['--reference',args.project/'reference/GRCh38_full_analysis_set_plus_decoy_hla.fa',
                '--catalog',source/'catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed',
                '--locus',locus,'--samtools',args.samtools]
        model=source/f'catalogs/{locus}.structure.v1.json'
        graph=['--graph-mode','structure','--locus-model',model,
               '--structure-policy',source/f'catalogs/{locus}.structure-policy.v1.json','--max-loss-bp','0']
        commands=[
            [sys.executable,source/'scripts/vntrstructure.py','--bam-dir',args.project/'bamdata_by_locus',
             '--out',args.out/'structures','--locus-model',model,'--workers','2',*common],
            [sys.executable,source/'scripts/vntrpathgraph_cli.py','--structures-json',
             args.out/f'structures/{locus}/{locus}.structures.json','--out',args.out/'import_audit',*common,*graph],
            [sys.executable,source/'scripts/vntrpathgraph_cli.py','--bam-dir',args.project/'bamdata_by_locus',
             '--out',args.out/'graphs',*common,*graph]]
        for command in commands:subprocess.run(list(map(str,command)),check=True)
        old=args.previous_b2/locus;direct=args.out/'graphs'/locus;imported=args.out/'import_audit'/locus
        assert source_signature(old)==source_signature(direct)==source_signature(imported), 'P/molecule membership changed'
        oldcuts,newcuts=cut_signature(old),cut_signature(direct)
        assert oldcuts==newcuts,'B2 core cuts or candidate pairs changed'
        assert digest(old/f'{locus}.structure.gfa')==digest(direct/f'{locus}.structure.gfa'), 'Direct graph changed'
        # Imported per-read annotations retain platform-dependent decomposition;
        # the authoritative B2 graph uses the pre-existing common HiFi profile.
        # Both routes must preserve the same raw sequences and membership.
        assert digest(direct/'paths.original.fa')==digest(imported/'paths.original.fa')
        summary=read(direct/'summary.json');validation=read(direct/'validation.json')
        reports.append(dict(locus=locus,raw_records=len(newcuts),paths=summary['paths'],molecules=summary['molecules'],
                            graph_stats=summary['graph_stats'],substituted_bases=summary['substituted_bases'],
                            validation=validation,previous_B2_cuts_identical=True,previous_B2_members_identical=True,
                            previous_B2_graph_identical=True,import_sequences_and_members_identical=True,
                            imported_graph_stats=read(imported/'summary.json')['graph_stats'],
                            import_and_direct_graph_identical=digest(direct/f'{locus}.structure.gfa')==digest(imported/f'{locus}.structure.gfa'),
                            extraction_status=dict(Counter(r.get('extraction_status',r['status']) for r in lines(direct/'records.jsonl'))),
                            core_status=dict(Counter(r.get('core_status',r['status']) for r in lines(direct/'records.jsonl')))))
        (args.out/'verification.json').write_text(json.dumps(reports,indent=2)+'\n')
        print('B2_FULL_PIPELINE_VERIFIED',locus,flush=True)

if __name__=='__main__':main()
