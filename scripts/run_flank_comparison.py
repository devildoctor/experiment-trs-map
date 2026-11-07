"""Run the frozen flank comparison against a separately replayed B0 baseline."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time

from vntrgraph import read_catalog, IndexedFastaReference
from vntrreads import write_table
from vntrflankbounds import load_snapshot, build_profile, discover, base_model, FETCH


def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,required=True)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--locus',choices=['ATXN2','D1S80'],required=True)
    args=parser.parse_args()
    source=Path(__file__).resolve().parents[1]
    catalog=source/'catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed'
    reference_path=args.project/'reference/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    locus=next(l for l in read_catalog(catalog) if l.name==args.locus)
    baseline=args.run/'comparison/B0'/locus.name
    records_path=baseline/'records.jsonl'
    records=load_snapshot(records_path)
    discovery=args.run/'models'/locus.name
    discovery.mkdir(parents=True,exist_ok=False)
    with IndexedFastaReference(reference_path) as reference:
        fl,fr=FETCH[locus.name];dna=reference.fetch(locus.chrom,fl,fr).upper()
    profile=build_profile(records,[fl,fr],dna)
    models,audits=discover(locus,profile)
    models={'B1':base_model(locus,[fl,fr],dna,100),'B2':base_model(locus,[fl,fr],dna),**models}
    profile_rows=[]
    for i in range(fr-fl):
        profile_rows.append(dict(position=fl+i,reference_base=dna[i],signal=float(profile['values'][i]),bin=int(profile['bins'][i]),
            coverage={s:int(v[i]) for s,v in profile['coverage'].items()},alts={s:float(v[i]) for s,v in profile['alts'].items()}))
    write_table(discovery/'profile.tsv',profile_rows,['position','reference_base','signal','bin','coverage','alts'])
    write_table(discovery/'boundary_candidates.tsv',audits,['group','side','distance','boundary','scores','anchor_status','anchor_counts','status'])
    hash_before=sha256(records_path.read_bytes()).hexdigest()
    for group,model in models.items():
        model['input_records_sha256']=hash_before
        model['algorithm_source_sha256']={name:sha256((source/'scripts'/name).read_bytes()).hexdigest() for name in ['vntrflankalign.py','vntrflankbounds.py']}
        dump(discovery/f'{group}.json',model)
        print(json.dumps(dict(stage='model',locus=locus.name,group=group,analysis_interval=model['analysis_interval'])),flush=True)
    timings={}
    for group in models:
        started=time.time()
        command=[sys.executable,str(source/'scripts/vntrpathgraph_cli.py'),
            '--bam-dir',str(args.project/'bamdata_by_locus'),'--reference',str(reference_path),
            '--catalog',str(catalog),'--out',str(args.run/'comparison'/group),'--locus',locus.name,
            '--samtools','/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools',
            '--preview-paths','10','--graph-mode','structure',
            '--locus-model',str(source/f'catalogs/{locus.name}.structure.v1.json'),
            '--structure-policy',str(source/f'catalogs/{locus.name}.structure-policy.v1.json'),
            '--max-loss-bp','0','--flank-boundary-model',str(discovery/f'{group}.json')]
        subprocess.run(command,check=True)
        timings[group]=time.time()-started
        dump(discovery/'timings.json',timings)
    assert sha256(records_path.read_bytes()).hexdigest()==hash_before
    print('COMPARISON_COMPLETE',locus.name,flush=True)


if __name__=='__main__':main()
