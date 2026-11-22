"""Run all B2-passing sequences through reference-anchored lossless graphs."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--previous-b2',type=Path,required=True)
    parser.add_argument('--samtools',required=True)
    parser.add_argument('--locus',action='append',choices=['ATXN2','D1S80'])
    args=parser.parse_args();source=Path(__file__).resolve().parents[1]
    args.out.mkdir(parents=True,exist_ok=False)
    reference=args.project/'reference/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    catalog=source/'catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed'
    reports=[]
    for locus in args.locus or ['ATXN2','D1S80']:
        command=[sys.executable,source/'scripts/vntrpathgraph_cli.py','--bam-dir',args.project/'bamdata_by_locus',
                 '--reference',reference,'--catalog',catalog,'--locus',locus,'--samtools',args.samtools,
                 '--out',args.out/'graphs','--graph-mode','structure',
                 '--locus-model',source/f'catalogs/{locus}.structure.v1.json',
                 '--structure-policy',source/f'catalogs/{locus}.structure-policy.v1.json',
                 '--max-loss-bp','0','--preview-paths','0',
                 '--previous-members',args.previous_b2/locus/'source_members.json']
        subprocess.run(list(map(str,command)),check=True)
        folder=args.out/'graphs'/locus
        subprocess.run([sys.executable,str(source/'tests/verify_all_flank_results.py'),str(folder),
                        str(args.previous_b2/locus),str(reference),str(catalog)],check=True)
        reports.append(json.loads((folder/'verification.json').read_text()))
        (args.out/'verification.json').write_text(json.dumps(reports,indent=2)+'\n')
        print('B2_ALL_PASS_PIPELINE_VERIFIED',locus,flush=True)


if __name__=='__main__':main()
