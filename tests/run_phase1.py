"""Reproducible first-stage acceptance; retain inputs, outputs and logs."""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import runpy
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'tests')]
from test_phase1 import Phase1Tests,reads,structure_rows
from vntrcore import Allele,decompose,motif_dictionary,build_graph
from vntrgraph import Locus,write_svg,write_gfa,write_tsv,source_hashes
from vntrstructure import seed_dictionary,learn_variants,annotate_rows,emit
from vntrmotifs import compound_annotations


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--samtools',required=True);args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    results=[]
    for name in unittest.defaultTestLoader.getTestCaseNames(Phase1Tests):
        result=unittest.TestResult();Phase1Tests(name).run(result)
        results.append(dict(case=name,status='PASS' if result.wasSuccessful() else 'FAIL',
                            details=[detail for _,detail in result.errors+result.failures]))
    (args.out/'acceptance.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    write_tsv(args.out/'acceptance.tsv',['case','status','details'],results)
    if any(r['status']!='PASS' for r in results):raise AssertionError(results)
    a='AACCGTACGATGAACC';b='AACCGTTCGATGAACC'
    examples=[('PURE','ACGT'*30,('ACGTACGT','ACGTACGTACGT')),
              ('MIXED','ACGT'*25+'GATAA'*20,('ACGT','GATAA')),
              ('INTERRUPTION',a+a+b+a+a,(a,)),
              ('AMBIGUOUS','ACGT',('ACAT','ACCT'))]
    summaries=[]
    for name,sequence,catalog in examples:
        destination=args.out/'examples'/name;destination.mkdir(parents=True)
        allele=Allele(name+'_A1','S',sequence,read_ids=tuple(r.molecule_key for r in reads(sequence)),support={'HiFi':3})
        records=motif_dictionary(sequence,[allele],reads(sequence),catalog,3)
        templates=[m['motif'] for m in records if m['level']=='base' and m['usable']]
        allele.tokens,_,allele.annotation=decompose(sequence,templates,return_details=True)
        allele.compound_details=compound_annotations(allele.tokens,[m['motif'] for m in records if m['source']=='catalog' and m['level']=='compound'])
        allele.uncertain=allele.annotation['annotation_ambiguity'] or any(not m['validated'] for m in records if m['usable'])
        if allele.uncertain:allele.flags=['annotation_or_seed_uncertainty']
        write_svg(destination/f'{name}.svg',name,[allele])
        write_gfa(destination/f'{name}.gfa',*build_graph([allele],'G'*25,'C'*25))
        write_tsv(destination/'motifs.tsv',list(records[0]),records)
        (destination/'annotation_alternatives.json').write_text(json.dumps({allele.name:allele.annotation},indent=2),encoding='utf-8')
        (destination/'input.json').write_text(json.dumps(dict(sequence=sequence,catalog=catalog,synthetic=True)),encoding='utf-8')
        (destination/'compound_blocks.json').write_text(json.dumps(allele.compound_details,indent=2),encoding='utf-8')
        structure=destination/'structure';structure.mkdir()
        dic=seed_dictionary(catalog);rows=structure_rows(sequence,dic)
        dic,candidates=learn_variants(rows,dic);annotate_rows(rows,dic,1)
        emit(structure,Locus(name,'chrSynthetic',0,len(sequence),catalog),rows,dic,
             dict(read_rows=3,rows_with_repeat_sequence=3,motif_candidates=candidates))
        metrics=allele.annotation['metrics'];structure_counts=[r['fragments'][0]['base_metrics']['full_units'] for r in rows]
        if name=='PURE':assert metrics['full_units']==30 and structure_counts==[30]*3
        if name=='MIXED':assert metrics['full_units']==45 and structure_counts==[45]*3
        if name=='INTERRUPTION':assert metrics['full_units']==5 and structure_counts==[5]*3
        if name=='AMBIGUOUS':assert allele.annotation['annotation_ambiguity']
        summaries.append(dict(case=name,base_metrics=metrics,structure_unit_counts=structure_counts,
                              compound_blocks=allele.compound_details,annotation_score_gap=allele.annotation['score_gap']))
    (args.out/'examples.json').write_text(json.dumps(summaries,indent=2),encoding='utf-8')
    # Run the three existing CLI integration programs, preserving their temporary directories.
    integration=args.out/'integration';integration.mkdir()
    for script in ('integration_synthetic.py','integration_all_reads.py','integration_structure.py'):
        output=integration/Path(script).stem;output.mkdir()
        class RetainedDirectory:
            def __init__(self,*a,**kw):self.name=tempfile.mkdtemp(prefix='inputs-and-results-',dir=output)
            def __enter__(self):return self.name
            def __exit__(self,*a):pass
        original=tempfile.TemporaryDirectory
        tempfile.TemporaryDirectory=RetainedDirectory
        sys.argv=[str(ROOT/'tests'/script),args.samtools];os.environ['SAMTOOLS']=args.samtools
        try:runpy.run_path(sys.argv[0],run_name='__main__')
        finally:tempfile.TemporaryDirectory=original
    (args.out/'run-manifest.json').write_text(json.dumps(dict(status='PASS',acceptance_cases=len(results),
        integration_entrypoints=3,source_hashes=source_hashes(),python=sys.version,samtools=args.samtools,
        evidence='synthetic sequences and synthetic BAM inputs; not validation against biological truth'),indent=2),encoding='utf-8')
    print(f'PHASE1_PASS: {len(results)} acceptance cases; 4 saved examples; 3 CLI integration entrypoints',flush=True)


if __name__=='__main__':main()
