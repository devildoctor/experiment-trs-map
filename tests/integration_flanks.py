"""Real BAM/CLI round trip for new extraction, both entry points and import guard."""
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from vntrgraph import Locus


def main():
    samtools=sys.argv[1];source=Path(__file__).resolve().parents[1]
    def run(*args):subprocess.run([str(a) for a in args],check=True)
    with tempfile.TemporaryDirectory(prefix='flank-integration-') as temporary:
        root=Path(temporary);rng=random.Random(81)
        dna=''.join(rng.choice('ACGT') for _ in range(500))+'CTG'*10+''.join(rng.choice('ACGT') for _ in range(570))
        reference=root/'reference.fa';reference.write_text('>chr1\n'+dna+'\n');run(samtools,'faidx',reference)
        catalog=root/'loci.bed';catalog.write_text('chr1\t500\t530\tID=ATXN2;MOTIFS=CTG;STRUC=<TR>;BOUNDARY=confirmed\n')
        locus=Locus('ATXN2','chr1',500,530,('CTG',))
        structure=json.loads((source/'catalogs/ATXN2.structure.v1.json').read_text())
        structure.update(chrom='chr1',start=500,end=530)
        structure_path=root/'structure.json';structure_path.write_text(json.dumps(structure))
        bam_folder=root/'bam/LCL5/HiFi';bam_folder.mkdir(parents=True)
        lines=['@HD\tVN:1.6\tSO:unsorted',f'@SQ\tSN:chr1\tLN:{len(dna)}']
        for name,a,b,flag,mq in [('full',0,1100,0,60),('full',0,1100,2048,60),('reverse',0,1100,16,60),
             ('core_only',250,780,0,60),('partial',250,515,0,60),('duplicate',0,1100,1024,60),('low',0,1100,0,0)]:
            sequence=dna[a:b]
            lines.append(f'{name}\t{flag}\tchr1\t{a+1}\t{mq}\t{b-a}=\t*\t0\t0\t{sequence}\t'+('I'*len(sequence)))
        sam=root/'test.sam';sam.write_text('\n'.join(lines)+'\n')
        bam=bam_folder/'test.ATXN2.bam';run(samtools,'sort','-o',bam,sam);run(samtools,'index',bam)
        common=['--reference',reference,'--catalog',catalog,'--locus','ATXN2','--samtools',samtools]
        run(sys.executable,source/'scripts/vntrpathgraph_cli.py','--bam-dir',root/'bam','--out',root/'graph',*common)
        graph=root/'graph/ATXN2';summary=json.loads((graph/'summary.json').read_text())
        assert (summary['molecules'],summary['paths'])==(3,1),summary
        audit=[json.loads(x) for x in (graph/'records.jsonl').read_text().splitlines()]
        assert len(audit)==7
        assert next(r for r in audit if r['read_id']=='core_only')['extraction_status']=='core_only'
        for r in audit:
            if r.get('core_status')=='success':assert r['sequence']=='CTG'*10
        assert (graph/'extraction_audit.jsonl').exists()
        assert not (graph/'contexts.html').exists()
        assert all('diagnostic_100bp' not in r and 'analysis_anchor_evidence' not in r for r in audit)
        run(sys.executable,source/'scripts/vntrstructure.py','--bam-dir',root/'bam','--out',root/'structure',
            '--locus-model',structure_path,'--workers','1',*common)
        payload_path=root/'structure/ATXN2/ATXN2.structures.json'
        payload=json.loads(payload_path.read_text())
        assert payload['audit']['boundary_model_sha256']==summary['metadata']['boundary_model_sha256']
        run(sys.executable,source/'scripts/vntrpathgraph_cli.py','--structures-json',payload_path,'--out',root/'imported',*common)
        imported=json.loads((root/'imported/ATXN2/source_members.json').read_text())
        assert {r['sha256'] for r in imported}=={r['sha256'] for r in json.loads((graph/'source_members.json').read_text())}
        payload['audit']['boundary_model_sha256']='wrong'
        payload_path.write_text(json.dumps(payload))
        bad=subprocess.run([str(a) for a in [sys.executable,source/'scripts/vntrpathgraph_cli.py','--structures-json',payload_path,
                         '--out',root/'rejected',*common]],capture_output=True,text=True)
        assert bad.returncode and 'Imported flank boundary model mismatch' in bad.stderr
    print('FLANK_BAM_BOTH_CLI_AND_IMPORT_PASS')


if __name__=='__main__':main()
