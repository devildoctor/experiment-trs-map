"""End-to-end real BAM format test; executed inside each Slurm job."""
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main():
    samtools=sys.argv[1];source=Path(__file__).resolve().parents[1]
    def run(*a):subprocess.run([str(v) for v in a],check=True)
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);rng=random.Random(41)
        seq=''.join(rng.choice('ACGT') for _ in range(500))+'ACGT'*10+''.join(rng.choice('ACGT') for _ in range(560))
        ref=root/'ref.fa';ref.write_text('>chr1\n'+seq+'\n',encoding='utf-8');run(samtools,'faidx',ref)
        bed=root/'loci.bed';bed.write_text('chr1\t500\t540\tID=IL1RN;MOTIFS=ACGT;STRUC=<TR>\n',encoding='utf-8')
        folder=root/'bam'/'S1'/'HiFi';folder.mkdir(parents=True)
        lines=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chr1\tLN:1100']
        for name,a,b,flag,mapq in [('fullA',250,790,0,60),('fullA',250,790,2048,60),('fullB',250,790,16,60),('left',250,520,0,60),('right',520,790,0,60),('internal',510,530,1024,0),('flank',50,100,0,60)]:
            bases=seq[a:b];lines.append(f'{name}\t{flag}\tchr1\t{a+1}\t{mapq}\t{b-a}M\t*\t0\t0\t{bases}\t'+('I'*len(bases)))
        lines.append('missing\t0\tchr1\t501\t60\t40M\t*\t0\t0\t*\t*')
        sam=root/'test.sam';sam.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        bam=folder/'test.IL1RN.bam';run(samtools,'sort','-o',bam,sam)
        run(samtools,'index',bam)
        run(sys.executable,source/'scripts'/'vntrgraph.py','--structure-view','--bam-dir',root/'bam','--reference',ref,
            '--catalog',bed,'--out',root/'out','--locus','IL1RN','--samtools',samtools,'--workers','2')
        out=root/'out'/'IL1RN';d=json.loads((out/'IL1RN.structures.json').read_text(encoding='utf-8'))
        rows={r['read_id']:r for r in d['rows']}
        assert len(rows)==7 and d['audit']['alignment_records']==8
        assert rows['left']['fragments'][0]['anchor']=='left_only'
        assert rows['right']['fragments'][0]['anchor']=='right_only'
        assert rows['internal']['fragments'][0]['anchor']=='cigar_only'
        assert rows['flank']['status']=='flank_only' and rows['missing']['status']=='sequence_unavailable'
        assert len(rows['fullA']['fragments'])==2
        n=next(n for n in d['graph']['structure_nodes'] if n['id']=='M1')
        assert n['support']==5,(n,rows.keys())
        assert n['occurrences']>n['support']
        assert d['audit']['self_loop_edges']>0
        run(sys.executable,source/'tests'/'verify_structure_results.py',out)
        print('STRUCTURE_BAM_INTEGRATION_PASS: all 7 reads/8 records; single anchors; shared nodes; support dedup; exact reconstruction',flush=True)


if __name__=='__main__':main()
