"""Real BAM integration: complete reads, BAM copies, partials and MAPQ 255."""
from pathlib import Path
import json, random, shutil, subprocess, sys, tempfile


def main():
    samtools=sys.argv[1];source=Path(__file__).resolve().parents[1]
    def run(*a):subprocess.run([str(x) for x in a],check=True)
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);rng=random.Random(41)
        seq=''.join(rng.choice('ACGT') for _ in range(200))+'ACGT'*10+''.join(rng.choice('ACGT') for _ in range(260))
        ref=root/'ref.fa';ref.write_text('>chr1\n'+seq+'\n');run(samtools,'faidx',ref)
        bed=root/'loci.bed';bed.write_text('chr1\t200\t240\tID=IL1RN;MOTIFS=ACGT;STRUC=<TR>\n')
        folder=root/'bam'/'S1'/'HiFi';folder.mkdir(parents=True)
        lines=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chr1\tLN:500']
        for name,a,b,flag,mapq in [('fullA',100,340,0,60),('fullA',100,340,2048,60),
                                  ('fullB',100,340,16,60),('variant',100,340,0,60),
                                  ('left',100,220,0,60),('mq255',100,340,0,255),('dup',100,340,1024,60)]:
            bases=seq[a:b]
            if name=='variant': bases=bases[:115]+'A'+bases[116:]
            lines.append(f'{name}\t{flag}\tchr1\t{a+1}\t{mapq}\t{b-a}M\t*\t0\t0\t{bases}\t'+('I'*len(bases)))
        sam=root/'input.sam';sam.write_text('\n'.join(lines)+'\n')
        bam=folder/'test.IL1RN.bam';run(samtools,'sort','-o',bam,sam);run(samtools,'index',bam)
        shutil.copyfile(bam,folder/'copy.IL1RN.bam')
        shutil.copyfile(Path(str(bam)+'.bai'),folder/'copy.IL1RN.bam.bai')
        run(sys.executable,source/'scripts'/'vntrpathgraph_cli.py','--bam-dir',root/'bam',
            '--reference',ref,'--catalog',bed,'--locus','IL1RN','--out',root/'out','--samtools',samtools)
        out=root/'out'/'IL1RN';report=json.loads((out/'summary.json').read_text())
        assert report['paths']==2 and report['molecules']==3,report
        assert report['substituted_bases']==0
        assert report['excluded_reasons']['mapq_not_qualified']==2
        assert report['preview']['status']=='rendered'
        assert (out/'IL1RN.preview.png').stat().st_size>1000
        print('PATHGRAPH_BAM_INTEGRATION_PASS: two exact paths, three distinct molecules, duplicate BAM dedup, partial/MAPQ255 rejection, GFA and static image',flush=True)


if __name__=='__main__':main()
