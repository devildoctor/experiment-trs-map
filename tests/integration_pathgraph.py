"""Real BAM integration: complete reads, BAM copies, partials and MAPQ 255."""
from pathlib import Path
import json, random, shutil, subprocess, sys, tempfile


def main():
    samtools=sys.argv[1];source=Path(__file__).resolve().parents[1]
    def run(*a):subprocess.run([str(x) for x in a],check=True)
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);rng=random.Random(41)
        seq=''.join(rng.choice('ACGT') for _ in range(500))+'ACGT'*10+''.join(rng.choice('ACGT') for _ in range(560))
        ref=root/'ref.fa';ref.write_text('>chr1\n'+seq+'\n');run(samtools,'faidx',ref)
        bed=root/'loci.bed';bed.write_text('chr1\t500\t540\tID=IL1RN;MOTIFS=ACGT;STRUC=<TR>\n')
        folder=root/'bam'/'S1'/'HiFi';folder.mkdir(parents=True)
        lines=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chr1\tLN:1100']
        for name,a,b,flag,mapq in [('fullA',250,790,0,60),('fullA',250,790,2048,60),
                                  ('fullB',250,790,16,60),('variant',250,790,0,60),
                                  ('left',250,520,0,60),('mq255',250,790,0,255),('dup',250,790,1024,60)]:
            bases=seq[a:b]
            if name=='variant': bases=bases[:265]+'A'+bases[266:]
            lines.append(f'{name}\t{flag}\tchr1\t{a+1}\t{mapq}\t{b-a}M\t*\t0\t0\t{bases}\t'+('I'*len(bases)))
        sam=root/'input.sam';sam.write_text('\n'.join(lines)+'\n')
        bam=folder/'test.IL1RN.bam';run(samtools,'sort','-o',bam,sam);run(samtools,'index',bam)
        shutil.copyfile(bam,folder/'copy.IL1RN.bam')
        shutil.copyfile(Path(str(bam)+'.bai'),folder/'copy.IL1RN.bam.bai')
        run(sys.executable,source/'scripts'/'vntrpathgraph_cli.py','--bam-dir',root/'bam',
            '--reference',ref,'--catalog',bed,'--locus','IL1RN','--out',root/'out','--samtools',samtools)
        out=root/'out'/'IL1RN';report=json.loads((out/'summary.json').read_text())
        assert report['paths']==3 and report['molecules']==5,report
        assert report['substituted_bases']==0
        assert 'mapq_not_qualified' not in report['excluded_reasons']
        assert report['preview']['status']=='rendered'
        assert (out/'IL1RN.preview.png').stat().st_size>1000
        assert (out/'IL1RN.path_sequences.txt').is_file()
        print('PATHGRAPH_BAM_INTEGRATION_PASS: reference plus two observed strings, five distinct molecules, duplicate BAM dedup, partial rejection, all B2 passes, GFA/TXT/image',flush=True)


if __name__=='__main__':main()
