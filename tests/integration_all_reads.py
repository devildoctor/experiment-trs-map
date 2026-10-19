"""Small real-BAM integration check, run on a Slurm compute node."""
import csv
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def run(*args):
    subprocess.run([str(v) for v in args],check=True)


def main():
    samtools=sys.argv[1]
    source=Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory() as name:
        root=Path(name)
        rng=random.Random(41)
        sequence=''.join(rng.choice('ACGT') for _ in range(200))+'ACGT'*10+''.join(rng.choice('ACGT') for _ in range(260))
        reference=root/'reference.fa'
        reference.write_text('>chr1\n'+sequence+'\n',encoding='utf-8')
        run(samtools,'faidx',reference)
        catalog=root/'loci.bed'
        catalog.write_text('chr1\t200\t240\tID=IL1RN;MOTIFS=ACGT\n',encoding='utf-8')
        folder=root/'bamdata'/'S1'/'HiFi'
        folder.mkdir(parents=True)
        sam=root/'input.sam'
        lines=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chr1\tLN:500']
        def record(read_id,start,end,flag=0,mapq=60):
            bases=sequence[start:end]
            lines.append(f'{read_id}\t{flag}\tchr1\t{start+1}\t{mapq}\t{end-start}M\t*\t0\t0\t{bases}\t'+('I'*len(bases)))
        record('full_forward',100,340)
        record('full_reverse',100,340,16)
        record('single_left',100,220)
        record('single_right',220,340)
        record('internal_low_mapq_duplicate',210,230,1024,0)
        record('flank_only',50,100)
        record('full_forward',100,340,2048)
        lines.append('missing_seq\t0\tchr1\t201\t60\t40M\t*\t0\t0\t*\t*')
        sam.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        bam=folder/'sample.IL1RN.bam'
        run(samtools,'sort','-o',bam,sam)
        run(samtools,'index',bam)
        run(sys.executable,source/'scripts'/'vntrgraph.py','--read-view','--bam-dir',root/'bamdata',
            '--reference',reference,'--catalog',catalog,'--out',root/'out','--locus','IL1RN',
            '--samtools',samtools,'--workers','2','--anchor-length','100')
        out=root/'out'/'IL1RN'
        data=json.loads((out/'IL1RN.reads.json').read_text(encoding='utf-8'))
        rows={r['read_id']:r for r in data['rows']}
        assert len(rows)==7,rows.keys()
        assert data['audit']['alignment_records']==8
        assert len(rows['full_forward']['fragments'])==2
        assert rows['single_left']['fragments'][0]['anchor']=='left_only'
        assert rows['single_right']['fragments'][0]['anchor']=='right_only'
        assert rows['internal_low_mapq_duplicate']['fragments'][0]['anchor']=='cigar_only'
        assert rows['flank_only']['status']=='flank_only'
        assert rows['missing_seq']['status']=='sequence_unavailable'
        assert rows['full_forward']['fragments'][0]['sequence']==rows['full_reverse']['fragments'][0]['sequence']
        nodes,paths={},{}
        for line in (out/'IL1RN.reads.gfa').read_text(encoding='utf-8').splitlines():
            fields=line.split('\t')
            if fields[0]=='S':nodes[fields[1]]=fields[2]
            if fields[0]=='P':paths[fields[1]]=fields[2].split(',')
        for row in rows.values():
            for f in row['fragments']:
                assert ''.join(nodes[n[:-1]] for n in paths[f['id']])==f['sequence']
        print('ALL_READS_INTEGRATION_PASS: 7/7 reads, 8/8 records, single anchors and partial paths verified',flush=True)


if __name__=='__main__':
    main()
