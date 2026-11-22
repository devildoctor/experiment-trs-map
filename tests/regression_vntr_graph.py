"""End-to-end >6 bp VNTR regression against a source snapshot, on Slurm."""
from pathlib import Path
import argparse
import os
import random
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SAMTOOLS = os.environ.get("SAMTOOLS","/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools")


def run(*args):
    subprocess.run([str(a) for a in args],check=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline",type=Path,required=True)
    args = parser.parse_args()
    rng = random.Random(83)
    left, right = ["".join(rng.choices("ACGT",k=180)) for _ in range(2)]
    motif = "GAGGACCACCAGGAAG"
    repeat = motif*4
    sequence = left+repeat+right
    with tempfile.TemporaryDirectory(dir=ROOT) as temp:
        folder = Path(temp)
        fasta = folder/"ref.fa"
        fasta.write_text(">chrT\n"+sequence+"\n")
        run(SAMTOOLS,"faidx",fasta)
        catalog = folder/"vntr.bed"
        catalog.write_text(f"chrT\t180\t244\tID=LONGVNTR;MOTIFS={motif}\n")
        bamdir = folder/"bam"/"S1"/"HiFi"
        bamdir.mkdir(parents=True)
        sam = folder/"input.sam"
        sam.write_text(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:{len(sequence)}\n"+
                       "".join(f"r{i}\t0\tchrT\t1\t60\t{len(sequence)}M\t*\t0\t0\t{sequence}\t{'I'*len(sequence)}\n"
                               for i in range(3)))
        bam = bamdir/"synthetic.S1.LONGVNTR.bam"
        run(SAMTOOLS,"sort","-o",bam,sam)
        run(SAMTOOLS,"index",bam)
        for source,label in ((args.baseline,"before"),(ROOT,"after")):
            run(sys.executable,source/"scripts"/"vntrgraph.py","run","--bam-dir",folder/"bam",
                "--reference",fasta,"--catalog",catalog,"--out",folder/label,
                "--samtools",SAMTOOLS,"--anchor-length","150")
        old, new = folder/"before"/"LONGVNTR", folder/"after"/"LONGVNTR"
        assert (new/"LONGVNTR.gfa").is_file()
        assert not (new/"LONGVNTR.detail.gfa").exists()
        assert not (new/"str_paths.json").exists()
        checked = 0
        for path in old.iterdir():
            if path.is_file() and path.name != "parameters.tsv":
                assert path.read_bytes() == (new/path.name).read_bytes(), path.name
                checked += 1
        assert {p.name for p in old.iterdir()} == {p.name for p in new.iterdir()}
        print(f"VNTR_16BP_REGRESSION_PASS: {checked} files byte-identical; parameters.tsv records distinct source/output paths")


if __name__ == "__main__":
    main()
