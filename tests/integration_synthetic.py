"""从合成 FASTA/SAM/BAM 到 GFA/SVG 的独立闭环测试。"""

from pathlib import Path
import random
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SAMTOOLS = Path("/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools")


def run(*command):
    subprocess.run([str(item) for item in command], check=True, capture_output=True, text=True)


def main():
    random.seed(41)
    left = "".join(random.choices("ACGT", k=150))
    right = "".join(random.choices("ACGT", k=150))
    repeat1 = "ACGTACGTACGTACAT"
    repeat2 = "ACGTACATACAT"
    with tempfile.TemporaryDirectory(dir=ROOT) as folder:
        folder = Path(folder)
        fasta = folder / "ref.fa"
        fasta.write_text(f">chrT\n{left + repeat1 + right}\n", encoding="ascii")
        run(SAMTOOLS, "faidx", fasta)
        catalog = folder / "candidate.bed"
        catalog.write_text(f"chrT\t150\t166\tID=TEST;MOTIFS=ACGT,ACAT\n", encoding="ascii")
        bamdir = folder / "bamdata"
        for sample, repeat, cigar in (("S1", repeat1, "316M"), ("S2", repeat2, "162M4D150M")):
            directory = bamdir / sample / "HiFi"
            directory.mkdir(parents=True)
            sam = directory / "test.sam"
            header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrT\tLN:316\n"
            records = []
            for number in range(3):
                seq = left + repeat + right
                records.append(f"r{number}\t0\tchrT\t1\t60\t{cigar}\t*\t0\t0\t{seq}\t{'I' * len(seq)}\n")
            sam.write_text(header + "".join(records), encoding="ascii")
            bam = directory / f"synthetic.{sample}.TEST.bam"
            run(SAMTOOLS, "sort", "-o", bam, sam)
            run(SAMTOOLS, "index", bam)
        output = folder / "results"
        run(sys.executable, ROOT / "scripts" / "vntrgraph.py", "run", "--bam-dir", bamdir,
            "--reference", fasta, "--catalog", catalog, "--out", output,
            "--samtools", SAMTOOLS, "--anchor-length", "150")
        result = output / "TEST"
        assert (result / "TEST.svg").is_file()
        gfa = (result / "TEST.gfa").read_text(encoding="utf-8")
        nodes = {}
        paths = []
        for line in gfa.splitlines():
            parts = line.split("\t")
            if parts[0] == "S":
                nodes[parts[1]] = parts[2]
            if parts[0] == "P":
                paths.append("".join(nodes[x[:-1]] for x in parts[2].split(",")))
        assert len(paths) == 2, paths
        assert set(paths) == {left[-25:] + repeat1 + right[:25], left[-25:] + repeat2 + right[:25]}
    print("synthetic integration OK")


if __name__ == "__main__":
    main()
