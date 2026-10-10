"""无需真实 BAM 的核心回归测试。"""

import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[1] / "scripts" / "vntrgraph.py"
spec = importlib.util.spec_from_file_location("vntrgraph", MODULE)
mod = importlib.util.module_from_spec(spec)
import sys
sys.modules["vntrgraph"] = mod
spec.loader.exec_module(mod)


class VntrGraphTests(unittest.TestCase):
    def test_catalog_and_coordinates(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[1]) as folder:
            path = Path(folder) / "loci.bed"
            path.write_text("chr1\t100\t112\tID=TEST;MOTIFS=ACGT\n", encoding="utf-8")
            self.assertEqual(mod.read_catalog(path)[0], mod.Locus("TEST", "chr1", 100, 112, ("ACGT",)))

    def test_reference_fai(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[1]) as folder:
            fasta = Path(folder) / "ref.fa"
            fasta.write_bytes(b">chr1\nACGTACGT\nTTTTAAAA\n")
            Path(str(fasta) + ".fai").write_text("chr1\t16\t6\t8\t9\n", encoding="utf-8")
            with mod.IndexedFastaReference(fasta) as reference:
                self.assertEqual(reference.fetch("chr1", 6, 12), "GTTTTT")

    def test_anchor_preserves_insertion(self):
        left = "CGATTCGATTCGATTCGATTCGATTCGATT"
        right = "TGGCATGGCATGGCATGGCATGGCATGGCA"
        sequence = left + "ACGTACATACAT" + right
        a = mod.flank_hit(left, sequence, 0, .9)
        b = mod.flank_hit(right, sequence, len(left) + 12, .9)
        self.assertEqual(sequence[a[1]:b[0]], "ACGTACATACAT")

    def test_two_alleles_reconstruct_exactly(self):
        reads = ([mod.Read(str(i), "S1", "HiFi", "ACGTACGTACGTACAT") for i in range(3)]
                 + [mod.Read(str(i), "S2", "HiFi", "ACGTACATACAT") for i in range(3)])
        alleles = mod.cluster_reads(reads, 3)
        motifs = ("ACGT", "ACAT")
        for allele in alleles:
            allele.tokens, explained = mod.decompose(allele.sequence, motifs)
            self.assertEqual(explained, 1.0)
        nodes, edges, paths = mod.build_graph(alleles, "TTGA", "GGCA")
        self.assertEqual(len(paths), 2)
        self.assertIn(("ACAT", "ACAT"), nodes)
        self.assertTrue(edges)

    def test_unknown_is_not_counted(self):
        tokens, explained = mod.decompose("ACGTTTTT", ("ACGT",))
        self.assertEqual("".join(text for _, text in tokens), "ACGTTTTT")
        self.assertLess(explained, 1.0)

    def test_same_length_variant_kept_separate(self):
        first = "ACGT" * 30
        second = first[:60] + "T" + first[61:]
        reads = ([mod.Read(f"a{i}", "S1", "HiFi", first) for i in range(3)]
                 + [mod.Read(f"b{i}", "S1", "HiFi", second) for i in range(3)])
        alleles = mod.cluster_reads(reads, 3)
        self.assertEqual({a.sequence for a in alleles}, {first, second})

    def test_cigar_projection(self):
        self.assertEqual(mod.reference_to_query(100, "10M8I10M", 110), 10)
        self.assertEqual(mod.reference_to_query(100, "10M8I10M", 115), 23)

    def test_motif_small_indel(self):
        tokens, explained = mod.decompose("ACGTTACGT", ("ACGT",))
        self.assertEqual("".join(text for _, text in tokens), "ACGTTACGT")
        self.assertEqual(explained, 1.0)

    def test_consensus_retains_multibase_insertion(self):
        reads = [mod.Read(str(i), "S1", "HiFi", "ACGTTTACGT") for i in range(3)]
        reads.append(mod.Read("error", "S1", "HiFi", "ACGTACGT"))
        self.assertEqual(mod.make_consensus(reads), "ACGTTTACGT")


if __name__ == "__main__":
    unittest.main()
