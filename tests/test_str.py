"""STR contract tests; run together with the existing suite on Slurm."""
import copy
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
import vntrstr as s
from vntrgraph import read_catalog


class STRTests(unittest.TestCase):
    def test_2_to_6_bp(self):
        for motif in ("AC", "CAG", "ACGT", "AACGT", "AACCGT"):
            a = s.annotate(motif*12, f"({motif})*")
            self.assertEqual(a["expression"], f"({motif})12")
            self.assertTrue(a["model_match"])
            self.assertEqual(s.reconstruct(a["blocks"]), motif*12)

    def test_fixed_interruption_and_independent_counts(self):
        a = s.annotate("CAG"*12+"TTGACC"+"CTG"*7, "(CAG)*TTGACC(CTG)*")
        self.assertEqual(a["expression"], "(CAG)12–TTGACC–(CTG)7")
        self.assertTrue(a["model_match"])
        self.assertEqual([b["count"] for b in a["blocks"] if b["kind"]=="repeat"], [12,7])

    def test_error_and_tail_not_counted(self):
        a = s.annotate("CAG"*5+"CAA"+"CAG"*7+"CA", "(CAG)*")
        self.assertEqual(a["expression"], "(CAG)5–CAA–(CAG)7–CA")
        self.assertEqual(a["literal_bases"], 5)
        self.assertFalse(a["model_match"])

    def test_single_unit_and_N_preserved(self):
        a = s.annotate("CAG"+"N"+"CTG"*4, "(CAG)*(CTG)*")
        self.assertEqual(a["expression"], "CAGN–(CTG)4")

    def test_reverse_complements_order_and_literals(self):
        seq="CAG"*3+"TTGACC"+"CTG"*7
        a=s.annotate(seq,"(CAG)*TTGACC(CTG)*")
        self.assertEqual(a["reverse_complement_expression"],"(CAG)7–GGTCAA–(CTG)3")
        self.assertEqual(s.reverse_blocks(a["reverse_complement_blocks"]),a["blocks"])

    def test_equal_cost_ambiguity(self):
        a=s.annotate("CAGCAGCAGC","(CAG)*(AGC)*")
        self.assertTrue(a["annotation_ambiguity"])
        self.assertEqual(a["alternative"]["cost"],a["cost"])
        self.assertEqual(s.reconstruct(a["alternative"]["blocks"]),"CAGCAGCAGC")

    def test_fixed_phase_not_rotated(self):
        a=s.annotate("GCT"*5,"(CTG)*")
        self.assertEqual(a["expression"],"G–(CTG)4–CT")
        self.assertFalse(a["model_match"])

    def test_compound_to_primitive_metadata(self):
        a=s.annotate("CAG"*6,"(CAGCAG)*")
        self.assertEqual(a["expression"],"(CAG)6")
        self.assertEqual(a["model"][0]["repeat_factor"],2)
        self.assertEqual(a["model"][0]["catalog_motif"],"CAGCAG")

    def test_scope_and_syntax(self):
        for model in ("(A)*","(AAAAAA)*","(ACGTACG)*","(CAG)+","(CAG){3}", "CAG", "(CAG)*?"):
            with self.assertRaises(ValueError):
                s.parse_structure(model)
        for motif in ("A","ACGTACG","CCTGGTGGTCCTCCTT","CAGCAGCAG"):
            self.assertEqual(s.validate_config("",(motif,)),"")
            with self.assertRaises(ValueError):
                s.validate_config("(CAG)*",(motif,))

    def test_catalog_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"catalog.bed"
            p.write_text("chr1\t1\t10\tID=VNTR;MOTIFS=CCTGGTGGTCCTCCTT\n",encoding="utf-8")
            self.assertEqual(read_catalog(p)[0].str_structure,"")
            p.write_text("chr1\t1\t10\tID=STR;MOTIFS=CAG;STR_STRUCTURE=(CAG)*\n",encoding="utf-8")
            self.assertEqual(read_catalog(p)[0].str_structure,"(CAG)*")

    def test_fragment_and_candidate_agree_without_mutating_sequence(self):
        sequence="CTG"*3+"CAT"+"CTG"*5
        rows=[dict(id="r",read_id="read",sample="s",platform="ONT",fragments=[
            dict(id="f",sequence=sequence,coverage="left_partial",anchor="left_only",qstart=8,qend=35)])]
        original=copy.deepcopy(rows)
        record=s.annotate_rows(rows,"(CTG)*")[0]
        self.assertEqual(record["coverage"],"left_partial")
        self.assertEqual(record["annotation"],s.path_record("c",sequence,"(CTG)*")["annotation"])
        del rows[0]["fragments"][0]["str_annotation"]
        self.assertEqual(rows,original)

    def test_lossless_mixed_random_paths_and_alternatives(self):
        rng=random.Random(1729)
        for _ in range(80):
            seq="CAG"*rng.randrange(8)+"".join(rng.choices("ACGTN",k=rng.randrange(8)))+"CTG"*rng.randrange(8)
            a=s.annotate(seq,"(CAG)*(CTG)*")
            self.assertEqual(s.reconstruct(a["blocks"]),seq)
            if a["alternative"]:
                self.assertEqual(s.reconstruct(a["alternative"]["blocks"]),seq)
            self.assertEqual(s.reconstruct(a["reverse_complement_blocks"]),
                             seq.translate(str.maketrans("ACGTN","TGCAN"))[::-1])

    def test_candidate_figure_keeps_support_and_uncertainty(self):
        record=s.path_record("c","CTG"*8,"(CTG)*",kind="candidate",sample="S",
                             read_support=53,platform_support={"HiFi":40,"ONT":13},
                             call_status="candidate_ploidy_conflict",ambiguous_molecules=4,
                             uncertain=True,flags=["annotation_ambiguity"])
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"figure.svg"
            s.write_paths_svg(p,"STR",[record])
            text=p.read_text(encoding="utf-8")
            self.assertIn("candidate_ploidy_conflict",text)
            self.assertIn("molecules=53",text)
            self.assertIn("ambiguous_molecules=4",text)
            self.assertIn("UNCERTAIN",text)


if __name__=="__main__":
    unittest.main()
