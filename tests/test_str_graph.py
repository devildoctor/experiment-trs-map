"""Graph semantics tests; execute on Slurm with the full existing suite."""
from pathlib import Path
import json
import random
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
import vntrstr as s
import vntrstrgraph as g


def record(name, sequence, model="(CTG)*"):
    return s.path_record(name, sequence, model, kind="candidate", uncertain=True)


class STRGraphTests(unittest.TestCase):
    def test_real_example_three_separate_contexts(self):
        seq = "CTG"*8+"TTG"+"CTG"*4+"TTG"+"CTG"*8
        graph = g.build_graph("ATXN2", [record("LCL5_C13", seq)], "AA", "CC")
        p = graph["paths"][0]
        self.assertEqual([r["count"] for r in p["runs"]], [1,8,1,4,1,8,1])
        self.assertEqual(len(graph["nodes"]), 7)
        self.assertEqual(len(graph["edges"]), 9)
        self.assertEqual(len({r["node"] for r in p["runs"]}), 7)
        self.assertEqual("".join(graph["nodes"][n]["sequence"] for n in p["walk"]), "AA"+seq+"CC")

    def test_counts_belong_to_paths(self):
        graph = g.build_graph("TEST", [record("A","CTG"*8), record("B","CTG"*12)], "AA","CC")
        self.assertEqual(len(graph["nodes"]), 3)
        self.assertEqual([p["runs"][1]["count"] for p in graph["paths"]], [8,12])
        self.assertEqual(graph["paths"][0]["runs"][1]["node"], graph["paths"][1]["runs"][1]["node"])
        self.assertNotIn("count", graph["nodes"]["N1"])

    def test_shared_suffix_and_distinct_literals(self):
        graph = g.build_graph("TEST", [record("A","CTG"*3+"TTG"+"CTG"*4),
                                       record("B","CTG"*5+"TCG"+"CTG"*6)], "AA","CC")
        paths = graph["paths"]
        self.assertEqual(paths[0]["runs"][-2]["node"], paths[1]["runs"][-2]["node"])
        self.assertNotEqual(paths[0]["runs"][2]["node"], paths[1]["runs"][2]["node"])
        self.assertNotEqual(paths[0]["runs"][1]["node"], paths[0]["runs"][-2]["node"])

    def test_molecule_support_not_multiplied_by_visits(self):
        graph = g.build_graph("TEST", [record("A","CTG"*8), record("B","CTG"*12)],
                              "AA","CC", {"A":["m1","m2"], "B":["m2","m3"]})
        self.assertEqual(graph["nodes"]["N1"]["molecule_support"], 3)
        loop = next(e for e in graph["edges"] if e["source"]==e["target"])
        self.assertEqual(loop["molecule_support"], 3)

    def test_deterministic_input_order(self):
        records = [record("A","CTG"*3+"TTG"+"CTG"*4), record("B","CTG"*6),
                   record("C","CTG"*5+"ATG"+"CTG"*8)]
        self.assertEqual(g.build_graph("TEST", records, "AA","CC"),
                         g.build_graph("TEST", list(reversed(records)), "AA","CC"))

    def test_unknown_tail_and_empty_internal(self):
        records = [record("A","CTG"*2+"NCT"), record("B","")]
        graph = g.build_graph("TEST", records, "AA","CC")
        self.assertIn("NCT", [n["sequence"] for n in graph["nodes"].values()])
        self.assertEqual(graph["paths"][1]["walk"], ["L","R"])

    def test_2_to_6_and_random_lossless(self):
        rng = random.Random(17)
        for motif in ("AC","CTG","ACGT","AACGT","AACCGT"):
            for i in range(8):
                seq = motif*rng.randrange(2,20)+"TTN"+motif*rng.randrange(2,20)+"A"
                graph = g.build_graph("TEST",[record("A",seq,f"({motif})*")],"AA","CC")
                self.assertEqual("".join(graph["nodes"][n]["sequence"] for n in graph["paths"][0]["walk"]),
                                 "AA"+seq+"CC")

    def test_gfa_standard_path_and_svg_counts(self):
        graph = g.build_graph("TEST",[record("A","CTG"*8)],"AA","CC")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"TEST.gfa"
            g.write_gfa(path,graph)
            fields = next(l for l in path.read_text().splitlines() if l.startswith("P\t")).split("\t")
            self.assertEqual(fields[2], "L+,"+",".join(["N1+"]*8)+",R+")
            self.assertEqual(fields[4], "rp:Z:L+,N1+*8,R+")
            self.assertEqual(json.loads(path.with_suffix(".json").read_text())["paths"][0]["walk"],
                             graph["paths"][0]["walk"])
            g.write_paths_svg(path.with_suffix(".svg"), graph)
            root = ET.parse(path.with_suffix(".svg")).getroot()
            self.assertEqual(len(root.findall(".//*[@class='graph-loop']")),1)
            self.assertEqual(len(root.findall(".//*[@class='graph-edge']")),2)
            self.assertEqual(len(root.findall(".//*[@class='graph-node']")),3)
            self.assertIn("×8", "".join(root.itertext()))
            self.assertNotIn("uncertain", "".join(root.itertext()).lower())

    def test_reject_corrupt_annotation_and_name_collision(self):
        bad = record("A","CTG"*8)
        bad["annotation"]["blocks"][0]["count"] = 9
        with self.assertRaises(ValueError):
            g.build_graph("TEST",[bad],"AA","CC")
        with self.assertRaises(ValueError):
            g.build_graph("TEST",[record("N1","CTG"*8)],"AA","CC")
        with self.assertRaises(ValueError):
            g.build_graph("TEST",[],"AA","CC")


if __name__ == "__main__":
    unittest.main()
