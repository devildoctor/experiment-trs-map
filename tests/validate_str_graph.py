#!/usr/bin/env python3
"""Independent GFA/SVG acceptance checks, without importing the graph writer."""
from collections import Counter
from hashlib import sha256
from itertools import groupby
from pathlib import Path
import argparse
import csv
import json
import xml.etree.ElementTree as ET


def validate(folder, before, locus):
    segments, edges, paths, rle = {}, set(), {}, {}
    for line in (folder/f"{locus}.gfa").read_text(encoding="ascii").splitlines():
        f = line.split("\t")
        if f[0] == "H":
            assert f == ["H", "VN:Z:1.0"]
        elif f[0] == "S":
            assert len(f) == 3 and f[1] not in segments and set(f[2]) <= set("ACGTN")
            segments[f[1]] = f[2]
        elif f[0] == "L":
            assert len(f) == 6 and f[2] == f[4] == "+" and f[5] == "0M"
            assert (f[1], f[3]) not in edges
            edges.add((f[1], f[3]))
        elif f[0] == "P":
            assert len(f) == 5 and f[1] not in paths and f[3] == "*" and f[4].startswith("rp:Z:")
            assert all(n.endswith("+") for n in f[2].split(","))
            paths[f[1]] = [n[:-1] for n in f[2].split(",")]
            rle[f[1]] = []
            for item in f[4][5:].split(","):
                node, count = item.split("*") if "*" in item else (item, "1")
                assert node.endswith("+") and int(count)>0
                rle[f[1]].extend([node[:-1]]*int(count))
        else:
            raise AssertionError(f"Unexpected GFA line: {f[0]}")
    assert not set(segments) & set(paths)
    data = json.loads((folder/f"{locus}.json").read_text())
    records = {r["id"]:r for r in json.loads((folder/"str_paths.json").read_text())["records"]}
    candidates = {r["allele"]:r for r in csv.DictReader((folder/"candidates.tsv").open(),delimiter="\t")}
    expected = {r["id"] for r in records.values() if r["kind"]=="reference" or r.get("supported")}
    assert set(paths) == expected == {p["id"] for p in data["paths"]}
    assert segments == {n:v["sequence"] for n,v in data["nodes"].items()}
    assert edges == {(e["source"],e["target"]) for e in data["edges"]}
    used = set()
    for p in data["paths"]:
        name = p["id"]
        walk = paths[name]
        assert walk == p["walk"] == rle[name]
        assert walk[0]=="L" and walk[-1]=="R"
        adjacent = set(zip(walk,walk[1:]))
        assert adjacent <= edges
        used.update(adjacent)
        internal = "".join(segments[n] for n in walk[1:-1])
        assert sha256(internal.encode()).hexdigest() == records[name]["annotation"]["sequence_sha256"]
        assert len(internal) == records[name]["annotation"]["length"]
        assert sha256("".join(segments[n] for n in walk).encode()).hexdigest() == p["full_sequence_sha256"]
        if name in candidates:
            assert internal == candidates[name]["sequence"]
            assert p["molecule_ids"] == sorted(set(json.loads(candidates[name]["molecule_ids"])))
            assert p["uncertain"] == records[name]["uncertain"]
        else:
            assert name == "REFERENCE" and not p["molecule_ids"]
        runs = [(n,sum(1 for _ in group)) for n,group in groupby(walk)]
        assert runs == [(r["node"],r["count"]) for r in p["runs"]]
        root = ET.parse(folder/"paths-svg"/f"{locus}.{name}.svg").getroot()
        boxes = root.findall(".//*[@class='graph-node']")
        assert [(b.attrib["data-node"],int(b.attrib["data-count"])) for b in boxes] == runs
        for box,(node,count) in zip(boxes,runs):
            texts = [e.text or "" for e in box.findall("{http://www.w3.org/2000/svg}text")]
            assert texts[0] == node
            assert "".join(texts[1:-1] if count>1 else texts[1:]) == segments[node]
            if count>1:
                assert texts[-1] == f"×{count}"
        arrows = root.findall(".//*[@class='graph-edge']")+root.findall(".//*[@class='graph-loop']")
        assert {(e.attrib["data-source"],e.attrib["data-target"]) for e in arrows} == adjacent
        assert len(arrows)==len(adjacent)
    assert used == edges
    for name,node in data["nodes"].items():
        support = set().union(*(set(p["molecule_ids"]) for p in data["paths"] if name in p["walk"]))
        assert node["molecule_ids"] == sorted(support)
        assert node["molecule_support"] == len(support)
    for edge in data["edges"]:
        support = set().union(*(set(p["molecule_ids"]) for p in data["paths"]
                             if (edge["source"],edge["target"]) in set(zip(p["walk"],p["walk"][1:]))))
        assert edge["molecule_ids"] == sorted(support)
        assert edge["molecule_support"] == len(support)
    topology = ET.parse(folder/f"{locus}.svg").getroot()
    ns = {"s":"http://www.w3.org/2000/svg"}
    titles = lambda cls: {e.find("s:title",ns).text for e in topology.findall(f".//*[@class='{cls}']")}
    assert titles("node") == set(segments)
    assert titles("edge") == {a+"->"+b for a,b in edges}
    all_paths = ET.parse(folder/f"{locus}.paths.svg").getroot()
    assert {e.attrib["data-path"] for e in all_paths.findall(".//*[@class='graph-path']")} == set(paths)
    unchanged = ["candidates.tsv","assignments.tsv","genotypes.tsv","str_paths.json","str_paths.tsv",
                 "tokens.tsv","motifs.tsv","status.json","molecules.tsv"]
    for filename in unchanged:
        assert (folder/filename).read_bytes() == (before/filename).read_bytes(), filename
    for ext in ("gfa","json"):
        assert (folder/f"{locus}.detail.{ext}").read_bytes() == (before/f"{locus}.{ext}").read_bytes(), ext
    assert (folder/f"{locus}.detail.svg").read_bytes() == (before/f"{locus}.detail.svg").read_bytes()
    def counts(path):
        result = Counter(line.split("\t")[0] for line in path.read_text().splitlines())
        return {k:result[k] for k in ("S","L","P")}
    summary = dict(status="PASS", before=counts(before/f"{locus}.gfa"),
                   after=counts(folder/f"{locus}.gfa"), self_loops=sum(a==b for a,b in edges),
                   exact_paths_checked=len(paths), candidate_paths=len(paths)-1,
                   per_path_svgs_checked=len(paths), unchanged_files=unchanged,
                   baseline_detail_pair_byte_identical=True,
                   full_sequence_reconstruction=True, rle_expansion=True,
                   svg_node_edge_and_count_match=True, molecule_support_deduplicated=True)
    (folder/"graph_validation.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--result",type=Path,required=True)
    parser.add_argument("--before",type=Path,required=True)
    parser.add_argument("--locus",default="ATXN2")
    args = parser.parse_args()
    validate(args.result,args.before,args.locus)
