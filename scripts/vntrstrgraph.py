#!/usr/bin/env python3
"""Exact STR unit GFA and static views from the same graph model."""
from __future__ import annotations

from collections import defaultdict, deque
from hashlib import sha256
from html import escape
from itertools import groupby
from pathlib import Path
import json
import re
import subprocess

import vntrstr

VERSION = "str_unit_graph_v1"
COLORS = {"repeat": "#dbeafe", "literal": "#fff0cc", "reference_flank": "#edf0f4"}


def compact_walk(walk):
    """Run lengths are visits to a segment, not traversals of its self-link."""
    return [dict(node=node, count=sum(1 for _ in visits))
            for node, visits in groupby(walk)]


def _token(block):
    return (block["kind"], block["motif"] if block["kind"] == "repeat"
            else block["exact_bases"])


def build_graph(locus, records, left, right, molecule_ids=None):
    """Share prefixes, then identical suffix states of count-blind block paths.

    Context remains part of node identity. Equal bases alone never justify
    merging two positions. Counts are stored on each explicit path.
    """
    if not records:
        raise ValueError("An STR graph requires at least one explicit path")
    if not left or not right:
        raise ValueError("Both reference flanks must contain exact bases")
    molecule_ids = molecule_ids or {}
    records = sorted(records, key=lambda r: (r["kind"] != "reference", r["id"]))
    identifiers = [r["id"] for r in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Duplicate path identifier")
    for name in identifiers:
        if not re.fullmatch(r"[A-Za-z0-9_.:-]+", name):
            raise ValueError(f"Unsafe path identifier: {name!r}")

    trie = [dict(token=None, children={}, terminal=False)]
    routes = {}
    for record in records:
        annotation = record["annotation"]
        seq = vntrstr.reconstruct(annotation["blocks"])
        if sha256(seq.encode()).hexdigest() != annotation["sequence_sha256"]:
            raise ValueError("Annotation does not reconstruct its recorded sequence")
        current, route = 0, []
        for block in annotation["blocks"]:
            token = _token(block)
            if token[0] == "repeat" and not 2 <= len(token[1]) <= 6:
                raise ValueError("STR graph repeat units must have 2-6 bp")
            if token not in trie[current]["children"]:
                trie[current]["children"][token] = len(trie)
                trie.append(dict(token=token, children={}, terminal=False))
            current = trie[current]["children"][token]
            route.append(current)
        trie[current]["terminal"] = True
        routes[record["id"]] = route

    # Children are always created after their parent; reverse order is postorder.
    canonical, interned, states = {}, {}, {}
    for index in range(len(trie)-1, 0, -1):
        state = trie[index]
        children = tuple(sorted((token, canonical[child])
                                for token, child in state["children"].items()))
        key = (state["token"], state["terminal"], children)
        if key not in interned:
            cid = len(interned)
            interned[key] = cid
            states[cid] = dict(token=state["token"], terminal=state["terminal"],
                               children=[child for _, child in children])
        canonical[index] = interned[key]

    roots = [canonical[i] for _, i in sorted(trie[0]["children"].items())]
    queue, names = deque(roots), {}
    while queue:
        cid = queue.popleft()
        if cid in names:
            continue
        names[cid] = f"N{len(names)+1}"
        queue.extend(states[cid]["children"])
    if set(identifiers) & ({"L", "R"} | set(names.values())):
        raise ValueError("GFA segment and path identifiers must be distinct")
    nodes = {"L": dict(kind="reference_flank", sequence=left)}
    nodes.update({name: dict(kind=states[cid]["token"][0],
                             sequence=states[cid]["token"][1])
                  for cid, name in names.items()})
    nodes["R"] = dict(kind="reference_flank", sequence=right)
    edges, paths = set(), []
    node_support, edge_support = defaultdict(set), defaultdict(set)
    for record in records:
        walk = ["L"]
        for index, block in zip(routes[record["id"]], record["annotation"]["blocks"]):
            walk.extend([names[canonical[index]]] *
                        (block["count"] if block["kind"] == "repeat" else 1))
        walk.append("R")
        sequence = "".join(nodes[node]["sequence"] for node in walk)
        expected = left + vntrstr.reconstruct(record["annotation"]["blocks"]) + right
        if sequence != expected:
            raise ValueError(f"Lossy graph path: {record['id']}")
        used_edges = set(zip(walk, walk[1:]))
        edges.update(used_edges)
        # A molecule supports a node/link once, however many times it is visited.
        support = set(molecule_ids.get(record["id"], ()))
        for node in set(walk):
            node_support[node].update(support)
        for edge in used_edges:
            edge_support[edge].update(support)
        paths.append(dict(record, walk=walk, runs=compact_walk(walk),
                          molecule_ids=sorted(support),
                          full_sequence_sha256=sha256(sequence.encode()).hexdigest()))
    return dict(version=VERSION, locus=locus, orientation="reference",
                flank_scope="reference_context_25bp; not observed molecule ends",
                count_scope="consecutive segment visits in one explicit path",
                limitation="Candidate paths are not confirmed genotypes. Topology alone "
                           "also admits walks/counts absent from the explicit P records.",
                nodes={name: dict(node, molecule_support=len(node_support[name]),
                                  molecule_ids=sorted(node_support[name]))
                       for name, node in nodes.items()},
                edges=[dict(source=a, target=b, molecule_support=len(edge_support[a,b]),
                            molecule_ids=sorted(edge_support[a,b])) for a,b in sorted(edges)],
                paths=paths)


def write_gfa(path, graph):
    """GFA1 P is authoritative; rp is an optional, ASCII, user-defined RLE tag."""
    lines = ["H\tVN:Z:1.0"]
    lines += [f"S\t{name}\t{node['sequence']}" for name,node in graph["nodes"].items()]
    lines += [f"L\t{e['source']}\t+\t{e['target']}\t+\t0M" for e in graph["edges"]]
    for p in graph["paths"]:
        expanded = ",".join(node+"+" for node in p["walk"])
        compact = ",".join(r["node"]+"+"+(f"*{r['count']}" if r["count"]>1 else "")
                           for r in p["runs"])
        lines.append(f"P\t{p['id']}\t{expanded}\t*\trp:Z:{compact}")
    Path(path).write_text("\n".join(lines)+"\n", encoding="ascii")
    Path(path).with_suffix(".json").write_text(json.dumps(graph, indent=2, ensure_ascii=False),
                                             encoding="utf-8")


def _lines(sequence, width=22):
    return [sequence[i:i+width] for i in range(0, len(sequence), width)] or [""]


def write_topology(path, graph):
    """Graphviz lays out the exact S/L topology; no path-specific counts here."""
    dot = ['digraph STR {', 'graph [rankdir=LR, bgcolor="white", pad=0.3, nodesep=0.45];',
           'node [shape=box, style="rounded,filled", fontname="monospace", fontsize=12];',
           'edge [arrowsize=0.7, color="#59687c"];']
    for name,node in graph["nodes"].items():
        label = "\n".join([name]+_lines(node["sequence"]))
        dot.append(f'{json.dumps(name)} [label={json.dumps(label)}, '
                   f'fillcolor="{COLORS[node["kind"]]}"];')
    for edge in graph["edges"]:
        dot.append(f'{json.dumps(edge["source"])} -> {json.dumps(edge["target"])};')
    dot.append("}")
    dotpath = Path(path).with_suffix(".dot")
    dotpath.write_text("\n".join(dot)+"\n", encoding="utf-8")
    subprocess.run(["dot", "-Tsvg", str(dotpath), "-o", str(path)], check=True)


def write_paths_svg(path, graph, paths=None):
    """One row per explicit path; boxes/arrows/loops keep the GFA identifiers."""
    paths = graph["paths"] if paths is None else paths
    body, top, total_width = [], 52, 900
    for p in paths:
        runs = compact_walk(p["walk"])
        boxes = []
        for run in runs:
            node = graph["nodes"][run["node"]]
            lines = _lines(node["sequence"])
            boxes.append((max(90, max(map(len, lines))*8+24), 44+len(lines)*18, lines))
        row_height = max(b[1] for b in boxes)+105
        y, x = top+47, 24
        body.append(f'<g class="graph-path" data-path="{escape(p["id"], quote=True)}">')
        body.append(f'<text x="24" y="{top+8}" font-size="17" font-weight="bold">{escape(p["id"])}</text>')
        previous = None
        for run,(width,height,lines) in zip(runs, boxes):
            name, count = run["node"], run["count"]
            node = graph["nodes"][name]
            mid = y+height/2
            if previous is not None:
                px, py, pn = previous
                body.append(f'<path class="graph-edge" data-source="{pn}" data-target="{name}" '
                            f'd="M {px} {py} L {x-5} {mid}" fill="none" stroke="#59687c" '
                            f'stroke-width="1.6" marker-end="url(#arrow)"/>')
            body.append(f'<g class="graph-node" data-node="{name}" data-count="{count}">'
                        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="8" '
                        f'fill="{COLORS[node["kind"]]}" stroke="#59687c"/>')
            body.append(f'<text x="{x+width/2}" y="{y+20}" text-anchor="middle" font-size="12" '
                        f'fill="#526175">{name}</text>')
            for j,line in enumerate(lines):
                body.append(f'<text x="{x+width/2}" y="{y+42+j*18}" text-anchor="middle" '
                            f'font-family="monospace" font-size="14">{escape(line)}</text>')
            if count > 1:
                body.append(f'<path class="graph-loop" data-source="{name}" data-target="{name}" '
                            f'd="M {x+width*.7} {y} C {x+width*1.2} {y-38}, '
                            f'{x-width*.2} {y-38}, {x+width*.3} {y-2}" fill="none" '
                            f'stroke="#2563a6" stroke-width="1.6" marker-end="url(#arrow)"/>'
                            f'<text x="{x+width/2}" y="{y-31}" text-anchor="middle" '
                            f'font-size="16" fill="#14518d">×{count}</text>')
            body.append("</g>")
            previous = x+width, mid, name
            x += width+52
        body.append("</g>")
        total_width = max(total_width, x-28)
        top += row_height
    header = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_width}" height="{top+20}" '
              f'viewBox="0 0 {total_width} {top+20}" font-family="sans-serif">'
              '<rect width="100%" height="100%" fill="white"/>'
              '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" '
              'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
              '<path d="M 0 0 L 10 5 L 0 10 z" fill="#59687c"/></marker></defs>'
              f'<text x="24" y="27" font-size="21">{escape(graph["locus"])} | paths</text>')
    Path(path).write_text(header+"".join(body)+"</svg>\n", encoding="utf-8")


def write_outputs(destination, locus, records, left, right, molecule_ids=None):
    graph = build_graph(locus, records, left, right, molecule_ids)
    write_gfa(destination/f"{locus}.gfa", graph)
    write_topology(destination/f"{locus}.svg", graph)
    write_paths_svg(destination/f"{locus}.paths.svg", graph)
    folder = destination/"paths-svg"
    folder.mkdir(exist_ok=True)
    for p in graph["paths"]:
        write_paths_svg(folder/f"{locus}.{p['id']}.svg", graph, [p])
    return graph
