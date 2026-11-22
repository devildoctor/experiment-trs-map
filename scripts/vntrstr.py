"""Opt-in, lossless STR expressions. This module does not call genotypes.

The catalog fixes motif phase and reference direction. Only literal bases and
(2..6 bp motif)* blocks are accepted; no approximate match contributes a copy.
"""
from __future__ import annotations
import csv
from functools import lru_cache
from hashlib import sha256
from html import escape
import json
import re

from vntrmotifs import primitive_unit

VERSION = "1.0-exact-str"


# 解析限定的 STR 语法，并把复合重复词约简为基本周期；不支持任意正则表达式。
@lru_cache(maxsize=128)
def parse_structure(structure):
    tokens, pos = [], 0
    while pos < len(structure):
        match = re.match(r"\(([ACGT]{2,6})\)\*|([ACGT]+)", structure[pos:])
        if not match:
            raise ValueError("STR_STRUCTURE supports only (2..6 bp ACGT motif)* and ACGT literals")
        motif, literal = match.groups()
        if motif:
            primitive = primitive_unit(motif)
            if len(primitive) == 1:
                raise ValueError("Mononucleotide STRs are outside this pilot")
            tokens.append(dict(kind="repeat", motif=primitive, catalog_motif=motif,
                               repeat_factor=len(motif)//len(primitive)))
        else:
            tokens.append(dict(kind="literal", exact_bases=literal))
        pos += len(match[0])
    if not any(t["kind"] == "repeat" for t in tokens):
        raise ValueError("STR_STRUCTURE requires a repeat block")
    return tuple(tokens)


# 先检查原始目录单元长度，再约简周期，防止长 VNTR 借约简进入 STR 模式。
def validate_config(structure, catalog_motifs):
    if structure:
        parse_structure(structure)
        # Check original catalog lengths BEFORE primitive reduction.
        if not catalog_motifs or any(not 2 <= len(m) <= 6 for m in catalog_motifs):
            raise ValueError("STR_STRUCTURE requires explicit 2..6 bp catalog motifs; VNTR is excluded")
    return structure


# 重复块按次数展开，其他块使用原始字面碱基，得到完整输入序列。
def reconstruct(blocks):
    return "".join(b["motif"]*b["count"] if b["kind"] == "repeat"
                   else b["exact_bases"] for b in blocks)


# 生成便于阅读的表达式；表达式次数只描述这条序列，不直接表示基因型。
def format_blocks(blocks):
    return "–".join(f'({b["motif"]}){b["count"]}' if b["kind"] == "repeat"
                    else b["exact_bases"] for b in blocks) or "(empty)"


# 反向互补时同时反转区块顺序、坐标和字面序列，不能只替换 motif。
def reverse_blocks(blocks):
    """Reverse-complement the WHOLE path, including literals and their order."""
    def rc(s):
        return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]
    n = max((b["end"] for b in blocks), default=0)
    result = []
    for b in reversed(blocks):
        item = dict(b, start=n-b["end"], end=n-b["start"], exact_bases=rc(b["exact_bases"]))
        if b["kind"] == "repeat":
            item["motif"] = rc(b["motif"])
        result.append(item)
    return result


# 把精确重复匹配之间的所有间隙和末端补回为字面序列。
def _blocks(sequence, path):
    blocks, pos = [], 0
    for start, end, motif in path:
        if start > pos:
            blocks.append(dict(kind="literal", start=pos, end=start, exact_bases=sequence[pos:start]))
        blocks.append(dict(kind="repeat", start=start, end=end, motif=motif,
                           count=(end-start)//len(motif), exact_bases=sequence[start:end]))
        pos = end
    if pos < len(sequence):
        blocks.append(dict(kind="literal", start=pos, end=len(sequence), exact_bases=sequence[pos:]))
    return blocks


# 固定模板相位，以字面碱基最少、重复块最少的顺序求解，并保留一个备选。
def annotate(sequence, structure):
    """Two best distinct exact decompositions; fixed configured phase.

    Minimize (literal bases, repeat blocks) lexicographically. Retain one
    alternative, including its cost; never interpret this cost as probability.
    A run needs >=2 copies. Single units/tails remain literal sequence.
    """
    model = parse_structure(structure)
    # 切分使用模型中的 motif 集合；完整模型的先后顺序和字面片段由 model_match 另行检查。
    motifs = sorted({t["motif"] for t in model if t["kind"] == "repeat"})
    n = len(sequence)
    dp = [[] for _ in range(n+1)]
    dp[n] = [((0, 0), ())]
    for i in range(n-1, -1, -1):
        choices = {(cost[0]+1, cost[1], path) for cost, path in dp[i+1]}
        for motif in motifs:
            end, count = i, 0
            while sequence.startswith(motif, end):
                count += 1
                end += len(motif)
                if count < 2:
                    continue
                for cost, path in dp[end]:
                    choices.add((cost[0], cost[1]+1, ((i, end, motif),)+path))
        ranked = sorted(choices)
        # The same repeat path can be reached through multiple literal edges.
        dp[i] = [((a, b), path) for a, b, path in ranked[:2]]
    cost, path = dp[0][0]
    blocks = _blocks(sequence, path)
    alternative = None
    if len(dp[0]) > 1:
        alt_cost, alt_path = dp[0][1]
        alt_blocks = _blocks(sequence, alt_path)
        alternative = dict(blocks=alt_blocks, expression=format_blocks(alt_blocks),
                           cost=list(alt_cost), equal_cost=alt_cost == cost)
        assert reconstruct(alt_blocks) == sequence
    assert reconstruct(blocks) == sequence
    pattern = "".join(f'(?:{t["catalog_motif"]})*' if t["kind"] == "repeat"
                      else t["exact_bases"] for t in model)
    reverse = reverse_blocks(blocks)
    return dict(version=VERSION, structure=structure, model=list(model),
                expression=format_blocks(blocks), blocks=blocks,
                reverse_complement_expression=format_blocks(reverse),
                reverse_complement_blocks=reverse,
                sequence_sha256=sha256(sequence.encode()).hexdigest(), length=n,
                model_match=bool(re.fullmatch(pattern, sequence)), literal_bases=cost[0],
                cost=list(cost), cost_definition="lexicographic_literal_bases_then_repeat_blocks",
                annotation_ambiguity=bool(alternative and alternative["equal_cost"]),
                alternative=alternative, alternatives_exhaustive=False,
                phase_policy="fixed_configured_motif_phase", orientation="reference",
                reconstruction_ok=True, copy_number_scope="exact_runs_in_this_sequence_only")


# 为注释绑定路径 ID 和来源元数据，保持候选与观测片段的身份区别。
def path_record(identifier, sequence, structure, **metadata):
    return dict(id=identifier, **metadata, annotation=annotate(sequence, structure))


# 给每个已有观测片段附加 STR 注释；不更改提取边界或原始序列。
def annotate_rows(rows, structure):
    records = []
    for row in rows:
        for fragment in row["fragments"]:
            record = path_record(fragment["id"], fragment["sequence"], structure,
                                 kind="observed_fragment", sample=row["sample"], platform=row["platform"],
                                 row_id=row["id"], read_id=row["read_id"],
                                 coverage=fragment["coverage"], anchor=fragment["anchor"],
                                 qstart=fragment["qstart"], qend=fragment["qend"])
            fragment["str_annotation"] = record["annotation"]
            records.append(record)
    return records


# 结构框不按 bp 等比例绘制；长字面序列通过换行显示而非省略。
def svg_blocks(annotation, x=20, y=0, width=1460):
    """Wrapped, non-proportional exact blocks; every literal base is displayed."""
    parts, cursor, top = [], x, y
    for b in annotation["blocks"]:
        label = f'{b["motif"]} × {b["count"]}' if b["kind"] == "repeat" else b["exact_bases"]
        # Long literals wrap within the chip, never turn into unknown placeholders.
        chars = max(8, int((width-24)/8))
        lines = [label[i:i+chars] for i in range(0, len(label), chars)] or [""]
        box_width = min(width, max(80, min(len(label), chars)*8+24))
        box_height = 34+18*(len(lines)-1)
        if cursor+box_width > x+width and cursor > x:
            cursor, top = x, top+52
        fill = "#286c96" if b["kind"] == "repeat" else "#f3e2b6"
        foreground = "white" if b["kind"] == "repeat" else "#573d12"
        title = f'[{b["start"]},{b["end"]}) {b["exact_bases"]}'
        parts.append(f'<g class="str-block" data-kind="{b["kind"]}" data-start="{b["start"]}" data-end="{b["end"]}"><title>{escape(title)}</title><rect x="{cursor}" y="{top}" width="{box_width}" height="{box_height}" rx="5" fill="{fill}"/>')
        for j, line in enumerate(lines):
            parts.append(f'<text x="{cursor+12}" y="{top+22+18*j}" font-family="monospace" font-size="13" fill="{foreground}">{escape(line)}</text>')
        parts.append("</g>")
        cursor += box_width+24
        if len(lines) > 1:
            cursor, top = x, top+box_height+18
    return "".join(parts), top-y+52


# 逐路径显示注释、覆盖和候选不确定性，保留支持信息的原有口径。
def write_paths_svg(path, title, records, reason=""):
    body, y = [], 110
    if reason:
        for start in range(0,len(reason),165):
            body.append(f'<text x="20" y="{y}" font-size="13" fill="#994444">{escape(reason[start:start+165])}</text>')
            y += 20
        y += 15
    for record in records:
        a = record["annotation"]
        label = f'{record["id"]} | {record.get("kind","")} | {record.get("sample","")} | {a["length"]} bp | {record.get("coverage","")}'
        body.append(f'<text x="20" y="{y}" font-size="13">{escape(label)}</text>')
        diagram, height = svg_blocks(a, y=y+14)
        body.append(diagram)
        y += height+20
        if record["kind"] == "candidate":
            support = (f'call_status={record["call_status"]}; molecules={record["read_support"]}; '
                       f'platforms={json.dumps(record["platform_support"],sort_keys=True)}; '
                       f'ambiguous_molecules={record["ambiguous_molecules"]}')
            body.append(f'<text x="20" y="{y}" font-size="12">{escape(support)}</text>')
            y += 20
        note = f'model_match={a["model_match"]}; literal_bases={a["literal_bases"]}; STR_split_ambiguous={a["annotation_ambiguity"]}'
        if record.get("uncertain"):
            note += "; UNCERTAIN: "+",".join(record.get("flags", []))
        for start in range(0,len(note),165):
            body.append(f'<text x="20" y="{y}" font-size="12" fill="#72501c">{escape(note[start:start+165])}</text>')
            y += 18
        y += 22
    header = f'<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="{y}" viewBox="0 0 1500 {y}" font-family="sans-serif"><rect width="100%" height="100%" fill="white"/><text x="20" y="32" font-size="24">{escape(title)}</text><text x="20" y="58" font-size="13">Exact STR runs + literal bases; fixed reference direction. Counts describe this sequence, not a genotype.</text><text x="20" y="80" font-size="12">Blue = exact repeated motif; gold = literal sequence. Chip widths are NOT proportional to base pairs.</text>'
    path.write_text(header+"".join(body)+"</svg>", encoding="utf-8")


# JSON 保存完整切分和备选，TSV 保存摘要，SVG 每页最多 40 条且不抽样。
def write_outputs(destination, locus, structure, records, figures=True):
    payload = dict(version=VERSION, locus=locus, structure=structure, records=records,
                   limitation="Exact configured-motif annotation only; no new genotype inference or boundary adjustment")
    (destination/"str_paths.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    fields = ["id", "kind", "sample", "platform", "coverage", "anchor", "call_status", "read_support", "length",
              "expression", "sequence_sha256", "model_match", "literal_bases", "annotation_ambiguity"]
    with (destination/"str_paths.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(dict(record, **record["annotation"]))
    if not figures:
        return
    folder = destination/"str-svg"
    folder.mkdir(exist_ok=True)
    links = []
    for start in range(0, len(records), 40):
        name = f"{locus}.str.{start//40+1:04d}.svg"
        write_paths_svg(folder/name, f"{locus} STR expressions | page {start//40+1}", records[start:start+40])
        links.append(f'<li><a href="{name}">{start+1}–{min(start+40,len(records))}</a></li>')
    (folder/"index.html").write_text('<!doctype html><meta charset="utf-8"><h2>All STR paths, no sampling</h2><ol>'+"".join(links)+"</ol>", encoding="utf-8")
