#!/usr/bin/env python3
"""从双侧翼跨越 read 构建可回溯的 VNTR 位点图。仅依赖 Python 标准库。"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from difflib import SequenceMatcher
from collections import Counter, defaultdict
from html import escape
from pathlib import Path
import re
import shutil
import subprocess


DNA = set("ACGT")
COLORS = ("#4477AA", "#EE7733", "#228833", "#AA3377", "#66CCEE", "#CCBB44")


@dataclass(frozen=True)
class Locus:
    name: str
    chrom: str
    start: int
    end: int
    motifs: tuple[str, ...]


@dataclass(frozen=True)
class Read:
    name: str
    sample: str
    platform: str
    sequence: str


@dataclass
class Allele:
    name: str
    sample: str
    sequence: str
    support: dict[str, int]
    read_ids: tuple[str, ...]
    identity: float
    tokens: list[tuple[str, str]]
    uncertain: bool = False
    dbg_k: int = 0
    graph_support: float = 0.0
    unitigs: int = 0


class IndexedFastaReference:
    """独立 FASTA/FAI 读取器，所有坐标均为 0-based 半开区间。"""

    def __init__(self, path: Path):
        self.path = path
        fai = Path(str(path) + ".fai")
        if not path.is_file() or not fai.is_file():
            raise ValueError("参考 FASTA 或同名 FAI 不存在")
        self.entries = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            name, length, offset, line_bases, line_width, *_ = line.split("\t")
            self.entries[name] = tuple(map(int, (length, offset, line_bases, line_width)))
        self.handle = path.open("rb")

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.handle.close()

    def contig_length(self, chrom: str) -> int:
        return self.entries[chrom][0]

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if chrom not in self.entries:
            raise ValueError(f"参考基因组无染色体 {chrom}")
        length, offset, line_bases, line_width = self.entries[chrom]
        if start < 0 or end > length or end < start:
            raise ValueError(f"参考坐标越界: {chrom}:{start}-{end}")
        parts = []
        while start < end:
            row, column = divmod(start, line_bases)
            size = min(end - start, line_bases - column)
            self.handle.seek(offset + row * line_width + column)
            part = self.handle.read(size)
            if len(part) != size:
                raise ValueError("FASTA/FAI 不匹配")
            parts.append(part)
            start += size
        return b"".join(parts).decode("ascii").upper()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def global_unique_probes(reference: Path, loci: list[Locus], fasta: IndexedFastaReference,
                         flank_length: int) -> dict[str, bool]:
    """流式扫描整套 FASTA，确认两侧至少有一个唯一的 31-mer。"""
    probes: dict[str, list[bytes]] = {}
    for locus in loci:
        for side, flank in (("L", fasta.fetch(locus.chrom, locus.start - flank_length, locus.start)),
                            ("R", fasta.fetch(locus.chrom, locus.end, locus.end + flank_length))):
            positions = sorted({0, max(0, (len(flank) - 31) // 2), len(flank) - 31})
            probes[f"{locus.name}:{side}"] = [flank[i:i + 31].encode("ascii") for i in positions]
    patterns = {p for group in probes.values() for p in group if len(p) == 31 and set(p) <= set(b"ACGT")}
    patterns |= {reverse_complement(p.decode("ascii")).encode("ascii") for p in patterns}
    counts = {p: 0 for p in patterns}

    def scan(chunk: bytes) -> None:
        for pattern in patterns:
            begin = 0
            while counts[pattern] <= 2:
                hit = chunk.find(pattern, begin)
                if hit < 0:
                    break
                counts[pattern] += 1
                begin = hit + 1

    with reference.open("rb") as handle:
        chunk = bytearray()
        for raw in handle:
            if raw.startswith(b">"):
                scan(bytes(chunk))
                chunk.clear()
                continue
            chunk.extend(raw.strip().upper())
            if len(chunk) >= 4_000_000:
                scan(bytes(chunk))
                chunk = chunk[-30:]
        scan(bytes(chunk))
    answer = {}
    for key, group in probes.items():
        answer[key] = any(counts.get(p, 0) + (0 if p == reverse_complement(p.decode("ascii")).encode("ascii")
                                             else counts.get(reverse_complement(p.decode("ascii")).encode("ascii"), 0)) == 1
                          for p in group if len(p) == 31)
    return answer


def read_catalog(path: Path) -> list[Locus]:
    """同时接受带 ID/MOTIFS 的 TRGT BED 与四列候选 BED。"""
    loci = []
    seen = set()
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"BED 第 {number} 行少于四列")
        chrom, a, b, meta = fields[:4]
        start, end = int(a), int(b)
        if start < 0 or end <= start:
            raise ValueError(f"BED 第 {number} 行坐标无效")
        info = dict(item.split("=", 1) for item in meta.split(";") if "=" in item)
        name = info.get("ID", meta.split("_")[0])
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", name) or name in seen:
            raise ValueError(f"BED 第 {number} 行 ID 无效或重复: {name}")
        motifs = tuple(x.upper() for x in info.get("MOTIFS", "").split(",") if x)
        if any(set(x) - DNA for x in motifs):
            raise ValueError(f"BED 第 {number} 行 motif 不是精确 ACGT")
        loci.append(Locus(name, chrom, start, end, motifs))
        seen.add(name)
    if not loci:
        raise ValueError("BED 中没有位点")
    return loci


def discover_bams(folder: Path, locus: str) -> list[tuple[str, str, Path]]:
    """从 bamdata/样本/平台/文件名.位点.bam 取得输入。"""
    found = []
    for bam in sorted(folder.glob("*/*/*.bam")):
        if bam.name.rsplit(".", 2)[-2] != locus:
            continue
        if not (Path(str(bam) + ".bai").exists() or bam.with_suffix(".bai").exists()):
            raise ValueError(f"缺少 BAM 索引: {bam}")
        sample, platform = bam.parent.parent.name, bam.parent.name
        if platform not in {"HiFi", "ONT", "ONTUL"}:
            raise ValueError(f"未知测序平台: {bam}")
        found.append((sample, platform, bam))
    return found


def reference_to_query(pos0: int, cigar: str, target: int) -> int | None:
    """CIGAR 投影仅用于侧翼搜索的初始范围，不决定重复边界。"""
    ref, query = pos0, 0
    for size_text, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        size = int(size_text)
        if op in "M=X":
            if ref <= target <= ref + size:
                return query + target - ref
            ref += size
            query += size
        elif op in "DN":
            if ref <= target <= ref + size:
                return query
            ref += size
        elif op in "IS":
            query += size
    return None


def flank_hit(flank: str, sequence: str, expected: int | None,
              minimum_identity: float) -> tuple[int, int, float] | None:
    """在 read 中独立寻找一个高分且位置唯一的侧翼。"""
    length = len(flank)
    starts = set()
    begin = sequence.find(flank)
    while begin >= 0:
        starts.add(begin)
        begin = sequence.find(flank, begin + 1)
    if not starts:
        # 先找短精确种子，再对整段侧翼评分；允许 ONT 的局部 indel。
        for offset in (0, max(0, length // 2 - 10), max(0, length - 20)):
            seed = flank[offset:offset + 20]
            begin = sequence.find(seed)
            while begin >= 0:
                starts.update(range(max(0, begin - offset - 5), min(len(sequence), begin - offset + 6)))
                begin = sequence.find(seed, begin + 1)
    if expected is not None:
        starts = {x for x in starts if abs(x - expected) <= max(200, length)}
    if not starts:
        return None
    scored = []
    for start in starts:
        for delta in (-8, -4, 0, 4, 8):
            end = start + length + delta
            if end > len(sequence) or end <= start:
                continue
            segment = sequence[start:end]
            score = SequenceMatcher(None, flank, segment, autojunk=False).ratio()
            if score >= minimum_identity:
                scored.append((score, start, end))
    if not scored:
        return None
    scored.sort(reverse=True)
    best = scored[0]
    # 第二个独立位置得分相近时，不能唯一定位边界。
    if any(abs(item[1] - best[1]) > 10 and item[0] >= best[0] - .02 for item in scored[1:]):
        return None
    return best[1], best[2], best[0]


def extract_from_bam(bam: Path, sample: str, platform: str, locus: Locus,
                     reference: IndexedFastaReference, samtools: str,
                     flank_length: int, min_mapq: int) -> tuple[list[Read], Counter]:
    """按位点流式读取 BAM，并按 read 名去重主/补充比对。"""
    if not (Path(str(bam) + ".bai").exists() or bam.with_suffix(".bai").exists()):
        raise ValueError(f"缺少 BAI: {bam}")
    left = reference.fetch(locus.chrom, locus.start - flank_length, locus.start)
    right = reference.fetch(locus.chrom, locus.end, locus.end + flank_length)
    region = f"{locus.chrom}:{max(1, locus.start - flank_length * 2 + 1)}-{locus.end + flank_length * 2}"
    command = [samtools, "view", "-h", str(bam), region]
    process = subprocess.Popen(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert process.stdout is not None
    best: dict[str, tuple[float, Read]] = {}
    audit = Counter()
    header_length = None
    for line in process.stdout:
        if line.startswith("@SQ"):
            tags = dict(item.split(":", 1) for item in line.rstrip().split("\t")[1:] if ":" in item)
            if tags.get("SN") == locus.chrom:
                header_length = int(tags["LN"])
            continue
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            audit["invalid_sam"] += 1
            continue
        name, flag_text, chrom, pos_text, mapq_text, cigar = fields[:6]
        flag, pos0, mapq = int(flag_text), int(pos_text) - 1, int(mapq_text)
        if flag & (0x4 | 0x100) or chrom != locus.chrom:
            audit["unmapped_or_secondary"] += 1
            continue
        if mapq < min_mapq or fields[9] == "*":
            audit["low_mapq_or_no_sequence"] += 1
            continue
        sequence = fields[9].upper()
        if set(sequence) - DNA:
            audit["ambiguous_bases"] += 1
            continue
        q_left = reference_to_query(pos0, cigar, locus.start - flank_length)
        q_right = reference_to_query(pos0, cigar, locus.end)
        identity = .90 if platform == "HiFi" else .72
        left_hit = flank_hit(left, sequence, q_left, identity)
        right_hit = flank_hit(right, sequence, q_right, identity)
        if left_hit is None or right_hit is None:
            audit["missing_or_ambiguous_anchor"] += 1
            continue
        if left_hit[1] > right_hit[0]:
            audit["wrong_anchor_order"] += 1
            continue
        repeat = sequence[left_hit[1]:right_hit[0]]
        if not repeat:
            audit["empty_repeat"] += 1
            continue
        score = min(left_hit[2], right_hit[2]) + .0001 * mapq
        read = Read(name, sample, platform, repeat)
        if name not in best or score > best[name][0]:
            best[name] = (score, read)
        audit["spanning_alignment"] += 1
    stderr = process.stderr.read() if process.stderr else ""
    if process.wait() != 0:
        raise RuntimeError(f"samtools view 失败: {bam}: {stderr.strip()}")
    if header_length != reference.contig_length(locus.chrom):
        raise ValueError(f"BAM 与参考基因组染色体长度不符: {bam}")
    audit["unique_spanning_reads"] = len(best)
    return [item[1] for item in best.values()], audit


def similarity(a: str, b: str) -> float:
    """限制长度差后计算全序列相似度；关闭 autojunk 以保留重复序列。"""
    if not a or not b or abs(len(a) - len(b)) > max(10, .20 * max(len(a), len(b))):
        return 0.0
    return SequenceMatcher(None, a, b, autojunk=False).ratio()


def make_consensus(reads: list[Read]) -> str:
    """以中心 read 为骨架，对每个位置和插入串作简单多数表决。"""
    seed = Counter(r.sequence for r in reads).most_common(1)[0][0]
    columns = [Counter() for _ in seed]
    insertions = [Counter() for _ in range(len(seed) + 1)]
    for read in reads:
        matcher = SequenceMatcher(None, seed, read.sequence, autojunk=False)
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag == "insert":
                insertions[i1][read.sequence[j1:j2]] += 1
            elif tag == "delete":
                for i in range(i1, i2):
                    columns[i]["-"] += 1
            elif tag == "equal":
                for i, j in zip(range(i1, i2), range(j1, j2)):
                    columns[i][read.sequence[j]] += 1
            else:
                width = min(i2 - i1, j2 - j1)
                for offset in range(width):
                    columns[i1 + offset][read.sequence[j1 + offset]] += 1
                for i in range(i1 + width, i2):
                    columns[i]["-"] += 1
                if j2 - j1 > width:
                    insertions[i2][read.sequence[j1 + width:j2]] += 1
    result = []
    for i in range(len(seed) + 1):
        if insertions[i]:
            insertion, count = insertions[i].most_common(1)[0]
            if count > len(reads) / 2:
                result.append(insertion)
        if i < len(seed):
            base = columns[i].most_common(1)[0][0] if columns[i] else seed[i]
            if base != "-":
                result.append(base)
    return "".join(result)


def cluster_reads(reads: list[Read], minimum: int) -> list[Allele]:
    """先分平台聚类，再用完整序列共识和回贴比例过滤错误簇。"""
    groups = defaultdict(list)
    for read in reads:
        groups[(read.sample, read.platform)].append(read)
    candidates = []
    for (sample, platform), source in sorted(groups.items()):
        clusters: list[list[Read]] = []
        threshold = .96 if platform == "HiFi" else .82
        exact_counts = Counter(r.sequence for r in source)
        stable_sequences = {sequence for sequence, count in exact_counts.items() if count >= minimum}
        for sequence in sorted(stable_sequences):
            clusters.append([r for r in source if r.sequence == sequence])
        remaining = [r for r in source if r.sequence not in stable_sequences]
        for read in sorted(remaining, key=lambda r: (len(r.sequence), r.name)):
            scores = [similarity(read.sequence, group[0].sequence) for group in clusters]
            if scores and max(scores) >= threshold:
                clusters[scores.index(max(scores))].append(read)
            else:
                clusters.append([read])
        for group in clusters:
            if len(group) < minimum:
                continue
            consensus = make_consensus(group)
            identities = [similarity(consensus, r.sequence) for r in group]
            if not consensus or sum(x >= threshold for x in identities) < minimum:
                continue
            candidates.append(Allele("", sample, consensus, {platform: len(group)},
                                     tuple(f"{platform}:{r.name}" for r in group),
                                     sum(identities) / len(group), []))
    # 跨平台只合并全长、全碱基一致的共识；单碱基差异先保留为独立候选。
    combined: list[Allele] = []
    for candidate in sorted(candidates, key=lambda a: (-sum(a.support.values()), a.sample)):
        match = next((x for x in combined if x.sample == candidate.sample
                      and x.sequence == candidate.sequence
                      and not (set(x.support) & set(candidate.support))), None)
        if match is None:
            combined.append(candidate)
        else:
            for platform, count in candidate.support.items():
                match.support[platform] = match.support.get(platform, 0) + count
            match.read_ids += candidate.read_ids
            match.identity = min(match.identity, candidate.identity)
    by_sample = Counter()
    for allele in combined:
        by_sample[allele.sample] += 1
        allele.name = f"{allele.sample}_allele{by_sample[allele.sample]}"
    for allele in combined:
        allele.uncertain = by_sample[allele.sample] > 2
    return combined


def dbg_evidence(allele: Allele, reads: list[Read]) -> tuple[int, float, int]:
    """建立 read 支持的有向局部 de Bruijn 图并核查共识路径。"""
    selected = [r for r in reads if r.sample == allele.sample
                and f"{r.platform}:{r.name}" in allele.read_ids]
    if not selected:
        return 0, 0.0, 0
    k = next((size for size in (21, 31, 51) if len(allele.sequence) >= size * 2
              and len({allele.sequence[i:i + size] for i in range(len(allele.sequence) - size + 1)})
              / (len(allele.sequence) - size + 1) >= .35), 51)
    k = min(k, max(3, len(allele.sequence) // 2))
    edges: dict[tuple[str, str], set[str]] = defaultdict(set)
    incoming: dict[str, set[str]] = defaultdict(set)
    outgoing: dict[str, set[str]] = defaultdict(set)
    for read in selected:
        sequence = read.sequence
        for i in range(max(0, len(sequence) - k)):
            a, b = sequence[i:i + k], sequence[i + 1:i + k + 1]
            edges[(a, b)].add(f"{read.platform}:{read.name}")
            outgoing[a].add(b)
            incoming[b].add(a)
    path = [(allele.sequence[i:i + k], allele.sequence[i + 1:i + k + 1])
            for i in range(max(0, len(allele.sequence) - k))]
    supported = sum(bool(edges.get(edge)) for edge in path) / max(1, len(path))
    # 一入一出链在图上可压成 unitig；分叉节点保留。
    nodes = set(incoming) | set(outgoing)
    unitigs = sum(len(incoming[node]) != 1 or len(outgoing[node]) != 1 for node in nodes)
    return k, supported, max(1, unitigs)


def discover_motifs(sequences: list[str], supplied: tuple[str, ...]) -> tuple[str, ...]:
    """根据自相关找短周期；仅强周期作为新 motif，避免臆测复杂序列。"""
    motifs = list(supplied)
    for sequence in sequences:
        best = None
        for period in range(2, min(100, len(sequence) // 2) + 1):
            total = len(sequence) - period
            agreement = sum(sequence[i] == sequence[i + period] for i in range(total)) / total
            if agreement >= .82 and (best is None or (agreement - .0005 * period) > best[0]):
                best = (agreement - .0005 * period, period)
        if best:
            period = best[1]
            for offset in range(min(period, len(sequence) // period)):
                motif = sequence[offset * period:(offset + 1) * period]
                if motif not in motifs:
                    motifs.append(motif)
        if len(motifs) >= 16:
            break
    return tuple(motifs[:16])


def bounded_edit_distance(a: str, b: str, limit: int) -> int:
    """只计算小差异；超过上限立即返回，包含替换、插入和缺失。"""
    if abs(len(a) - len(b)) > limit:
        return limit + 1
    previous = {j: j for j in range(min(len(b), limit) + 1)}
    for i, base in enumerate(a, 1):
        current = {}
        for j in range(max(0, i - limit), min(len(b), i + limit) + 1):
            if j == 0:
                current[j] = i
            else:
                current[j] = min(previous.get(j, limit + 1) + 1,
                                 current.get(j - 1, limit + 1) + 1,
                                 previous.get(j - 1, limit + 1) + (base != b[j - 1]))
        if min(current.values(), default=limit + 1) > limit:
            return limit + 1
        previous = current
    return previous.get(len(b), limit + 1)


def decompose(sequence: str, motifs: tuple[str, ...]) -> tuple[list[tuple[str, str]], float]:
    """DP 分解；状态含前一 motif，允许小 indel、切换惩罚和 UNKNOWN。"""
    n = len(sequence)
    states: list[dict[str | None, tuple[float, tuple | None]]] = [dict() for _ in range(n + 1)]
    states[0][None] = (0.0, None)
    for i in range(n):
        if not states[i]:
            continue
        for previous_family, (cost, _) in list(states[i].items()):
            unknown_score = cost + 1.15
            old = states[i + 1].get(previous_family)
            if old is None or unknown_score < old[0]:
                states[i + 1][previous_family] = (
                    unknown_score, (i, previous_family, "UNKNOWN", sequence[i:i + 1]))
            for motif in motifs:
                length = len(motif)
                allowed = 1 if length <= 6 else min(3, max(2, int(length * .10)))
                for width in range(max(1, length - allowed), length + allowed + 1):
                    if i + width > n:
                        continue
                    exact = sequence[i:i + width]
                    differences = bounded_edit_distance(exact, motif, allowed)
                    if differences > allowed:
                        continue
                    switch = .25 if previous_family not in (None, motif) else 0.0
                    score = cost + .15 + differences + switch
                    old = states[i + width].get(motif)
                    if old is None or score < old[0]:
                        states[i + width][motif] = (score, (i, previous_family, motif, exact))
    family = min(states[n], key=lambda item: states[n][item][0])
    tokens = []
    position = n
    while position:
        step = states[position][family][1]
        assert step is not None
        start, previous_family, token_family, exact = step
        tokens.append((token_family, exact))
        position, family = start, previous_family
    tokens.reverse()
    merged = []
    for family, exact in tokens:
        if family == "UNKNOWN" and merged and merged[-1][0] == "UNKNOWN":
            merged[-1] = ("UNKNOWN", merged[-1][1] + exact)
        else:
            merged.append((family, exact))
    explained = sum(len(exact) for family, exact in merged if family != "UNKNOWN") / max(1, n)
    assert "".join(exact for _, exact in merged) == sequence
    return merged, explained


def build_graph(alleles: list[Allele], left: str, right: str):
    """共享精确序列节点；边只来自已观察到的等位路径。"""
    nodes = {("LEFT", left): "L", ("RIGHT", right): "R"}
    edges = set()
    paths = []
    for allele in alleles:
        ids = ["L"]
        for family, exact in allele.tokens:
            key = (family, exact)
            if key not in nodes:
                nodes[key] = f"N{len(nodes) - 1}"
            ids.append(nodes[key])
        ids.append("R")
        edges.update(zip(ids, ids[1:]))
        paths.append((allele, ids))
        reconstructed = "".join(next(seq for (_, seq), node in nodes.items() if node == item) for item in ids)
        if reconstructed != left + allele.sequence + right:
            raise ValueError(f"路径无法重建: {allele.name}")
    return nodes, edges, paths


def write_gfa(path: Path, nodes, edges, paths) -> None:
    with path.open("w", encoding="utf-8") as out:
        out.write("H\tVN:Z:1.0\n")
        for (family, exact), node in nodes.items():
            out.write(f"S\t{node}\t{exact}\tFT:Z:{family}\n")
        for start, end in sorted(edges):
            out.write(f"L\t{start}\t+\t{end}\t+\t0M\n")
        for allele, ids in paths:
            support = ",".join(f"{key}:{value}" for key, value in sorted(allele.support.items()))
            status = "uncertain" if allele.uncertain else "supported"
            out.write(f"P\t{allele.name}\t{','.join(x + '+' for x in ids)}\t*"
                      f"\tRC:i:{sum(allele.support.values())}\tPS:Z:{support}"
                      f"\tRI:f:{allele.identity:.4f}\tST:Z:{status}\n")


def write_svg(path: Path, locus: str, alleles: list[Allele], reason: str = "") -> None:
    height = max(150, 85 + 95 * len(alleles))
    families = list(dict.fromkeys(f for a in alleles for f, _ in a.tokens))
    rows = []
    width = 1200
    for allele in alleles:
        runs = []
        for token in allele.tokens:
            if runs and runs[-1][0] == token:
                runs[-1][1] += 1
            else:
                runs.append([token, 1])
        boxes = []
        row_width = 80
        for (family, exact), count in runs:
            label = exact if family != "UNKNOWN" else f"UNKNOWN {len(exact)} bp"
            label += f" ×{count}" if count > 1 else ""
            boxwidth = max(76, len(label) * 8 + 16)
            boxes.append((family, label, boxwidth))
            row_width += boxwidth + 24
        width = max(width, row_width + 30)
        rows.append((allele, boxes))
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<defs><marker id="arrow" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto"><path d="M0,0 L8,3 L0,6" fill="#555"/></marker></defs>',
             f'<text x="30" y="38" font-family="sans-serif" font-size="24">{escape(locus)} VNTR</text>']
    if reason:
        lines.append(f'<text x="30" y="95" font-family="sans-serif" font-size="18" fill="#994444">{escape(reason)}</text>')
    for index, (allele, boxes) in enumerate(rows):
        y = 80 + index * 95
        lines.append(f'<text x="30" y="{y}" font-family="sans-serif" font-size="15">{escape(allele.name)}</text>')
        x = 30
        lines.append(f'<text x="{x}" y="{y + 35}" font-size="13">左</text>')
        x += 28
        for family, label, boxwidth in boxes:
            color = "#999999" if family == "UNKNOWN" else COLORS[families.index(family) % len(COLORS)]
            lines.append(f'<line x1="{x - 17}" y1="{y + 29}" x2="{x - 3}" y2="{y + 29}" stroke="#555" marker-end="url(#arrow)"/>')
            lines.append(f'<rect x="{x}" y="{y + 12}" width="{boxwidth}" height="34" rx="5" fill="{color}"/>')
            lines.append(f'<text x="{x + 8}" y="{y + 35}" font-family="monospace" font-size="13" fill="white">{escape(label)}</text>')
            x += boxwidth + 24
        lines.append(f'<line x1="{x - 17}" y1="{y + 29}" x2="{x - 3}" y2="{y + 29}" stroke="#555" marker-end="url(#arrow)"/>')
        lines.append(f'<text x="{x}" y="{y + 35}" font-size="13">右</text>')
        copies = sum(1 for family, _ in allele.tokens if family != "UNKNOWN")
        support = ", ".join(f"{p}:{n}" for p, n in sorted(allele.support.items()))
        note = f"{support}; {copies} 完整单元" if not allele.uncertain else f"{support}; 状态待定，拷贝数不报告"
        lines.append(f'<text x="30" y="{y + 65}" font-family="sans-serif" font-size="13" fill="#444">{escape(note)}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def process_locus(locus: Locus, args, reference: IndexedFastaReference,
                  unique: dict[str, bool]) -> None:
    destination = args.out / locus.name
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / "parameters.tsv").open("w", encoding="utf-8") as out:
        out.write("parameter\tvalue\n")
        for key in ("catalog", "reference", "min_support", "anchor_length", "min_mapq", "samtools"):
            out.write(f"{key}\t{getattr(args, key)}\n")
        out.write(f"left_anchor_unique\t{unique[locus.name + ':L']}\n")
        out.write(f"right_anchor_unique\t{unique[locus.name + ':R']}\n")
    left = reference.fetch(locus.chrom, locus.start - 25, locus.start)
    right = reference.fetch(locus.chrom, locus.end, locus.end + 25)
    bams = discover_bams(args.bam_dir, locus.name)
    reads = []
    failures = []
    per_bam = []
    if not bams:
        failures.append("没有与该位点匹配的 BAM")
    if unique[locus.name + ':L'] and unique[locus.name + ':R']:
        for sample, platform, bam in bams:
            try:
                extracted, audit = extract_from_bam(bam, sample, platform, locus, reference,
                                                     str(args.samtools), args.anchor_length,
                                                     args.min_mapq)
                reads.extend(extracted)
                per_bam.append((str(bam), len(extracted), dict(audit)))
            except (ValueError, RuntimeError) as error:
                failures.append(f"{bam.name}: {error}")
    else:
        failures.append("左右参考锚点没有通过全基因组唯一性检查")
    reads = list({(r.sample, r.platform, r.name): r for r in reads}.values())
    alleles = cluster_reads(reads, args.min_support)
    if alleles:
        ref_sequence = reference.fetch(locus.chrom, locus.start, locus.end)
        motifs = discover_motifs([ref_sequence] + [a.sequence for a in alleles], locus.motifs)
        for allele in alleles:
            allele.tokens, explained = decompose(allele.sequence, motifs)
            allele.dbg_k, allele.graph_support, allele.unitigs = dbg_evidence(allele, reads)
            allele.uncertain |= (explained < 1.0 or allele.graph_support < .90
                                 or allele.identity < .90)
        nodes, edges, paths = build_graph(alleles, left, right)
        write_gfa(destination / f"{locus.name}.gfa", nodes, edges, paths)
        write_svg(destination / f"{locus.name}.svg", locus.name, alleles)
    else:
        write_svg(destination / f"{locus.name}.svg", locus.name, [],
                  failures[0] if failures else "双侧翼跨越 read 支持不足，未推断等位序列")
    with (destination / "audit.tsv").open("w", encoding="utf-8") as out:
        out.write("allele\tsample\tsupport\tidentity\tdbg_k\tgraph_support\tunitigs\tuncertain\tread_ids\n")
        for allele in alleles:
            out.write(f"{allele.name}\t{allele.sample}\t{allele.support}\t{allele.identity:.4f}\t"
                      f"{allele.dbg_k}\t{allele.graph_support:.4f}\t{allele.unitigs}\t"
                      f"{allele.uncertain}\t{','.join(allele.read_ids)}\n")
        for failure in failures:
            out.write(f"ERROR\t{failure}\t\t\t\t\t\t\t\n")
    with (destination / "bam_audit.tsv").open("w", encoding="utf-8") as out:
        out.write("bam\tunique_spanning_reads\treasons\n")
        for bam, count, audit in per_bam:
            out.write(f"{bam}\t{count}\t{audit}\n")
    print(f"{locus.name}: {len(bams)} BAM, {len(reads)} spanning reads, {len(alleles)} alleles")


def main() -> int:
    parser = argparse.ArgumentParser(description="VNTR 位点 GFA/SVG 生成")
    parser.add_argument("run", nargs="?", default="run")
    parser.add_argument("--bam-dir", type=Path, default=Path("bamdata"))
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--locus", action="append", default=[])
    parser.add_argument("--samtools", default=shutil.which("samtools") or "samtools")
    parser.add_argument("--min-support", type=int, default=3)
    parser.add_argument("--anchor-length", type=int, default=150)
    parser.add_argument("--min-mapq", type=int, default=20)
    args = parser.parse_args()
    if args.min_support < 1 or args.anchor_length < 100 or args.anchor_length > 300 or args.min_mapq < 0:
        parser.error("min-support 必须 >=1，anchor-length 必须在 100–300 之间，min-mapq 必须 >=0")
    if not Path(args.samtools).exists() and not shutil.which(str(args.samtools)):
        parser.error("找不到 samtools")
    loci = read_catalog(args.catalog)
    if args.locus:
        loci = [x for x in loci if x.name in args.locus]
        if len(loci) != len(set(args.locus)):
            parser.error("--locus 在 BED 中不存在")
    args.out.mkdir(parents=True, exist_ok=True)
    with IndexedFastaReference(args.reference) as reference:
        for locus in loci:
            if locus.chrom not in reference.entries:
                raise ValueError(f"参考基因组缺少 {locus.chrom}")
            if locus.start < args.anchor_length or locus.end + args.anchor_length > reference.contig_length(locus.chrom):
                raise ValueError(f"位点侧翼超出参考边界: {locus.name}")
        unique = global_unique_probes(args.reference, loci, reference, args.anchor_length)
        for locus in loci:
            process_locus(locus, args, reference, unique)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
