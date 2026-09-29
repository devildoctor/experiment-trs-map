"""对覆盖重复位点的比对进行可配置、可审计的证据评分。"""

from __future__ import annotations

from dataclasses import dataclass

from .catalog import Locus
from .cigar import query_interval_for_reference, reference_end, soft_clip_lengths


# 包含 IUPAC 简并碱基的互补表。
DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    """返回 DNA 序列的反向互补，并兼容 IUPAC 简并碱基。"""

    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def motif_orientations(motif: str) -> set[str]:
    """生成 motif 两条链方向上的全部循环移位形式。

    例如 CAG、AGC、GCA 描述的是同一周期；反向链上的 CTG、TGC、GCT
    也应视为同一个重复家族。目录仍保留用户输入的原始 motif，这里的
    等价展开只用于 read 证据匹配。
    """

    motif = motif.upper()
    variants: set[str] = set()
    for oriented in (motif, reverse_complement(motif)):
        variants.update(oriented[index:] + oriented[:index] for index in range(len(oriented)))
    return variants


def motif_coverage(sequence: str, motifs: tuple[str, ...]) -> tuple[float, str]:
    """估计序列中可由带噪声 motif 分段解释的碱基比例。

    短 motif 必须精确匹配，因为三联体若允许一个错配会失去区分度；较长
    motif 最多容忍 15% 的替换。动态规划允许跳过单个碱基，从而在测序 indel
    后恢复匹配，而不要求整条重复序列始终保持同一个相位。
    """

    sequence = sequence.upper()
    if not sequence or not motifs:
        return 0.0, "."
    # 每个目录 motif 对应一组允许匹配的方向和相位，但最终 best_motif
    # 仍报告目录中的原始名称，便于跨样本汇总。
    variants = {
        motif: sorted(motif_orientations(motif), key=lambda item: (len(item), item))
        for motif in motifs
    }
    length = len(sequence)
    # scores[i]：解释 sequence[:i] 时，能由 motif 匹配覆盖的最多碱基数。
    # paths[i]：达到该最优分数时，每种原始 motif 分别贡献多少碱基。
    # 最终 motif_fraction = best_score / len(sequence)。
    scores = [-1] * (length + 1)
    paths: list[dict[str, int]] = [{} for _ in range(length + 1)]
    scores[0] = 0

    for index in range(length):
        if scores[index] < 0:
            continue
        # 转移 1：跳过一个无法解释的碱基，得分不增加。这允许测序 indel、
        # repeat interruption 或局部结构变化后重新进入正确相位。
        if scores[index] > scores[index + 1]:
            scores[index + 1] = scores[index]
            paths[index + 1] = paths[index].copy()
        # 转移 2：尝试用任一 motif 的任一等价方向覆盖当前位置。
        for motif, rotations in variants.items():
            motif_length = len(motif)
            if index + motif_length > length:
                continue
            observed = sequence[index : index + motif_length]
            mismatches = min(
                sum(base != expected for base, expected in zip(observed, rotation))
                for rotation in rotations
            )
            # 短 motif 若允许一个错配会产生过高随机命中率，例如三联体
            # 1/3 错配几乎没有区分度，因此长度 <=6 时要求精确匹配。
            # 长 motif 允许 15% 替换误差；indel 由上面的跳过转移吸收。
            allowed = 0 if motif_length <= 6 else max(1, int(motif_length * 0.15))
            if mismatches > allowed:
                continue
            covered = motif_length - mismatches
            target = index + motif_length
            candidate = scores[index] + covered
            if candidate > scores[target]:
                scores[target] = candidate
                paths[target] = paths[index].copy()
                paths[target][motif] = paths[target].get(motif, 0) + covered

    best_index = max(range(length + 1), key=lambda item: scores[item])
    best_score = max(0, scores[best_index])
    contributions = paths[best_index]
    best_motif = max(contributions, key=contributions.get) if contributions else "."
    return best_score / length, best_motif


@dataclass(frozen=True)
class SelectionConfig:
    """所有会改变 read 选择结果的阈值，集中放置便于复现和扩展。"""

    anchor_length: int = 100
    min_mapq: int = 20
    min_motif_fraction: float = 0.70
    min_repeat_bases: int = 12
    max_nm_fraction: float | None = None
    rescue_partial: bool = True


@dataclass(frozen=True)
class AlignmentEvidence:
    """单条 alignment 的完整判定结果，而不是只有一个 True/False。"""

    accepted: bool
    evidence_class: str
    reason: str
    ref_end: int
    left_anchor_bases: int
    right_anchor_bases: int
    repeat_query_start: int | None
    repeat_query_end: int | None
    repeat_query_length: int
    motif_fraction: float
    best_motif: str
    left_soft_clip: int
    right_soft_clip: int


def evaluate_alignment(
    *,
    pos: int,
    cigar: str,
    mapq: int,
    sequence: str,
    nm: int | None,
    locus: Locus,
    config: SelectionConfig,
) -> AlignmentEvidence:
    """计算锚点和 motif 证据，并返回可审计的 read 分类。

    判定顺序是有意固定的：基础质量 -> 区间相交 -> 可选旧 NM 过滤 ->
    双侧锚点 -> 单侧锚点加 motif -> 明确拒绝原因。改变顺序可能改变
    audit.tsv 中的首要拒绝原因，因此需要通过版本号记录。
    """

    end = reference_end(pos, cigar)
    # 锚点长度按参考基因组计算。负值截断为 0，表示 alignment 没有越过
    # 对应的重复边界。它衡量的是可用于定位的侧翼，不是 read 总长度。
    left_anchor = max(0, locus.start1 - pos)
    right_anchor = max(0, end - locus.end1)
    query_interval = query_interval_for_reference(
        pos, cigar, locus.start1, locus.end1
    )
    # 通过 CIGAR 把参考重复区投影到 query，并保留区间内部的插入。
    # 对大型扩增，repeat_sequence 可能远长于参考区间。
    repeat_sequence = (
        sequence[query_interval[0] : query_interval[1]]
        if query_interval is not None and sequence != "*"
        else ""
    )
    fraction, best_motif = motif_coverage(repeat_sequence, locus.motifs)
    left_clip, right_clip = soft_clip_lengths(cigar)
    interval_start = query_interval[0] if query_interval else None
    interval_end = query_interval[1] if query_interval else None

    common = dict(
        ref_end=end,
        left_anchor_bases=left_anchor,
        right_anchor_bases=right_anchor,
        repeat_query_start=interval_start,
        repeat_query_end=interval_end,
        repeat_query_length=len(repeat_sequence),
        motif_fraction=fraction,
        best_motif=best_motif,
        left_soft_clip=left_clip,
        right_soft_clip=right_clip,
    )

    # 1. MAPQ 是 read 在全基因组中定位是否唯一的基础证据。
    if mapq < config.min_mapq:
        return AlignmentEvidence(False, "rejected", "low_mapq", **common)
    # samtools fetch 返回的是与较大检索窗口相交的 reads，因此其中一部分
    # 可能只在侧翼，完全没有覆盖重复区。
    if end < locus.start1 or pos > locus.end1:
        return AlignmentEvidence(False, "rejected", "no_repeat_overlap", **common)

    # 2. NM/read_length 只为复现旧流程而保留。默认禁用，因为 NM 同时把
    # 真实重复扩增、重复中断和测序错误计算为与参考的差异。
    if config.max_nm_fraction is not None:
        if nm is None or not sequence or sequence == "*":
            return AlignmentEvidence(False, "rejected", "missing_nm_or_sequence", **common)
        if nm / len(sequence) > config.max_nm_fraction:
            return AlignmentEvidence(False, "rejected", "nm_fraction", **common)

    left_ok = left_anchor >= config.anchor_length
    right_ok = right_anchor >= config.anchor_length
    # 3. 双侧锚点是最强证据。此时不强制 motif 阈值，因为复杂位点可能
    # 包含 interruption 或目录中尚未列出的 motif，但 read 的位置可信。
    if left_ok and right_ok:
        return AlignmentEvidence(True, "spanning", "accepted", **common)

    # 4. 只有一个侧翼可靠时，必须同时满足序列长度和 motif 比例，防止
    # 普通侧翼 read 因少量偶然 motif 命中而被救回。
    motif_ok = (
        len(repeat_sequence) >= config.min_repeat_bases
        and fraction >= config.min_motif_fraction
    )
    if config.rescue_partial and motif_ok and (left_ok or right_ok):
        evidence_class = "left_anchor_motif" if left_ok else "right_anchor_motif"
        return AlignmentEvidence(True, evidence_class, "accepted", **common)

    # 5. 为每条拒绝记录一个稳定、互斥的主原因，便于比较参数变化。
    if not left_ok and not right_ok:
        reason = "missing_both_anchors"
    elif not motif_ok:
        reason = "partial_without_motif_support"
    else:
        reason = "not_fully_spanning"
    return AlignmentEvidence(False, "rejected", reason, **common)
