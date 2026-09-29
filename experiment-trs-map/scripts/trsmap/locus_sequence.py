"""提取重复位点附近带固定参考侧翼、按参考方向排列的局部序列。

候选 BAM 有意保留完整比对，便于追溯；本模块另行生成仅包含以下部分的
紧凑分析序列：

    left flank + observed repeat sequence + right flank

通过 CIGAR 将参考坐标投影到 query 坐标，因此重复区内部的插入不会丢失，
样本特异的插入缺失也能被真实表示。
"""

from __future__ import annotations

from dataclasses import dataclass

from .catalog import Locus
from .cigar import query_interval_for_reference, reference_end
from .evidence import AlignmentEvidence


LOCUS_SEQUENCE_COLUMNS = (
    # 这些列的顺序也是 TSV 的稳定输出协议，下游脚本可按表头解析。
    "read_name",
    "locus_id",
    "evidence_class",
    "strand",
    "orientation",
    "requested_flank_length",
    "left_complete",
    "right_complete",
    "slice_query_start",
    "slice_query_end",
    "repeat_start_in_slice",
    "repeat_end_in_slice",
    "left_length",
    "repeat_length",
    "right_length",
    "left_sequence",
    "repeat_sequence",
    "right_sequence",
    "locus_sequence",
)


@dataclass(frozen=True)
class LocusSequenceSlice:
    """一条比对中实际观察到的左侧翼、重复区和右侧翼切片。"""

    requested_flank_length: int
    left_complete: bool
    right_complete: bool
    query_start: int
    query_end: int
    repeat_start_in_slice: int
    repeat_end_in_slice: int
    left_sequence: str
    repeat_sequence: str
    right_sequence: str
    sequence: str
    quality: str


def extract_locus_slice(
    *,
    pos: int,
    cigar: str,
    sequence: str,
    quality: str,
    locus: Locus,
    evidence: AlignmentEvidence,
    flank_length: int = 100,
) -> LocusSequenceSlice | None:
    """按参考方向返回“左侧翼 + 重复区 + 右侧翼”。

    请求的侧翼长度按参考坐标计算。样本存在插入或缺失时，query 中的实际
    序列长度可能不同；例如 100 bp 参考侧翼可能对应 99 或 101 个 query
    碱基。这里保留真实观察值，不通过补齐或截断制造固定长度。

    SAM 对负链比对的 SEQ 和 QUAL 以比对方向存储，因此直接切片后，正负链
    结果都能得到所需顺序：参考左侧翼、重复区、参考右侧翼。

    被单侧证据救回的 read 可能缺少一侧完整侧翼，程序仍输出已有序列，并用
    left_complete 和 right_complete 明确标记。下游构建等位基因共识时，可只
    选择两者均为真的记录作为种子 read。
    """

    if flank_length <= 0:
        raise ValueError("flank_length must be positive")
    if (
        sequence == "*"
        or evidence.repeat_query_start is None
        or evidence.repeat_query_end is None
    ):
        return None

    target_start = max(1, locus.start1 - flank_length)
    target_end = locus.end1 + flank_length
    # 目标区间按参考坐标定义，再通过 CIGAR 映射为 query 的半开区间。
    interval = query_interval_for_reference(
        pos,
        cigar,
        target_start,
        target_end,
    )
    if interval is None:
        return None

    query_start, query_end = interval
    # 防御性截断可避免异常边界产生负数或越界切片，同时不修改原始 BAM。
    query_start = max(0, min(query_start, len(sequence)))
    query_end = max(query_start, min(query_end, len(sequence)))
    repeat_start = max(
        query_start,
        min(evidence.repeat_query_start, query_end),
    )
    repeat_end = max(
        repeat_start,
        min(evidence.repeat_query_end, query_end),
    )

    left_sequence = sequence[query_start:repeat_start]
    repeat_sequence = sequence[repeat_start:repeat_end]
    right_sequence = sequence[repeat_end:query_end]
    locus_sequence = sequence[query_start:query_end]

    if quality == "*":
        # 缺失质量值时使用最高质量占位，保证 FASTQ 的序列和质量长度一致。
        locus_quality = "I" * len(locus_sequence)
    else:
        locus_quality = quality[query_start:query_end]

    alignment_end = reference_end(pos, cigar)
    return LocusSequenceSlice(
        requested_flank_length=flank_length,
        left_complete=pos <= target_start,
        right_complete=alignment_end >= target_end,
        query_start=query_start,
        query_end=query_end,
        repeat_start_in_slice=repeat_start - query_start,
        repeat_end_in_slice=repeat_end - query_start,
        left_sequence=left_sequence,
        repeat_sequence=repeat_sequence,
        right_sequence=right_sequence,
        sequence=locus_sequence,
        quality=locus_quality,
    )
