"""处理 SAM CIGAR 字符串的轻量级工具，不依赖第三方库。"""

from __future__ import annotations

import re


# SAM CIGAR 操作：M/= /X 同时消耗参考和 query；I/S 只消耗 query；
# D/N 只消耗参考；H/P 两边都不消耗。
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING = frozenset("MDN=X")
QUERY_CONSUMING = frozenset("MIS=X")


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """把 ``10M3I5M`` 转成 ``[(10, 'M'), (3, 'I'), (5, 'M')]``。"""

    if cigar == "*":
        return []
    operations = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not operations or "".join(f"{n}{op}" for n, op in operations) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar!r}")
    return operations


def reference_end(pos: int, cigar: str) -> int:
    """根据 1-based 起点和 CIGAR 计算 1-based inclusive 参考终点。"""

    consumed = sum(n for n, op in parse_cigar(cigar) if op in REF_CONSUMING)
    return pos + consumed - 1 if consumed else pos - 1


def soft_clip_lengths(cigar: str) -> tuple[int, int]:
    """返回左右软剪切长度，并兼容最外侧同时存在 hard clip 的情况。"""

    operations = parse_cigar(cigar)
    if not operations:
        return 0, 0
    left_index = 1 if operations[0][1] == "H" and len(operations) > 1 else 0
    right_index = -2 if operations[-1][1] == "H" and len(operations) > 1 else -1
    left = operations[left_index][0] if operations[left_index][1] == "S" else 0
    right = operations[right_index][0] if operations[right_index][1] == "S" else 0
    return left, right


def query_interval_for_reference(
    pos: int, cigar: str, start: int, end: int
) -> tuple[int, int] | None:
    """返回落在指定参考区间内的 query 半开区间 ``[start, end)``。

    ``pos``、``start``、``end`` 都是 1-based inclusive。若插入的断点位于
    重复区内部，插入序列也会被纳入 query 区间；这一步是保留真实扩增
    序列的关键，不能只按参考坐标长度截取。
    """

    if end < start:
        raise ValueError("end must be greater than or equal to start")
    ref_cursor = pos
    query_cursor = 0
    query_start: int | None = None
    query_end: int | None = None

    # ref_cursor 指向当前 CIGAR 操作的第一个参考位置（1-based）；
    # query_cursor 指向当前操作的第一个 query 下标（0-based）。
    for length, op in parse_cigar(cigar):
        if op in "M=X":
            op_ref_start = ref_cursor
            op_ref_end = ref_cursor + length - 1
            overlap_start = max(start, op_ref_start)
            overlap_end = min(end, op_ref_end)
            if overlap_start <= overlap_end:
                current_start = query_cursor + overlap_start - op_ref_start
                current_end = query_cursor + overlap_end - op_ref_start + 1
                query_start = (
                    current_start if query_start is None else min(query_start, current_start)
                )
                query_end = current_end if query_end is None else max(query_end, current_end)
            ref_cursor += length
            query_cursor += length
        elif op == "I":
            # 插入发生在 ref_cursor-1 与 ref_cursor 之间。若断点位于目标区间，
            # 必须把整个插入加入重复序列，否则会系统性低估扩增等位基因。
            if start < ref_cursor <= end + 1:
                query_start = query_cursor if query_start is None else min(query_start, query_cursor)
                inserted_end = query_cursor + length
                query_end = inserted_end if query_end is None else max(query_end, inserted_end)
            query_cursor += length
        elif op in "DN":
            ref_cursor += length
        elif op == "S":
            query_cursor += length
        elif op in "HP":
            continue

    if query_start is None or query_end is None or query_end <= query_start:
        return None
    return query_start, query_end
