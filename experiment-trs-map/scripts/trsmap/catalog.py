"""解析并校验与 TRGT 兼容的串联重复位点目录。"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


# 支持 IUPAC 简并碱基，而不只限制 A/C/G/T。例如 N 表示任意碱基，
# R 表示 A/G。目录校验只检查字母是否合法，不在这里展开简并碱基。
IUPAC_DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)
# 文件名仅保留字母、数字、下划线、点和连字符；其余字符统一替换。
SAFE_ID_RE = re.compile(r"[^A-Za-z0-9_.-]+")


@dataclass(frozen=True)
class Locus:
    """一个串联重复位点，内部统一使用 BED 的 0-based、half-open 坐标。"""

    chrom: str
    start0: int
    end0: int
    locus_id: str
    motifs: tuple[str, ...]
    structure: str = "<TR>"
    metadata: dict[str, str] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        # 在数据模型入口处集中校验，使后续算法可以假定位点一定合法，
        # 避免在 read 循环中反复处理同样的异常分支。
        if not self.chrom:
            raise ValueError("chromosome must not be empty")
        if self.start0 < 0 or self.end0 <= self.start0:
            raise ValueError(
                f"invalid BED interval for {self.locus_id}: "
                f"{self.chrom}:{self.start0}-{self.end0}"
            )
        if not self.locus_id:
            raise ValueError("repeat ID must not be empty")
        if not self.motifs:
            raise ValueError(f"repeat {self.locus_id} has no MOTIFS")
        for motif in self.motifs:
            if not motif or not IUPAC_DNA_RE.fullmatch(motif):
                raise ValueError(
                    f"repeat {self.locus_id} has invalid IUPAC motif: {motif!r}"
                )

    @property
    def start1(self) -> int:
        """返回 SAM 和 samtools 区域表达式使用的 1-based 闭区间起点。"""

        return self.start0 + 1

    @property
    def end1(self) -> int:
        """返回 SAM 和 samtools 区域表达式使用的 1-based 闭区间终点。"""

        return self.end0


def parse_info(text: str, line_number: int) -> dict[str, str]:
    """解析 BED 第四列的 ``KEY=VALUE;KEY=VALUE`` 结构。"""

    info: dict[str, str] = {}
    # 空项目允许出现在末尾；非空项目必须严格使用 KEY=VALUE 格式。
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(
                f"catalog line {line_number}: expected KEY=VALUE, found {item!r}"
            )
        key, value = item.split("=", 1)
        key = key.strip().upper()
        value = value.strip()
        if not key or not value:
            raise ValueError(f"catalog line {line_number}: empty key or value")
        if key in info:
            raise ValueError(f"catalog line {line_number}: duplicate field {key}")
        info[key] = value
    return info


def parse_repeat_catalog(path: str | Path) -> list[Locus]:
    """解析 TRGT 风格的四列 BED，并检查 ID 唯一性和 motif 合法性。"""

    catalog_path = Path(path)
    loci: list[Locus] = []
    seen_ids: set[str] = set()

    with catalog_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            # BED 中的空行和注释行不参与位点解析。
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    f"catalog line {line_number}: expected at least 4 tab-separated columns"
                )
            chrom, start_text, end_text, info_text = fields[:4]
            try:
                start0 = int(start_text)
                end0 = int(end_text)
            except ValueError as error:
                raise ValueError(
                    f"catalog line {line_number}: BED start/end must be integers"
                ) from error

            info = parse_info(info_text, line_number)
            # ID 用于跨文件关联，MOTIFS 用于后续 read 证据匹配，二者不可缺少。
            missing = {"ID", "MOTIFS"} - info.keys()
            if missing:
                raise ValueError(
                    f"catalog line {line_number}: missing {', '.join(sorted(missing))}"
                )
            locus_id = info["ID"]
            if locus_id in seen_ids:
                raise ValueError(
                    f"catalog line {line_number}: duplicate repeat ID {locus_id!r}"
                )
            # dict.fromkeys 在保留用户输入顺序的同时去重。顺序需要稳定，
            # 因为将来复杂 STRUC 解析和结果复现都可能依赖目录顺序。
            motifs = tuple(
                dict.fromkeys(m.strip().upper() for m in info["MOTIFS"].split(","))
            )
            locus = Locus(
                chrom=chrom,
                start0=start0,
                end0=end0,
                locus_id=locus_id,
                motifs=motifs,
                structure=info.get("STRUC", "<TR>"),
                metadata={
                    key: value
                    for key, value in info.items()
                    if key not in {"ID", "MOTIFS", "STRUC"}
                },
            )
            loci.append(locus)
            seen_ids.add(locus_id)

    if not loci:
        raise ValueError(f"repeat catalog is empty: {catalog_path}")
    return loci


def safe_locus_id(locus_id: str) -> str:
    """把目录 ID 转换成安全文件名；目录中的原始 ID 不会被修改。"""

    safe = SAFE_ID_RE.sub("_", locus_id).strip("._")
    if not safe:
        raise ValueError(f"repeat ID cannot be converted to a filename: {locus_id!r}")
    return safe
