#!/usr/bin/env python3
"""从一个或多个位点中筛选用于串联重复分析的长 read。

命令既接受与 TRGT 兼容的重复位点 BED，也兼容旧版单个位点参数。筛选基于
可解释证据：跨越两侧唯一侧翼的 read 直接保留；只有单侧锚定的比对若在
位点对应序列中含有预期重复 motif，也可被救回。这样可避免把真实重复扩增
误判成整条 read 的比对错误。
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from trsmap.catalog import Locus, parse_repeat_catalog, safe_locus_id
from trsmap.evidence import (
    AlignmentEvidence,
    SelectionConfig,
    evaluate_alignment,
    reverse_complement,
)
from trsmap.locus_sequence import (
    LOCUS_SEQUENCE_COLUMNS,
    extract_locus_slice,
)


VERSION = "2.0.0"


@dataclass(frozen=True)
class Candidate:
    """一条最终候选 alignment，以及选择它所需的全部证据。"""

    fields: list[str]
    sam_line: str
    evidence: AlignmentEvidence
    nm: int | None
    read_group: str
    audit_row: str

    @property
    def rank(self) -> tuple[int, int, float, int]:
        """同一 read 有多个 alignment 时使用的确定性优先级。

        优先顺序依次为：证据类型、MAPQ、motif 比例、重复 query 长度。
        tuple 的字典序保证相同输入每次选择同一个 alignment。
        """

        evidence_rank = {
            "spanning": 3,
            "left_anchor_motif": 2,
            "right_anchor_motif": 2,
        }.get(self.evidence.evidence_class, 0)
        return (
            evidence_rank,
            int(self.fields[4]),
            self.evidence.motif_fraction,
            self.evidence.repeat_query_length,
        )


def parse_args() -> argparse.Namespace:
    """定义命令行接口并解析参数；参数之间的组合约束稍后统一校验。"""

    parser = argparse.ArgumentParser(
        description=(
            "Select repeat-overlapping long reads using flank and motif evidence. "
            "Catalog coordinates follow BED/TRGT convention (0-based, half-open); "
            "--start/--end follow SAM convention (1-based, inclusive)."
        )
    )
    parser.add_argument("--bam", required=True, help="Coordinate-sorted, indexed BAM")
    parser.add_argument(
        "--repeats",
        help="TRGT-compatible BED with ID, MOTIFS and optional STRUC fields",
    )
    parser.add_argument(
        "--repeat-id",
        action="append",
        default=[],
        help="Only process this catalog ID; may be supplied more than once",
    )
    parser.add_argument("--chrom", help="Single-locus reference contig")
    parser.add_argument("--start", type=int, help="Single-locus start, 1-based inclusive")
    parser.add_argument("--end", type=int, help="Single-locus end, 1-based inclusive")
    parser.add_argument(
        "--motif",
        action="append",
        default=[],
        help="Expected motif for single-locus mode; may be supplied more than once",
    )
    parser.add_argument("--locus-id", help="ID for a single locus; defaults to prefix")
    parser.add_argument(
        "--anchor-length",
        "--flank",
        dest="anchor_length",
        type=int,
        default=100,
        help="Required aligned reference bases on each side (default: 100)",
    )
    parser.add_argument(
        "--fetch-padding",
        type=int,
        default=1000,
        help="Additional reference window queried around each locus (default: 1000)",
    )
    parser.add_argument(
        "--output-flank-length",
        type=int,
        default=100,
        help=(
            "Reference bases emitted on each side of the repeat in the separate "
            "locus-sequence outputs (default: 100)"
        ),
    )
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument(
        "--min-motif-fraction",
        type=float,
        default=0.70,
        help="Minimum motif-like fraction for rescuing a one-sided read",
    )
    parser.add_argument(
        "--min-repeat-bases",
        type=int,
        default=12,
        help="Minimum locus-aligned query length used for motif rescue",
    )
    parser.add_argument(
        "--max-nm-fraction",
        type=float,
        default=None,
        help=(
            "Optional legacy whole-read NM/read-length filter. Disabled by default "
            "because true repeat variation contributes to NM."
        ),
    )
    parser.add_argument(
        "--no-partial-rescue",
        action="store_true",
        help="Keep only reads that span both anchors",
    )
    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Consider supplementary alignments and keep the best alignment per read",
    )
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--samtools",
        default=os.environ.get("SAMTOOLS", "samtools"),
        help="samtools executable",
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--prefix",
        default=None,
        help="Output prefix; locus IDs are appended when processing multiple loci",
    )
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def tag_value(fields: list[str], name: str) -> str | None:
    """从 SAM optional fields 中读取 NM、RG 等标签。"""

    prefix = name + ":"
    for field in fields[11:]:
        if field.startswith(prefix):
            parts = field.split(":", 2)
            return parts[2] if len(parts) == 3 else None
    return None


def write_fastq_record(
    handle, qname: str, flag: int, sequence: str, quality: str
) -> None:
    """按原始测序方向写入一条 FASTQ 记录。"""

    # SAM/BAM 中负链 read 以参考方向存储；FASTQ 必须恢复测序仪产生的
    # 原始方向，所以序列做反向互补，质量字符串只做反转。
    if flag & 16:
        sequence = reverse_complement(sequence)
        quality = quality[::-1] if quality != "*" else quality
    if quality == "*":
        quality = "I" * len(sequence)
    handle.write(f"@{qname}\n{sequence}\n+\n{quality}\n")


def run_checked(command: list[str], **kwargs) -> subprocess.CompletedProcess:
    """运行外部命令；返回码非零时立即抛出异常，避免继续使用残缺输出。"""

    return subprocess.run(command, check=True, text=True, **kwargs)


def resolve_samtools(value: str) -> str:
    """解析 samtools 路径，并在启动主流程前确认文件确实存在。"""

    resolved = shutil.which(value) if os.sep not in value else value
    if not resolved or not Path(resolved).is_file():
        raise SystemExit(f"samtools not found: {value}")
    return str(resolved)


def load_loci(args: argparse.Namespace) -> list[Locus]:
    """把目录模式或旧版单个位点参数统一转换成 ``list[Locus]``。"""

    using_catalog = args.repeats is not None
    using_single = any(
        value is not None for value in (args.chrom, args.start, args.end)
    ) or bool(args.motif)
    if using_catalog and using_single:
        raise SystemExit(
            "Use either --repeats or --chrom/--start/--end/--motif, not both"
        )
    if using_catalog:
        try:
            loci = parse_repeat_catalog(args.repeats)
        except ValueError as error:
            raise SystemExit(str(error)) from error
        if args.repeat_id:
            requested = set(args.repeat_id)
            available = {locus.locus_id for locus in loci}
            missing = requested - available
            if missing:
                raise SystemExit(
                    "repeat ID(s) not found: " + ", ".join(sorted(missing))
                )
            loci = [locus for locus in loci if locus.locus_id in requested]
        return loci

    missing_args = [
        name
        for name, value in (
            ("--chrom", args.chrom),
            ("--start", args.start),
            ("--end", args.end),
        )
        if value is None
    ]
    if missing_args or not args.motif:
        suffix = " and at least one --motif" if not args.motif else ""
        raise SystemExit(
            "Single-locus mode requires " + ", ".join(missing_args) + suffix
        )
    locus_id = args.locus_id or args.prefix or "candidate_repeat"
    try:
        return [
            Locus(
                chrom=args.chrom,
                start0=args.start - 1,
                end0=args.end,
                locus_id=locus_id,
                motifs=tuple(dict.fromkeys(motif.upper() for motif in args.motif)),
            )
        ]
    except ValueError as error:
        raise SystemExit(str(error)) from error


def validate_args(args: argparse.Namespace) -> None:
    """集中检查数值范围和参数间约束，尽早给出明确错误。"""

    if args.anchor_length < 0:
        raise SystemExit("--anchor-length must be non-negative")
    if args.fetch_padding < args.anchor_length:
        raise SystemExit("--fetch-padding must be at least --anchor-length")
    if args.output_flank_length < 1:
        raise SystemExit("--output-flank-length must be positive")
    if not 0 <= args.min_motif_fraction <= 1:
        raise SystemExit("--min-motif-fraction must be between 0 and 1")
    if args.min_repeat_bases < 1:
        raise SystemExit("--min-repeat-bases must be positive")
    if args.max_nm_fraction is not None and not 0 <= args.max_nm_fraction <= 1:
        raise SystemExit("--max-nm-fraction must be between 0 and 1")
    if args.threads < 1:
        raise SystemExit("--threads must be positive")


def output_prefix(args: argparse.Namespace, locus: Locus, total_loci: int) -> str:
    """保证多位点运行不会因相同全局 prefix 互相覆盖结果。"""

    locus_name = safe_locus_id(locus.locus_id)
    if args.prefix is None:
        return locus_name
    return f"{args.prefix}.{locus_name}" if total_loci > 1 else args.prefix


def format_optional(value: object | None) -> str:
    """把缺失值统一写成点号，保持 TSV 字段可稳定解析。"""

    return "." if value is None else str(value)


def process_locus(
    *,
    args: argparse.Namespace,
    locus: Locus,
    total_loci: int,
    samtools: str,
    out_dir: Path,
) -> dict[str, str | int]:
    """完成一个位点的 fetch、判定、去重、输出和一致性校验。"""

    prefix = output_prefix(args, locus, total_loci)
    bam_out = out_dir / f"{prefix}.bam"
    bai_out = Path(str(bam_out) + ".bai")
    fastq_out = out_dir / f"{prefix}.fastq.gz"
    tsv_out = out_dir / f"{prefix}.tsv"
    audit_out = out_dir / f"{prefix}.audit.tsv"
    summary_out = out_dir / f"{prefix}.summary.tsv"
    # 完整候选 BAM 作为可审计的数据源保留；以位点为中心的紧凑序列则
    # 单独写入分析文件，避免混淆“原始证据”和“下游分析切片”。
    locus_tsv_out = out_dir / f"{prefix}.locus_sequences.tsv"
    locus_fastq_out = out_dir / f"{prefix}.locus_sequences.fastq.gz"
    repeat_fasta_out = out_dir / f"{prefix}.repeat_sequences.fasta"

    config = SelectionConfig(
        anchor_length=args.anchor_length,
        min_mapq=args.min_mapq,
        min_motif_fraction=args.min_motif_fraction,
        min_repeat_bases=args.min_repeat_bases,
        max_nm_fraction=args.max_nm_fraction,
        rescue_partial=not args.no_partial_rescue,
    )
    fetch_start = max(1, locus.start1 - args.fetch_padding)
    fetch_end = locus.end1 + args.fetch_padding
    fetch_region = f"{locus.chrom}:{fetch_start}-{fetch_end}"
    # SAM flag 位掩码：始终排除 unmapped(0x4)、secondary(0x100)、
    # QC-fail(0x200)、duplicate(0x400)。默认也排除 supplementary(0x800)；
    # 开启 --include-supplementary 后保留它，以发现被拆分的大扩增 alignment。
    exclude_flags = 1796 if args.include_supplementary else 3844
    view_command = [
        samtools,
        "view",
        "-h",
        "-F",
        str(exclude_flags),
        str(args.bam),
        fetch_region,
    ]

    counts: Counter[str] = Counter()
    headers: list[str] = []
    candidates: dict[str, Candidate] = {}
    audit_rows: list[str] = []

    process = subprocess.Popen(
        view_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    assert process.stderr is not None

    audit_header = (
        "read_name\tlocus_id\tchrom\tref_start\tref_end\tmapq\tflag\tstrand\t"
        "evidence_class\taccepted\treason\tleft_anchor_bases\t"
        "right_anchor_bases\trepeat_query_start\trepeat_query_end\t"
        "repeat_query_length\tmotif_fraction\tbest_motif\tleft_soft_clip\t"
        "right_soft_clip\tread_length\tNM\tnm_fraction\tread_group"
    )

    # 流式读取区域 SAM，避免把 BAM 或整个区域的序列一次性装入内存。
    for line in process.stdout:
        if line.startswith("@"):
            headers.append(line)
            continue
        counts["region_alignments"] += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            counts["rejected_malformed"] += 1
            continue
        try:
            qname = fields[0]
            flag = int(fields[1])
            pos = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            sequence = fields[9]
            quality = fields[10]
            nm_text = tag_value(fields, "NM")
            nm = int(nm_text) if nm_text is not None else None
            evidence = evaluate_alignment(
                pos=pos,
                cigar=cigar,
                mapq=mapq,
                sequence=sequence,
                nm=nm,
                locus=locus,
                config=config,
            )
        except (ValueError, IndexError):
            counts["rejected_malformed"] += 1
            continue

        read_length = 0 if sequence == "*" else len(sequence)
        nm_fraction = (
            nm / read_length if nm is not None and read_length > 0 else None
        )
        read_group = tag_value(fields, "RG") or "."
        strand = "-" if flag & 16 else "+"
        audit_row = "\t".join(
                map(
                    str,
                    (
                        qname,
                        locus.locus_id,
                        fields[2],
                        pos,
                        evidence.ref_end,
                        mapq,
                        flag,
                        strand,
                        evidence.evidence_class,
                        int(evidence.accepted),
                        evidence.reason,
                        evidence.left_anchor_bases,
                        evidence.right_anchor_bases,
                        format_optional(evidence.repeat_query_start),
                        format_optional(evidence.repeat_query_end),
                        evidence.repeat_query_length,
                        f"{evidence.motif_fraction:.6f}",
                        evidence.best_motif,
                        evidence.left_soft_clip,
                        evidence.right_soft_clip,
                        read_length,
                        format_optional(nm),
                        "." if nm_fraction is None else f"{nm_fraction:.6f}",
                        read_group,
                    ),
                )
        )
        # 无论接受还是拒绝都进入审计表；这使参数调优可以直接比较每条 read
        # 的中间量，而不必反复从 BAM 重建判定过程。
        audit_rows.append(audit_row)
        counts[f"alignment_{evidence.reason}"] += 1
        if not evidence.accepted or sequence == "*":
            continue
        candidate = Candidate(fields, line, evidence, nm, read_group, audit_row)
        # 一个 read 可能有多个 supplementary alignment。最终 BAM/FASTQ 只
        # 保留证据最强的一条，避免重复计数和等位基因聚类偏倚。
        previous = candidates.get(qname)
        if previous is None or candidate.rank > previous.rank:
            candidates[qname] = candidate
            if previous is not None:
                counts["duplicate_alignment_replaced"] += 1
        else:
            counts["duplicate_alignment_ignored"] += 1

    stderr = process.stderr.read()
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"samtools view failed ({return_code}):\n{stderr}")

    counts["candidate_reads"] = len(candidates)
    for candidate in candidates.values():
        counts[f"candidate_{candidate.evidence.evidence_class}"] += 1

    # 先写临时 SAM/BAM，再排序并原子地生成最终文件。临时目录正常结束后
    # 自动清理，失败时不会把半成品误当成有效结果。
    with tempfile.TemporaryDirectory(prefix="trs_select_", dir=out_dir) as temp_dir:
        sam_path = Path(temp_dir) / f"{prefix}.sam"
        unsorted_bam = Path(temp_dir) / f"{prefix}.unsorted.bam"
        with (
            sam_path.open("w", encoding="utf-8") as sam_handle,
            gzip.open(fastq_out, "wt", encoding="utf-8") as fastq_handle,
            tsv_out.open("w", encoding="utf-8") as tsv_handle,
            locus_tsv_out.open("w", encoding="utf-8") as locus_tsv_handle,
            gzip.open(locus_fastq_out, "wt", encoding="utf-8") as locus_fastq_handle,
            repeat_fasta_out.open("w", encoding="utf-8") as repeat_fasta_handle,
        ):
            sam_handle.writelines(headers)
            tsv_handle.write(audit_header + "\n")
            locus_tsv_handle.write("\t".join(LOCUS_SEQUENCE_COLUMNS) + "\n")

            for candidate in candidates.values():
                # 完整比对写入 BAM 作为审计来源；下方的独立文件仅保存
                # “左侧翼 + 重复区 + 右侧翼”，供等位基因分类使用。
                sam_handle.write(candidate.sam_line)
                write_fastq_record(
                    fastq_handle,
                    candidate.fields[0],
                    int(candidate.fields[1]),
                    candidate.fields[9],
                    candidate.fields[10],
                )
                tsv_handle.write(candidate.audit_row + "\n")

                qname = candidate.fields[0]
                strand = "-" if int(candidate.fields[1]) & 16 else "+"
                locus_slice = extract_locus_slice(
                    pos=int(candidate.fields[3]),
                    cigar=candidate.fields[5],
                    sequence=candidate.fields[9],
                    quality=candidate.fields[10],
                    locus=locus,
                    evidence=candidate.evidence,
                    flank_length=args.output_flank_length,
                )
                if locus_slice is None:
                    # 不凭空补造缺失碱基；投影失败会计入汇总，供后续检查。
                    counts["locus_sequence_unavailable"] += 1
                    continue

                row = (
                    qname,
                    locus.locus_id,
                    candidate.evidence.evidence_class,
                    strand,
                    "reference",
                    locus_slice.requested_flank_length,
                    int(locus_slice.left_complete),
                    int(locus_slice.right_complete),
                    locus_slice.query_start,
                    locus_slice.query_end,
                    locus_slice.repeat_start_in_slice,
                    locus_slice.repeat_end_in_slice,
                    len(locus_slice.left_sequence),
                    len(locus_slice.repeat_sequence),
                    len(locus_slice.right_sequence),
                    locus_slice.left_sequence,
                    locus_slice.repeat_sequence,
                    locus_slice.right_sequence,
                    locus_slice.sequence,
                )
                locus_tsv_handle.write("\t".join(map(str, row)) + "\n")

                # 记录名后缀用于区分局部片段和整条 read 的 FASTQ；无论 SAM
                # 中是哪条链，这里的局部序列都按参考方向输出。
                record_name = (
                    f"{qname}|locus={locus.locus_id}"
                    f"|flank={args.output_flank_length}"
                    f"|left_complete={int(locus_slice.left_complete)}"
                    f"|right_complete={int(locus_slice.right_complete)}"
                )
                locus_fastq_handle.write(
                    f"@{record_name}\n{locus_slice.sequence}\n+\n"
                    f"{locus_slice.quality}\n"
                )
                repeat_fasta_handle.write(
                    f">{qname}|locus={locus.locus_id}|orientation=reference\n"
                    f"{locus_slice.repeat_sequence}\n"
                )
                counts["locus_sequence_records"] += 1
                counts["repeat_sequence_records"] += 1

        # 即使输入 BAM 是坐标排序的，按 read 去重和替换后顺序也可能变化，
        # 因此显式 sort，随后建立索引并用 quickcheck 检查文件完整性。
        run_checked(
            [samtools, "view", "-b", "-o", str(unsorted_bam), str(sam_path)]
        )
        run_checked(
            [
                samtools,
                "sort",
                "-@",
                str(args.threads),
                "-o",
                str(bam_out),
                str(unsorted_bam),
            ]
        )
        run_checked([samtools, "index", "-@", str(args.threads), str(bam_out)])
        run_checked([samtools, "quickcheck", "-v", str(bam_out)])

    with audit_out.open("w", encoding="utf-8") as handle:
        handle.write(audit_header + "\n")
        handle.write("\n".join(audit_rows))
        if audit_rows:
            handle.write("\n")

    bam_count = int(
        run_checked(
            [samtools, "view", "-c", str(bam_out)], capture_output=True
        ).stdout.strip()
    )
    # 最终硬校验：内存中的候选数必须等于 BAM alignment 数。FASTQ 和 TSV
    # 都从同一个 candidates 字典生成，所以三种输出共享同一 read 集合。
    if bam_count != len(candidates):
        raise RuntimeError(
            f"output validation failed: BAM has {bam_count} reads, "
            f"expected {len(candidates)}"
        )

    with summary_out.open("w", encoding="utf-8") as handle:
        # 先记录可复现实验的全部设置，再写入按名称排序的运行计数。
        handle.write("metric\tvalue\n")
        settings = {
            "algorithm_version": VERSION,
            "input_bam": args.bam,
            "locus_id": locus.locus_id,
            "bed_interval": f"{locus.chrom}:{locus.start0}-{locus.end0}",
            "sam_interval": f"{locus.chrom}:{locus.start1}-{locus.end1}",
            "motifs": ",".join(locus.motifs),
            "structure": locus.structure,
            "fetch_region": fetch_region,
            "anchor_length": args.anchor_length,
            "output_flank_length": args.output_flank_length,
            "min_mapq": args.min_mapq,
            "min_motif_fraction": args.min_motif_fraction,
            "min_repeat_bases": args.min_repeat_bases,
            "max_nm_fraction": format_optional(args.max_nm_fraction),
            "partial_rescue": int(not args.no_partial_rescue),
            "supplementary_alignments": int(args.include_supplementary),
        }
        for key, value in settings.items():
            handle.write(f"{key}\t{value}\n")
        for key in sorted(counts):
            handle.write(f"{key}\t{counts[key]}\n")
        handle.write(f"validated_bam_read_count\t{bam_count}\n")

    print(
        f"{locus.locus_id}: {len(candidates)} candidate reads "
        f"({counts['candidate_spanning']} spanning, "
        f"{counts['candidate_left_anchor_motif'] + counts['candidate_right_anchor_motif']} rescued)"
    )
    return {
        "locus_id": locus.locus_id,
        "candidate_reads": len(candidates),
        "bam": str(bam_out),
        "bai": str(bai_out),
        "fastq": str(fastq_out),
        "tsv": str(tsv_out),
        "audit": str(audit_out),
        "summary": str(summary_out),
        "locus_sequences_tsv": str(locus_tsv_out),
        "locus_sequences_fastq": str(locus_fastq_out),
        "repeat_sequences_fasta": str(repeat_fasta_out),
    }


def main() -> int:
    """命令入口：参数校验 -> 位点循环 -> 汇总 manifest。"""

    args = parse_args()
    validate_args(args)
    bam = Path(args.bam)
    if not bam.is_file():
        raise SystemExit(f"BAM not found: {bam}")
    samtools = resolve_samtools(args.samtools)
    loci = load_loci(args)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = out_dir / f"{args.prefix or 'repeat_selection'}.manifest.tsv"
    results = [
        process_locus(
            args=args,
            locus=locus,
            total_loci=len(loci),
            samtools=samtools,
            out_dir=out_dir,
        )
        for locus in loci
    ]
    with manifest.open("w", encoding="utf-8") as handle:
        # manifest 是一次多位点运行的总索引，便于下游逐行定位所有产物。
        columns = [
            "locus_id",
            "candidate_reads",
            "bam",
            "bai",
            "fastq",
            "tsv",
            "audit",
            "summary",
            "locus_sequences_tsv",
            "locus_sequences_fastq",
            "repeat_sequences_fasta",
        ]
        handle.write("\t".join(columns) + "\n")
        for result in results:
            handle.write("\t".join(str(result[column]) for column in columns) + "\n")
    print(f"Manifest: {manifest}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as error:
        print(f"Command failed: {error}", file=sys.stderr)
        raise
