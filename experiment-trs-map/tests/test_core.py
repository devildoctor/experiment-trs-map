from __future__ import annotations

import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
import sys

# 项目未安装成系统包时，把 scripts 加入搜索路径以便测试直接导入源码。
sys.path.insert(0, str(SCRIPTS))

from trsmap.catalog import parse_repeat_catalog  # noqa: E402
from trsmap.cigar import (  # noqa: E402
    query_interval_for_reference,
    reference_end,
    soft_clip_lengths,
)
from trsmap.evidence import (  # noqa: E402
    SelectionConfig,
    evaluate_alignment,
    motif_coverage,
)
from trsmap.locus_sequence import extract_locus_slice  # noqa: E402


class CatalogTests(unittest.TestCase):
    """验证 TRGT 风格目录的坐标转换、元数据解析和错误检测。"""

    def test_parses_trgt_bed_coordinates_and_metadata(self) -> None:
        content = (
            "# comment\n"
            "chr19\t45770204\t45770264\t"
            "ID=DMPK;MOTIFS=CTG,CAG;STRUC=<TR>;GENE=DMPK\n"
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "repeats.bed"
            path.write_text(content, encoding="utf-8")
            loci = parse_repeat_catalog(path)
        self.assertEqual(len(loci), 1)
        self.assertEqual(loci[0].start1, 45770205)
        self.assertEqual(loci[0].end1, 45770264)
        self.assertEqual(loci[0].motifs, ("CTG", "CAG"))
        self.assertEqual(loci[0].metadata["GENE"], "DMPK")

    def test_rejects_duplicate_ids(self) -> None:
        content = (
            "chr1\t10\t20\tID=x;MOTIFS=CAG;STRUC=<TR>\n"
            "chr2\t10\t20\tID=x;MOTIFS=CTG;STRUC=<TR>\n"
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "repeats.bed"
            path.write_text(content, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate repeat ID"):
                parse_repeat_catalog(path)


class CigarTests(unittest.TestCase):
    """验证 CIGAR 的参考终点、插入保留和软剪切边界处理。"""

    def test_reference_end(self) -> None:
        self.assertEqual(reference_end(100, "10M3I5M2D4M"), 120)

    def test_query_interval_includes_internal_insertion(self) -> None:
        interval = query_interval_for_reference(100, "10M6I10M", 105, 114)
        self.assertEqual(interval, (5, 21))

    def test_soft_clips_with_hard_clips(self) -> None:
        self.assertEqual(soft_clip_lengths("5H7S100M9S3H"), (7, 9))


class EvidenceTests(unittest.TestCase):
    """验证 motif 评分、双侧锚定、单侧救回和旧版 NM 过滤。"""

    def setUp(self) -> None:
        self.locus = parse_repeat_catalog_from_text(
            "chr1\t1000\t1012\tID=test;MOTIFS=CAG;STRUC=<TR>\n"
        )

    def test_motif_coverage_handles_rotation_and_interruption(self) -> None:
        perfect, motif = motif_coverage("AGCAGCAGCAGC", ("CAG",))
        interrupted, _ = motif_coverage("CAGCAGTAACAGCAG", ("CAG",))
        random, _ = motif_coverage("ATTTACGGATTA", ("CAG",))
        self.assertEqual(perfect, 1.0)
        self.assertEqual(motif, "CAG")
        self.assertGreater(interrupted, 0.70)
        self.assertLess(random, 0.50)

    def test_accepts_fully_spanning_read_without_motif_gate(self) -> None:
        evidence = evaluate_alignment(
            pos=900,
            cigar="213M",
            mapq=60,
            sequence="A" * 213,
            nm=20,
            locus=self.locus,
            config=SelectionConfig(anchor_length=100),
        )
        self.assertTrue(evidence.accepted)
        self.assertEqual(evidence.evidence_class, "spanning")

    def test_rescues_one_sided_read_with_motif_support(self) -> None:
        sequence = "A" * 101 + "CAG" * 4 + "T" * 60
        evidence = evaluate_alignment(
            pos=900,
            cigar="113M60S",
            mapq=60,
            sequence=sequence,
            nm=15,
            locus=self.locus,
            config=SelectionConfig(anchor_length=100),
        )
        self.assertTrue(evidence.accepted)
        self.assertEqual(evidence.evidence_class, "left_anchor_motif")
        self.assertEqual(evidence.repeat_query_length, 12)

    def test_legacy_nm_filter_is_optional(self) -> None:
        base_args = dict(
            pos=900,
            cigar="213M",
            mapq=60,
            sequence="A" * 213,
            nm=20,
            locus=self.locus,
        )
        accepted = evaluate_alignment(
            **base_args, config=SelectionConfig(anchor_length=100)
        )
        rejected = evaluate_alignment(
            **base_args,
            config=SelectionConfig(anchor_length=100, max_nm_fraction=0.01),
        )
        self.assertTrue(accepted.accepted)
        self.assertFalse(rejected.accepted)
        self.assertEqual(rejected.reason, "nm_fraction")


class LocusSequenceTests(unittest.TestCase):
    """验证局部序列切片会保留真实插入，且不会人工填充缺失侧翼。"""

    def setUp(self) -> None:
        self.locus = parse_repeat_catalog_from_text(
            "chr1\t1000\t1012\tID=test;MOTIFS=CAG;STRUC=<TR>\n"
        )

    def test_extracts_left_repeat_right_and_internal_insertion(self) -> None:
        # 参考位置 991-1022 包含 10 bp 左侧翼、12 bp 重复区和 10 bp 右侧翼；
        # 位于重复区内部的 6 bp 插入必须保留在 repeat_sequence 中。
        sequence = "A" * 10 + "CAGCA" + "G" * 6 + "GCAGCAG" + "T" * 10
        evidence = evaluate_alignment(
            pos=991,
            cigar="15M6I17M",
            mapq=60,
            sequence=sequence,
            nm=6,
            locus=self.locus,
            config=SelectionConfig(anchor_length=10),
        )
        locus_slice = extract_locus_slice(
            pos=991,
            cigar="15M6I17M",
            sequence=sequence,
            quality="I" * len(sequence),
            locus=self.locus,
            evidence=evidence,
            flank_length=10,
        )

        self.assertIsNotNone(locus_slice)
        assert locus_slice is not None
        self.assertEqual(locus_slice.left_sequence, "A" * 10)
        self.assertEqual(locus_slice.right_sequence, "T" * 10)
        self.assertEqual(len(locus_slice.repeat_sequence), 18)
        self.assertEqual(
            locus_slice.sequence,
            locus_slice.left_sequence
            + locus_slice.repeat_sequence
            + locus_slice.right_sequence,
        )
        self.assertTrue(locus_slice.left_complete)
        self.assertTrue(locus_slice.right_complete)

    def test_partial_read_is_not_padded_with_artificial_bases(self) -> None:
        sequence = "A" * 10 + "CAG" * 4 + "T" * 10
        evidence = evaluate_alignment(
            pos=991,
            cigar="32M",
            mapq=60,
            sequence=sequence,
            nm=0,
            locus=self.locus,
            config=SelectionConfig(anchor_length=10),
        )
        locus_slice = extract_locus_slice(
            pos=991,
            cigar="32M",
            sequence=sequence,
            quality="I" * len(sequence),
            locus=self.locus,
            evidence=evidence,
            flank_length=100,
        )

        self.assertIsNotNone(locus_slice)
        assert locus_slice is not None
        self.assertFalse(locus_slice.left_complete)
        self.assertFalse(locus_slice.right_complete)
        self.assertNotIn("N", locus_slice.sequence)


def parse_repeat_catalog_from_text(content: str):
    """借助临时 BED 文件，把一段目录文本解析为单个测试位点。"""

    with tempfile.TemporaryDirectory() as temp_dir:
        path = Path(temp_dir) / "one.bed"
        path.write_text(content, encoding="utf-8")
        return parse_repeat_catalog(path)[0]


if __name__ == "__main__":
    unittest.main()
