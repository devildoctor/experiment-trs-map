# experiment-trs-map

从长读长 BAM 中为指定串联重复（TR）位点选择候选 reads。当前阶段只负责**可审计的候选 read 选择**，为后续局部重比对、等位基因聚类和重复结构解析提供输入。

## 当前结果基线

原始 DMPK 测试使用 `chr19:45770205-45770264`、两侧各 1,000 bp、`MAPQ >= 20` 和全 read `NM/read_length <= 0.01`：

| 指标 | 数量 |
| --- | ---: |
| 区域内主比对 | 59 |
| 候选 reads | 33 |
| 因全 read NM 比例被过滤 | 14 |
| 因未跨越两侧 1 kb 锚点被过滤 | 12 |

全 read NM 会混合侧翼测序错误、重复单元插入缺失以及真实等位基因差异，因此它不再是默认过滤条件。旧行为仍可通过 `--max-nm-fraction` 开启。

v2 在相同输入 BAM 和位点上的验证结果：59 条区域比对中保留 50 条完整跨越 reads，9 条因没有覆盖重复区本身而拒绝。相比旧结果新增 17 条，其中 14 条来自关闭全 read NM 硬过滤，3 条来自把固定 1 kb 锚点改为 100 bp 可配置锚点。这 17 条仍属于候选证据，需要在下一阶段经过局部重比对和等位基因聚类，不能直接视为已确认扩增。

## 输入格式

推荐使用与 TRGT 兼容的四列 BED。前三列为 **0-based、half-open** 坐标，第四列至少包含全局唯一的 `ID` 和逗号分隔的 `MOTIFS`；`STRUC` 可保留复杂结构描述。

```text
chr19	45770204	45770264	ID=DMPK;MOTIFS=CTG;STRUC=<TR>;GENE=DMPK
```

多个可能 motif 可以写为 `MOTIFS=CTG,CAG,CCG`。应尽量使用最小且不冗余的 motif 集合，避免同一序列有多种等价分解。

项目示例见 `catalogs/example.repeats.bed`。

## 选择逻辑

每条比对先经过基础检查，再分为以下证据类型：

1. `spanning`：比对在重复区左右都覆盖至少 `--anchor-length` 个参考碱基，直接保留。
2. `left_anchor_motif` / `right_anchor_motif`：只有一侧达到锚点长度，但重复区对应的 query 序列长度和 motif 覆盖比例达到阈值，作为可能的大扩增或分裂比对保留。
3. 其余比对拒绝，但全部写入 `*.audit.tsv`，包含明确原因和中间指标。

默认参数为 100 bp 锚点、`MAPQ >= 20`、motif 覆盖比例至少 0.70、重复证据至少 12 bp。短 motif 必须精确匹配；允许跳过少量碱基以在 HiFi indel 或重复中断后恢复相位。motif 会同时检查循环移位和反向互补形式。

## 运行

```bash
python3 scripts/extract_candidate_reads.py \
  --bam sample.bam \
  --repeats catalogs/example.repeats.bed \
  --samtools /path/to/samtools \
  --output-dir results/run1
```

只运行目录中的某个位点：

```bash
python3 scripts/extract_candidate_reads.py \
  --bam sample.bam \
  --repeats catalog.bed \
  --repeat-id DMPK \
  --output-dir results/DMPK
```

旧的单个位点命令仍可使用，但现在必须提供 motif：

```bash
python3 scripts/extract_candidate_reads.py \
  --bam sample.bam \
  --chrom chr19 --start 45770205 --end 45770264 \
  --motif CTG --locus-id DMPK \
  --output-dir results/DMPK
```

主要可调参数：

| 参数 | 作用 |
| --- | --- |
| `--anchor-length` | 每侧必须可靠覆盖的参考序列长度 |
| `--output-flank-length` | 局部序列输出使用的左右参考侧翼长度，默认固定为 100 bp |
| `--min-mapq` | 最低比对质量 |
| `--min-motif-fraction` | 单侧救回所需 motif 覆盖比例 |
| `--min-repeat-bases` | motif 证据的最短 query 序列 |
| `--no-partial-rescue` | 关闭单侧 motif 救回，只保留完整跨越 reads |
| `--include-supplementary` | 纳入 supplementary 比对，并按 read 名称选择最佳证据 |
| `--max-nm-fraction` | 可选的旧版全 read NM 比例过滤，默认关闭 |

## 输出

每个位点生成：

- `*.bam` / `*.bam.bai`：选中的最佳比对，已排序并建立索引；
- `*.fastq.gz`：候选 read 原始测序方向的序列；
- `*.tsv`：入选 read 的证据指标；
- `*.audit.tsv`：所有被检查比对及拒绝原因；
- `*.summary.tsv`：参数、计数和输出校验；
- `*.locus_sequences.tsv`：固定 100 bp 左侧翼、实际重复序列和固定 100 bp 右侧翼，并分别保存三段序列及完整性标记；
- `*.locus_sequences.fastq.gz`：参考方向的“左侧翼 + 重复区 + 右侧翼”局部序列；
- `*.repeat_sequences.fasta`：每条候选 read 的重复区序列；
- `*.manifest.tsv`：一次运行内所有位点及文件路径。

局部序列的 100 bp 是参考坐标长度。样本存在 indel 时，实际输出的 query 侧翼长度可能不是正好 100 bp；程序保留真实序列，不使用 `N` 人工补齐。单侧救回 read 通过 `left_complete` 和 `right_complete` 标明缺失侧翼。

## 代码结构和扩展点

```text
scripts/
  extract_candidate_reads.py  # CLI、samtools I/O 和结果落盘
  trsmap/
    catalog.py                 # TRGT BED 数据模型与严格校验
    cigar.py                   # 参考坐标到 query 序列的映射
    evidence.py                # motif 和锚点证据策略
    locus_sequence.py          # CIGAR 投影和固定 100 bp 局部序列提取
tests/test_core.py             # 不依赖真实 BAM 的核心测试
```

新目录字段放在 `Locus.metadata` 中，不会影响现有解析。新的 read 选择策略可以在 `evidence.py` 中实现并从 CLI 注入；BAM/FASTQ/审计输出不需要随策略重写。下一阶段建议在候选 read 选择之后增加独立的局部重比对模块：按 read 长度和序列相似度生成候选单倍型，使用 POA 生成簇共识，再把 reads 重比对到带侧翼的候选单倍型。这与 LongTR 的候选单倍型和 HMM 思路一致，同时保持第一阶段结果可单独复查。

## 文献依据

- [TRGT repeat definition](https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md)：`ID`、`MOTIFS`、`STRUC` 的目录表达及最小非冗余 motif 集合。
- [Dolzhenko et al., Nature Biotechnology, 2024](https://www.nature.com/articles/s41587-023-02057-3)：TRGT 从 HiFi reads 解析指定 TR 的共识序列、等位基因支持、甲基化和嵌合性。
- [Chiu et al., Genome Biology, 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02447-3)：Straglr 用侧翼锚定、motif 证据和分裂/软剪切比对救回扩增候选。
- [Ziaei Jam et al., Genome Biology, 2024](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03319-2)：LongTR 对重复区加上下文窗口，生成候选单倍型并局部重比对，而不是只依赖全局参考比对。
- [Weisburd et al., American Journal of Human Genetics, 2026](https://doi.org/10.1016/j.ajhg.2026.03.020)：相邻重复和周边多态性可能应作为 variation cluster 做序列级联合分析。

## 测试

```bash
python3 -m unittest discover -s tests -v
```
