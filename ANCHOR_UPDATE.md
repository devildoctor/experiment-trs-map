# 点位输入整理与锚点实现（2026-10-19）

## 输入与目录

`bamdata_by_locus/点位/样本/平台/*.bam` 是当前可用的点位 BAM；`data_by_locus` 是兼容旧路径的目录链接。现有数据按 10 个点位分目录存放，当前目录中没有 `inputs.tsv`、`locus.bed` 或总清单 `inventory.tsv`；位点定义由 `catalogs/` 中的 BED 提供。不裁剪、不重比对本次输入数据。缺少的平台文件不会补造。

每个 BAM 已经是用户提供的点位区间加左右各 500 bp 的数据集。程序使用 `samtools view -h BAM` 顺序读取其全部记录；不再增加区域过滤，不访问原始全基因组 BAM。这个 500 bp 限制是参考区间取数范围，不是对已保存的 read SEQ 截断。

点位目录沿用原有 BED 定义；D1S80 为 chr1，候选窗口 [2459200,2459800)，目录标注的重复区间 [2459279,2459743)，motif 种子 CCTGGTGGTCCTCCTT。`confirmed` 是继承的目录标注，本次没有重新验证生物学边界。

## 主点位图流程

1. 全参考扫描每侧 150、250、350、500 bp 模板中的全部规范化 31-mer，包含正反链及所有 contig。唯一 31-mer 的起点至少相距 31 bp 才算不重叠；至少两个为 preferred，仅一个为 limited，没有则 nonunique。左右模板独立选择，内侧端点不移动。完整支持要求两侧实际比对各覆盖至少两个不重叠的唯一片段；并不要求这些 31-mer 在 read 上精确匹配。
2. 所有输入文件按 sample、namespace（显式清单，否则 LB）、QNAME 汇总。同一 read 的原始方向坐标依据 FLAG、hard clip 计算。不同记录的重叠序列或总长度冲突时保留为 unresolved；没有完整中间序列时不拼接。
3. 对完整模板做 infix Levenshtein 搜索，枚举 read 上的局部最优端点及同分位置，再用 Edlib 回溯。替换、插入和缺失都计入编辑代价。精确命中不会跳过其他近似位置，也不依赖三个固定的 20-mer 种子。超过 100 个候选位置直接标记 unresolved。
4. 匹配相似度为 M/(M+X+I+D)：HiFi >=0.90，ONT/ONTUL >=0.80；模板实际观测比例 (M+X)/模板长度 >=0.90。每个用于识别唯一片段的 31-mer 至少观测 28 个模板位置。内侧端点落在比对缺口时保留为边界不确定。
5. 定位需要 MAPQ >=20 且不为 255 的主比对或补充比对；模板至少 80% 的参考碱基与该比对在原始 read 坐标上相差不超过 50 bp。低 MAPQ、次级和补充记录均保留审计。次级比对不能单独提供可信定位，但同一 query 片段的竞争定位会阻止完整支持。补充比对只贡献定位证据，不贡献可拼接的新 read。
6. 每个模板长度组合联合比较左右锚点的编辑代价。切点同分歧义保留；不同位置（相差 >10 bp）的次优对与最优对差距必须达到 max(3,ceil(0.01*(左长+右长)))。合格的不同长度模板若给出不同切点，则不自动选最长模板。切点一致时选择长度和最短的组合。
7. 提取提供序列的那条记录的实际 `SEQ[qL:qR]`。正常目标比对的 BAM SEQ 不因为 FLAG 0x10 再反向互补。N 在非截取部分不导致整条 read 丢弃；重复区含 N 则不作为完整碱基支持。零长度候选单列，不强行当作普通 VNTR 序列。
8. 相同 molecule 最多计一份完整支持，再进入原有候选生成、联合分配、motif 标注与 GFA/SVG 输出；min_support=3，ploidy=2。已有候选和共识算法未在本次重新设计。

上述相似度、MAPQ、比例与容差是工程默认值，不是经过校准的错误概率。锚点上下文可靠也不等同于每个 motif 注释或二倍体基因型已得到验证。

## 可追溯输出

- `anchor_reference.json`：所有长度模板、唯一片段坐标、分级及每个探针的全参考计数/contig。每个探针最多展示前 100 个出现位置，计数及 contig 计数不截断，并记录 positions_truncated。
- `anchor_decisions.jsonl`：每个 molecule 的全部来源、每个序列提供记录、模板组合、切点、相似度、定位来源及失败原因。
- `bam_audit.tsv`：原始记录计数、保留的标记类型、各 BAM 涉及的分子决策；不同 BAM 可能贡献同一个 molecule，不能直接把这些列相加当作独立支持。
- `molecules.tsv`、`dedup.tsv`：实际完整序列及去重；`assignments.tsv`、`genotypes.tsv`、`candidates.tsv` 和 `audit.tsv` 保留原候选流程的证据。
- GFA 中两端的 25 bp 节点属于参考展示上下文，不是 150–500 bp 的定位模板。

## 运行与兼容

计算通过 Slurm 执行。主程序新增 `vntranchors.py`，依赖固定在 `requirements-anchors.txt`，部署到项目 `.deps/vntr-anchors`；Slurm 脚本设置 PYTHONPATH。`--anchor-length` 现在表示最大模板长度，默认为 500，保留不超过它的各档模板。

主点位图使用上述严格流程。`--read-view` / `--structure-view` 的展示和分组逻辑未改为主图的支持判定，不能把这些展示行数与完整分子支持混用。

## D1S80 的 ALT 对照模式

本次全参考扫描发现，左右 500 bp 的 940 个候选 31-mer 各出现两次：一次在 chr1，一次在 chr1_KI270762v1_alt。严格 FASTA 唯一性结果原样保留。

[NCBI/GRC](https://www.ncbi.nlm.nih.gov/grc/human) 把 KI270762.1 放置在 chr1 的 2,448,811–2,791,270 区域（1-based），覆盖 D1S80。`catalogs/grch38_d1s80_alt_equivalence_v1.tsv` 仅声明本次两侧各 500 bp 的一一对应区段。启用 `--reference-equivalence` 后，运行时逐碱基验证每个区段完全相同；任一差异都会报错，不会按 contig 名称自动忽略 ALT。

原始 `count` 与 `raw_unique_count` 保留；`locus_count` 和模板 `unique_count` 按这张显式映射表折算。其他位置出现的同一探针仍算重复。此模式不修改 MAPQ，也不恢复提供的 BAM 之外的记录，不把所有 ALT contig 都视为可忽略。Slurm 可通过 `REFERENCE_EQUIVALENCE` 显式启用。
