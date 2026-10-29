# VNTR 点位图程序：中文代码功能解读

> **版本说明（2026-10-30）**：下文保留2026-10-29注释交付时的源码解读、行号与缺陷复现，属于修复前记录。当前分支已修复R1–R6；新行为、兼容性、回归测试及当前函数索引见 [修复说明](CODE_FIXES_20261029.zh-CN.md)。不要把下文“本轮仅加注释”“仍未修复”等历史描述用于当前源码。

本说明对应 2026-10-29 从服务器获取的 `vntr-point-graph` 源码，基线提交为 `5be507688fe7d1f8da4c1786fd85d6df8dfccb0e`，分支为 `feature/d1s80-terminal-unit-evidence`。本轮仅向源码添加中文注释，并补充功能解读与检查报告；算法、参数、输出格式和已有科学结果不变。

这里描述的是当前实现。此前讨论的“保护长度、排列、区块及已知重要变异，允许其他小差异被代表单元替换”的有损聚合方案尚未实现。现有模板近似匹配、抽象 motif 图和显示折叠，不能当作已实现这套有损方案。

## 1. 先理解程序在解决什么问题

程序从按位点整理的长读段 BAM、参考 FASTA 和 BED 位点目录出发，提取实际序列，组织候选序列，解释重复单元，并输出可以核对来源的图及证据表。它不是一个已经经过人群准确率验证的临床分型器。

最重要的三个对象是：

| 对象 | 含义 | 可以支持什么结论 |
|---|---|---|
| REF | 参考目录区间的序列或参考上下文 | 提供参照，不增加观测分子支持 |
| OBS | 实际读取并保存的分子序列或局部片段 | 支持该分子真正覆盖的范围；局部片段不能证明整条候选 |
| CAND | 程序形成并保存的候选序列 | 是待证据检验的假设；不能给自己投票 |

节点、边、单元注释和完整路径是不同层次。一个节点有很多分子经过，不表示经过该节点的任意完整路径都成立。某个 k-mer 被同一条 read 多次访问，也不增加独立分子数。

## 2. 从输入到结果的整体流程

```text
BAM + FASTA/FAI + BED + 可选分子来源/参考等价表
  │
  ├─主点图：vntrgraph → vntranchors → vntrcore / vntrmotifs
  │    参考锚点分级 → 实际分子切片与去重 → 候选与分配
  │    → motif 字典和分解 → 明确候选路径 → GFA / SVG / TSV
  │
  ├─可选 STR 表达：vntrstr → vntrstrgraph
  │    显式 STR_STRUCTURE → 精确单元/中断表达 → 单元图及逐路径视图
  │
  ├─可选 VNTR dBG：vntrdbg_inputs → vntrdbg → vntrdbg_report
  │    冻结 REF/OBS/CAND → 连续 k+1-mer 边 → unitig 压缩
  │    → 完整证据图 / 支持图 / 逐路径视图 / 汇总报告
  │        └─D1S80 分析：vntrdbg_analysis → vntrterminal（显式开启）
  │             连续跨分支证据 + 16 bp 单元分解 + 右端模板验证
  │
  └─观测展示支线：vntrreads → vntrstructure → vntrsummary
       read 切片和状态 → 精确/抽象结构图 → 全量与概要 HTML
```

这不是每次运行都必须顺序执行的单一流水线。主点图、read-view、structure-view 和 dBG 有不同入口、输入契约和证据范围，不能把一条支线的支持定义直接套到另一条支线。

## 3. 入口与阅读顺序

| 入口 | 实际职责 | 建议先看的源码 |
|---|---|---|
| `python scripts/vntrgraph.py run ...` | 主点图：提取、候选、分配、注释、明确路径输出 | `vntrgraph.py` 的 `main`、`process_locus` |
| `python scripts/vntrgraph.py --read-view ...` | 保留并显示各类 read 覆盖与提取状态 | `vntrreads.py` |
| `python scripts/vntrgraph.py --structure-view ...` | 全量 read 的 motif 精确图与抽象结构展示 | `vntrstructure.py` |
| `python scripts/vntrgraph.py --dbg-view ...` | 由已有证据构建多 k compact dBG | `vntrdbg.py`、`vntrdbg_inputs.py` |
| `python scripts/vntrdbg_analysis.py ...` | 冻结 D1S80 输入的联合支持和基本单元分析 | `vntrdbg_analysis.py` |
| 上一入口加 `--terminal-units` | 独立分子验证右端模板；开启后 CN 全部保留 NA | `vntrterminal.py` |

`vntrstr.py`、`vntrstrgraph.py`、`vntrsummary.py` 是被其他模块调用的函数库，没有独立 `main` 或 `--help`。`vntrdbg_report.py` 则接受已有运行目录作为命令行位置参数。具备命令行入口的模块应分别查看其参数，不能把主点图的全部参数传给任意展示模式。

初读建议：先看本说明的主流程和证据定义，再看 `vntrgraph.process_locus`，随后顺着 `vntranchors`、`vntrcore`、`vntrmotifs` 阅读。熟悉主流程后，再看 STR、dBG 与 structure 支线。不要从最长的 HTML 拼接函数开始读。

## 4. 输入契约与坐标

| 输入 | 当前约定 | 常见误解 |
|---|---|---|
| BAM / BAI | 按位点、样本、平台整理；主点图读取该位点 BAM 全部记录 | 展示行数不等于完整跨越分子数 |
| FASTA / FAI | 与 BAM 参考身份及染色体命名匹配；通过索引按区间读取 | 只有序列文件名一致不证明参考相同 |
| BED | 0-based 半开区间；名称和 MOTIFS 存在第 4 列属性中 | SAM 的 POS 是 1-based，解析时才减 1 |
| REPEAT_START / REPEAT_END | 显式重复区须位于搜索窗口内部 | 搜索窗口不自动等于已确认生物学重复边界 |
| BOUNDARY=confirmed | 输入目录的边界声明 | 不是程序独立完成了实验验证 |
| sources.tsv | 主图中用 BAM 路径映射 namespace，相对路径相对该表解析 | 不同平台标签不自动代表不同分子 |
| reference-equivalence | 显式、可核对的 ALT 与主参考等价区间 | `_alt` 后缀本身不是等价证明 |

主锚定流程保存 BAM 存储 SEQ 与同记录 CIGAR 的对应关系，不能看到负链 FLAG 就再次无条件反向互补。原始 query 区间、硬剪切和方向另有记录。主/补充记录可以帮助定位，但不能拼出输入中没有完整保存的序列。

当前主图使用 `[sample, namespace, read_name]` 分子键；旧 read / structure / summary 展示支线仍有 `[sample, QNAME]` 的去重路径。后文和检查报告会明确两者差别。

## 5. 四类图不能混读

| 图 | 节点含义 | 边与路径 | 主要限制 |
|---|---|---|---|
| 默认点图 | 候选序列中不重叠的真实 token 切片，带路径上下文 | 0M 连接；两侧可有明确参考上下文 | 不是所有同名 motif 都合并成一个全局节点 |
| STR 单元图 | 显式结构中的精确单元与中断 | 可有自环；P 路径保存真实次数和顺序 | 图上任意游走不能当成已经观察到的等位序列 |
| compact dBG | 压缩的精确定向 k-mer 链，即 unitig | k−1 bp 重叠；边来自输入中连续 k+1-mer | unitig 不是生物学重复单元；节点数不是 CN |
| structure 抽象图 | motif 模板节点；unknown/partial 仍按实际切片摘要区分 | 反映抽象访问关系；精确切片另存 | `S=*` 的抽象 GFA 不能直接重建实际碱基；仅家族概览将未知/部分片段汇入 U/P 显示桶 |

支持视图过滤低支持边时必须把路径拆为片段，不能跨越被删除的区域补桥。全量证据图与过滤后的支持图应同时理解。图中存在一条连通路线，只能说明拓扑允许该路线，完整证据仍要看明确保存的路径和连续分子支持。

## 6. 参数分别在控制什么

| 参数 / 规则 | 当前典型值 | 正确解释 |
|---|---|---|
| 主图 `--min-support` | 3 | 候选或证据处理中的分子支持条件，具体看调用点 |
| `--anchor-length` | 500，上限可选 150/250/350/500 | 允许的最大侧翼模板长度；保留不超过上限的多档长度 |
| `--min-mapq` | 20 | 定位证据门槛之一；不是碱基准确率 |
| `--score-margin` | 0.5 | 描述性分配代价的最优/次优差值；不是后验概率 |
| `--ploidy` | 2，允许 1–4 | 枚举候选组合时的倍性假设；超过倍性的候选仍保留冲突状态 |
| 锚点匹配 identity | HiFi 0.90；ONT/ONTUL 0.80 | 锚点匹配规则；不是 VNTR 单元有损聚合阈值 |
| dBG 的 k | D1S80 已有 17/21/31/33/49/65 等比较入口 | 序列上下文长度；不是重复周期 |
| dBG 支持阈值 | 默认至少一个样本内达到 3 个独立分子 | 按原始节点/边判断；不是对整条路径的准确率保证 |
| `--terminal-min-support` | 默认 3，不能小于 3 | 右端模板接受条件的一部分，还需实际切片和内外上下文验证 |

阈值属于不同阶段。把这些值统一称作“90%误差标准”会混淆定位、候选分配、图过滤和单元解释。

## 7. 输出文件的阅读顺序

1. 先读 `status.json`、`manifest.json`、参数和校验信息，确定运行模式、输入身份、是否证据不足。
2. 再读 `molecules.tsv`、锚定决策和拒绝原因，判断真实观测覆盖范围。
3. 看 `candidates.tsv`、`assignments.tsv`、`genotypes.tsv`，区分候选、描述性分配及不确定组合。
4. 看 `motifs.tsv`、`tokens.tsv` 或 `primitive_units.tsv`，核对模板、实际碱基、位置、编辑操作和分解歧义。
5. 看 GFA / SVG / HTML；它们负责组织信息，不能替代前面的证据判断。
6. 对跨分支问题看 `joint_queries.tsv`、`joint_support.tsv`、`joint_molecules.tsv`，追到同一个分子的完整连续区间。

特别注意：`support_fraction = support / (support + conflict)` 排除了无信息分子。它是已建模状态中的条件比例，不是等位频率、测序正确率或基因型概率。完整单元数、长度除周期、图访问次数与最终 CN 也不能混用。

下面按模块展开当前实现。


## 8. 核心提取、候选与单元分解


审阅基线：`5be507688fe7d1f8da4c1786fd85d6df8dfccb0e`。本节仅解读现有实现，不将讨论中的有损 VNTR 聚合当作已经实现的功能。

本节函数表行号指向 `annotated/scripts/` 中加注释后的文件；问题定位另给出固定基线行号，便于与原提交核对。

### 1. 主流程与数据模型

```text
vntrgraph.main
  ├─ read_catalog：位点目录、搜索窗口、重复边界、motif/STR 配置
  ├─ IndexedFastaReference：按 FAI 随机读取参考序列
  ├─ vntranchors.load_equivalences：可选 ALT 等价表逐碱基验证
  ├─ vntranchors.reference_profiles：完整参考的侧翼唯一性检查
  └─ process_locus：逐位点处理
       ├─ discover_bams → vntranchors.extract_bams → deduplicate_reads
       ├─ call_candidates：生成候选、联合分配、残余再发现
       ├─ motif_dictionary：提出、验证及分层 motif 字典
       ├─ decompose → compound_annotations：基础与复合注释分开
       ├─ dbg_evidence：候选的局部 k-mer 辅助证据
       └─ build_graph / write_gfa / write_svg / TSV / JSON
```

`--dbg-view`、`--structure-view`、`--read-view` 是入口分派，分别交给独立模块；未指定这些参数时才执行上图的基础点图流程。因此不能用基础点图的一项规则概括所有视图。

核心对象有以下含义：

| 对象 | 主要字段 | 含义和限制 |
|---|---|---|
| `Locus`（vntrgraph.py:33） | name/chrom/start/end/motifs、search_start/end、boundary_status、str_structure | 目录声明的位点；重复边界与搜索窗口分开，`confirmed` 不等于本次运行完成实验验证 |
| `Record`（vntranchors.py:254） | SAM 标志、CIGAR、SEQ/QUAL、namespace、原始 read 长度及硬剪切 | 一条实际比对记录，可与同分子的其他记录核对位置，但不能拿各记录 SEQ 拼接 |
| `Read`（vntrcore.py:33） | sequence、query_sequence、qL/qR、anchor_scores、records/platforms、original_query_interval | 通过锚定的连续分子切片；分子键是 `[sample, namespace, read_name]` |
| `Allele`（vntrcore.py:82） | sequence、seed_molecules、read_ids、ambiguous_ids、flags、tokens | 候选序列及其分配证据，类名不表示已经确认基因型 |
| `Token`（vntrcore.py:63） | family、exact_bases、start/end、phase/orientation、edit_operations、status | 对输入序列区间的注释；`full` 表示匹配完整模板，可含小 indel，不等于与模板逐碱基相同 |

### 2. vntrgraph.py：命令入口与结果组织

| 函数 / 行号 | 作用 |
|---|---|
| `IndexedFastaReference.fetch` / 67 | 以 FAI 的 offset、line_bases、line_width 定位实际字节；输入坐标统一为 0-based 半开区间，越界报错 |
| `read_catalog` / 87 | 解析 BED 第四列的 ID、MOTIFS、REPEAT_START/END、BOUNDARY、STR_STRUCTURE；核验重复区位于搜索窗口内 |
| `discover_bams` / 122 | 搜索位点/样本/平台或旧目录布局，检查文件名中的位点、BAI 和平台标签 |
| `reference_to_query` / 139、`flank_hit` / 158 | 兼容辅助函数；前者做 CIGAR 投影，后者使用精确种子及 SequenceMatcher。当前 `process_locus` 不用这套旧锚定算法 |
| `extract_sam_records` / 207 | 兼容入口，实际转交 `vntranchors.extract_records`，并非另一套主提取流程 |
| `write_tsv` / 218、`source_hashes` / 226 | TSV 中复杂字段编码为确定性 JSON；记录各 vntr*.py 来源哈希 |
| `write_gfa` / 231 | 输出 S/L/P 及 JSON 证据；基础点图连接是 `0M`，不是 compact dBG 的 k−1 重叠 |
| `compress_tokens` / 249 | 只合并连续且 family、exact_bases、status 均相同的显示单元；用断言确认 xN 展开后仍完全恢复序列 |
| `write_svg` / 262 | 按已保存候选逐行显示；节点等宽，横向几何不是 bp 比例；U1 等标签仅在当前位点内有意义 |
| `process_locus` / 316 | 将提取、候选、单元注释、图、参数与审计表串联起来 |
| `main` / 447 | 校验命令行、读取目录、检查已有输出、防止覆盖非空位点目录，再执行参考扫描和逐位点流程 |

输入准备需要参考 FASTA+FAI、按位点整理的 BAM+BAI、BED 位点目录。可选 `--sources` 显式给每个 BAM 指定分子 namespace；相对路径按该 TSV 所在目录解析。`--reference-equivalence` 提供经过说明和序列验证的 ALT 坐标等价表。

`process_locus` 先保存 `parameters.tsv`、`anchor_decisions.jsonl`、`molecules.tsv` 和 `dedup.tsv`，再生成候选与分配结果。只有分配后独立分子数达到门槛的候选进入主图，但 `candidates.tsv` 继续保留全部候选。未进入主图者不一定已有单元注释：注释循环只遍历达标候选，因此不应把它们空的 `unit_counts` 解释为“没有重复”。

随后输出 `assignments.tsv`、`genotypes.tsv`、`motifs.tsv`、`tokens.tsv`、`compound_blocks.tsv`、`annotation_alternatives.json`、`bam_audit.tsv`、`status.json`。若目录具有显式 STR 结构，通用图移到 `.detail.*`，再由 STR 模块写专用表达和图。

### 3. vntranchors.py：锚定、定位与分子证据

#### 3.1 为什么既要序列匹配，又要比对定位

一条 read 内找到了相似侧翼，仅能说明其 SEQ 含有相似片段，不能单独证明该片段来自目标位点。当前实现同时检查：参考侧翼的区分能力、read 中侧翼匹配质量、同一分子的基因组比对定位，以及竞争位置。补充比对可以证明位置，但重复区碱基仍必须由一条记录的连续 SEQ 提供。

| 步骤 / 函数 | 实际规则 |
|---|---|
| `load_equivalences` / 50 | 检查 evidence_url、方向、ALT 与主参考逐碱基相等，以及 ALT 块不重叠；只有明确声明的等价块可归一 |
| `reference_profiles` / 89 | 为每侧 150/250/350/500 bp 模板提取 canonical 31-mer，在整个 FASTA 所有 contig 中扫描；保留完整计数和最多100个位置 |
| 参考分级 | 至少2个不重叠的唯一31-mer为 preferred，只有少量唯一证据为 limited，无唯一证据为 nonunique；可选显式 ALT 等价归一 |
| `endpoint_costs` / 164 | Myers 位向量计算每个 read 终点的 infix 编辑代价；目标前缀免费，故与全局编辑距离递推不同 |
| `flank_candidates` / 199 | 枚举局部极小终点，用 edlib 反向 SHW 找起点、NW 回溯；近似位置不因已有一个精确命中就被忽略 |
| `parse_records` / 305 | 解析全输入 BAM 的 SAM 文本；保存低 MAPQ、secondary、supplementary 等状态，核验 CIGAR/SEQ/QUAL 长度及目标 contig 长度 |
| `consistent_sequences` / 356 | 按原始 query 坐标比较同分子各记录，检测长度或重叠序列冲突；不拼接硬剪切片段 |
| `localization` / 374 | 由符合 MAPQ 与标志要求的同分子记录，核验侧翼 query 位置与参考位置的一致性；允许 supplementary 提供定位 |
| `competing_mapping` / 400 | 发现同一 query 侧翼在其他基因组位置的竞争，即使是低 MAPQ secondary 也保留歧义；read 其他不相关区段不构成否决 |
| `choose_pair` / 435 | 左右锚点必须正确排序；切点并列、不同位置的接近竞争对均返回未解 |
| `extract_parsed` / 461 | 各合格模板长度须给同一对切点；单条 SEQ 连续切出 `[qL,qR)`；重复内部有 N、空区间或冲突不输出完整分子 |
| `extract_bams` / 597 | 先将所有已提供 BAM 的记录按样本合并，再判断分子；命令为 `samtools view -h BAM`，不另加区域或 MAPQ/FLAG 预过滤 |

关键默认值：

| 参数 | 默认值 | 分母/语义 |
|---|---|---|
| 侧翼模板长度 | 150、250、350、500 bp | `--anchor-length` 为最大模板长度，保留不超过上限的所有档位 |
| HiFi 侧翼 identity | 0.90 | 匹配碱基数 / 比对列数 |
| ONT、ONTUL 侧翼 identity | 0.80 | 同上，不是平台准确率声明 |
| observed fraction | 0.90 | 同时消耗侧翼和 query 的位置数 / 侧翼长度 |
| localized fraction | 0.80 | 得到合格定位支持的模板位置数 / 模板长度 |
| 定位容差 | 50 bp | 原 query 坐标之间的允许差值 |
| 单侧候选上限 | 100 | 超出返回未解，不截断后宣称唯一 |
| 最少独立侧翼证据 | 每侧2个不重叠31-mer | 每个31-mer至少28个模板位置出现在锚点映射中 |
| 定位 MAPQ | 默认20，255排除 | 仅限制定位证据，并非提前删除所有低 MAPQ 输入记录 |

该实现保守而有代价：完整 BAM 的 SAM 文本由 `capture_output=True` 一次读入内存，并继续保留记录对象。应使用已经按位点整理的 BAM；它不是为直接吞入全基因组大 BAM 而设计的流式实现。本次未运行真实数据内存基准，不据此给出具体容量或性能结论。

### 4. vntrcore.py：候选发现与单元注释

#### 4.1 候选和分配是两个阶段

`deduplicate_reads`（107）按样本、namespace、名称归组，同分子出现不同完整序列时整组进入冲突审计；相同序列的副本只保留一个代表，但合并来源和平台标签。代表选择优先 HiFi，再比较锚点和质量总分。

`edit_distance`（138）计算完整 Levenshtein 编辑距离；没有“长度相同就认为相同”的捷径。`align_to_seed`（179）逐步增加回溯带宽，并与这个精确距离核对；默认最多600万个 traceback 单元，达到限制显式返回 `AlignmentUnresolved`。

`_candidate_groups`（289）先收集达到 `minimum` 的精确序列种子，其余读段按平台顺序加入最相近且非并列的种子：HiFi 相似度阈值0.96，其他平台0.82。这是候选启发式分组门槛，不是后续要设计的显示聚合阈值。`split_stable_variants`（310）再检查比对列中的稳定替换/插入/缺失模式；两侧均达到独立支持后递归拆分，降低多数共识抹掉稳定小变异的风险。

`make_consensus`（247）从真实读段选择 medoid 骨架；有 HiFi 时骨架来自 HiFi，但中心代价和逐列投票使用全部去重分子。插入只有严格多数才加入，普通列未过半则保留骨架碱基，删除多数时去掉相应碱基。共识可能是原来没有逐条观察到的组合，所以结果是候选序列，不是直接观察分子。

`generate_candidates`（356）逐样本生成共识，合并相同共识并保留种子来源。`assign_reads`（392）将每条 read 与同样本全部候选重新比较：

- 先用0.96/0.82相似度门槛判断 poor_fit。
- 比对代价为替换、普通插入、普通删除与同聚物 indel 的和；HiFi 同聚物权重1.0，其他平台0.5。
- 最佳与次佳之差小于 `min_margin`（默认0.5）时可能 ambiguous；任何候选因资源限制未能完成评分，会让该 read 保持 unresolved。
- 输出 assigned、ambiguous、unexplained、unresolved 及每候选分数，而非只存一次最佳命中。

`call_candidates`（434）还给 unexplained 残余一次新候选发现机会。平台支持数允许重叠，独立分子数取 `read_ids`。它对每个平台单独构建共识，标记跨平台不一致；在支持候选不超过12个时，枚举大小1到倍性的候选组合，累加描述性拟合代价。状态明确区分没有支持、一个观察序列但不能判纯合、未分相支持集合、倍性冲突，不能把 `genotypes.tsv` 文件名理解为完成临床分型。

#### 4.2 模板可以近似，输出序列仍是精确切片

`discover_motifs`（532）结合局部精确串联片段、滑窗周期自相关和目录模板提出候选；周期搜索2–200 bp，自相关门槛0.70，普通发现字典数量约束为至少32或目录大小。长序列按600 bp 窗口、300 bp 步长处理。这是发现策略，不能保证列出全部真实单元。

`motif_dictionary`（575）先由参考与候选提名，再要求原始分子中的串联或双侧邻接上下文验证。`substring_read_support` 只计子串出现，`read_support` 计有重复上下文的独立分子，不能混用。目录模板即使无支持仍可能 `usable=True`，同时 `validated=False`。此函数的模板验证在全位点分子集合上计数，可以跨样本；与候选生成的逐样本门槛不同。若一个可用模板能由更短模板的形式精确拼成，会被标为 compound。

`canonical_motif`（519）将循环旋转及反向互补归到一个家族名；`motif_orientation`（523）记录给定模板相对家族的方向与相位。家族归一不等于可以任意旋转输入序列后重新计数。

`decompose`（631）在原序列的碱基坐标上做切片 DAG 动态规划：

| 转移 | 当前代价与限制 |
|---|---|
| UNKNOWN 单碱基 | 通常2.0；某些符合模板端部特征的单碱基为0.8，但仍为 unknown |
| 完整模板匹配 | 0.2 + 1.5×编辑距离；单碱基模板不允许编辑，长度≤6允许1次，较长模板约10%但最多3次 |
| 端部部分模板 | 0.8，只用于符合模板前/后缀的端部片段，不计完整单元 |
| 右端专用模板 | 保留给定方向/相位，新增模板只能用于结束于序列右边界的切片 |

每个位置保留最佳和一个次佳方案；分数相同则标注歧义，但不枚举全部并列方案。转换为 Token 后保存 exact_bases、坐标、相对模板的编辑 CIGAR，并断言拼接精确还原输入。`full_unit_fraction` 是被完整模板状态覆盖的输入碱基比例，编辑匹配也计入；它不是“与模板相同碱基比例”或准确率。

#### 4.3 基础点图为何不一定简洁

`build_graph`（704）的节点键包含完整 token 路径的哈希及当前位置。因此两条不同路径即使共享一段相同序列，也不会随意共享该段节点；只有完整 token 路径相同才复用内部节点。代价是图可能像若干并列链，但收益是避免由共享局部节点产生未观察的重组路径。SVG 的 `compress_tokens` 单独做连续精确折叠。

`graph_evidence`（730）把候选已分配分子集合映射到点图节点和边。因此 RC/read_support 不是“每条分子都与该节点完全相同”的保证，而是这些分子被分配给经过该节点的候选。两侧25 bp来自参考上下文，不是给每条 read 追加观测证据。

`dbg_evidence`（754）是另一项局部辅助检查：默认 k=min(31,max(1,候选长度//2))；在该候选已分配 reads 中收集实际相邻 k-mer 连边，统计候选路径上各边出现位置是否有至少一个分子支持。`graph_support<0.90` 触发标记。这个0.90既不是单元聚合相似度，也不是整条候选的联合支持比例。

### 5. vntrmotifs.py：共用的重复语义

| 函数 / 行号 | 作用与边界 |
|---|---|
| `primitive_unit` / 12 | 求最短精确字符串周期，不自动确立生物学单元 |
| `forms` / 21、`canonical` / 26、`definition` / 30 | 枚举旋转/反向互补，给家族名、primitive、repeat_factor及基础/复合层级 |
| `exact_runs` / 38 | 枚举局部精确串联重复并按跨度贪心选择不重叠命中；不是全部可能分解 |
| `exact_occurrences` / 63 | 收集模板各形式在序列中的所有精确出现坐标 |
| `context_support` / 74 | 要求相邻同家族或左右两边都有直接相邻单元；单次裸子串命中不够 |
| `TwoPaths` / 109 | 每终点保留两个不同前缀状态，回溯最佳/次佳；代价不是概率 |
| `base_metrics` / 150 | 分别汇总完整精确、完整编辑、部分、未知碱基及完整单元数 |
| `compound_annotations` / 170 | 目录复合块必须对齐完整基础 token；不同模式可重叠，各类计数不能相加 |
| `observed_compounds` / 198 | 从实际有序 token 中找含至少两种单元的重复模式，最多12个 token 为一个模式，不跨 unknown/partial |

### 6. 本轮确认的问题

#### C1：备选分解的 token 级歧义标记与全局结果不一致

- 定位：基线 `vntrcore.py:647`、`:654`、`:665`；加注释版对应 `:671`、`:678`、`:690`。
- 原因：内部 `convert()` 总与 `alternate_set` 比较。转换备选方案本身时，每个 step 都属于自己的集合，因此 token 级 `annotation_ambiguity` 总为 false。
- 无文件写入复现：`primary, _, audit = decompose('ACGT', ['ACAT','ACCT'], return_details=True)`。
- 结果：全局 `annotation_ambiguity=True`、score_gap=0；主方案 ACAT 的 token 标记为 true，备选 ACCT 的 token 标记为 false。
- 影响：`annotation_alternatives.json` 的备选 token 字段可能让消费方误读为该 token 不存在切分歧义。全局标记和主方案仍保留歧义；不能据此声称整个程序已经把该例判成确定基因型。D1S80 专用 `vntrdbg_analysis.py` 另行按主/备选互比重算标记，不应将这个问题无条件外推到该输出。
- 建议：转换主方案时与备选集合比，转换备选时与主集合比；补充同时检查两套 token 标记的针对性测试。
- 本轮仅记录，不修行为。

#### C2：位点 ID 的路径分量校验未排除 `.` 和 `..`

- 定位：基线 `vntrgraph.py:97`、`:307–309`；加注释版 `:99`、`:316–318`。
- 原因：`[A-Za-z0-9_.-]+` 接受 `..`，随后 `args.out/locus.name` 会解析到输出目录的父目录。
- 无副作用复现：向 `read_catalog` 提供一个仅实现 `read_text()` 的对象，内容为 `chr1\t100\t112\tID=..;MOTIFS=ACGT\n`；它会返回 `name='..'`。计算 `(Path('/tmp/new-output')/name)` 的规范路径为 `/tmp`，不需要真正创建目录或运行提取。
- 影响与边界：这是本地目录配置的校验缺口，可能使误填 ID 的输出位置不符合预期；现有非空输出目录检查会拦住部分情况。没有证据表明该程序暴露了接收远程不可信 BED 的服务，不能据此宣称远程入侵或已发生数据覆盖。
- 建议：拒绝 `.`/`..` 等特殊分量，并确认计算出的位点目录严格位于指定输出根下；加无写入参数测试。
- 本轮仅记录，不修行为。

### 7. 检查范围与注释验证

已通读上述四个模块，阅读相关核心回归测试及 Phase1、revision、anchor 测试。已有测试覆盖编辑距离独立DP核对、回溯重建、跨平台去重、硬剪切不拼接、参考等价、歧义分配、稳定小变异拆分、基础/复合计数分层、图重建等关键契约。C1 的现有用例只检查全局歧义和两个家族不同，没有检查备选 token 字段；C2 未见特殊路径分量用例。

四个文件共新增72条独立中文 `#` 注释：vntrcore 28、vntranchors 19、vntrmotifs 9、vntrgraph 16。AST（忽略行号）完全相同，原有行没有删除或修改，所有新增行都是独立注释，保持LF换行。本轮未修改参数默认值、字符串、现有docstring或算法。完整测试结果见第11节，本节最小复现不替代全量验证。

针对真实数据的提取成功率、候选准确率、内存消耗和疾病事件保留能力，本次均没有新增实测结论。


## 9. dBG、连续联合证据与终端模板


本节基于服务器提交 `5be507688fe7d1f8da4c1786fd85d6df8dfccb0e`。下列行号指本次 `annotated/scripts/` 中仅添加中文注释后的文件。检查与注释不实现此前讨论的 VNTR 有损聚合，也不改变计算结果。

### 1. 这组模块在整个程序中的位置

```text
点图提取结果 + 原始位点 BAM + 匹配参考/目录
    │
    └─ vntrdbg_inputs.prepare：冻结 REF / OBS / CAND 输入与来源
          │
          └─ vntrdbg.build_all：每个 k 生成全证据图、支持图、路径与指标
                │
                ├─ vntrdbg_report：读取既有验收状态，渲染预览并打包
                │
                └─ vntrdbg_analysis：对冻结 D1S80 输入增加两类分析
                      ├─ 分支块、相邻块联合区间、整条候选的分子证据
                      └─ 所有保存路径的 16 bp 单元注释
                            └─ 可选 vntrterminal：独立验证右端专用模板
```

这里有三种不同对象：REF 为参考目录区间；OBS 为真实分子完整切片或局部片段；CAND 为点图流程给出的候选序列。三者均可决定图形拓扑，但只有 OBS 给观测支持投票。图中可走通的组合不自动成为实际观测或推断成功的等位序列。

### 2. `vntrdbg_inputs.py`：准备冻结输入

#### 主要函数

| 函数（注释后起始行） | 作用 |
|---|---|
| `rows`，20 | 读取 TSV 为字典列表。 |
| `intervals`，26 | 将参考区间投影到单条 BAM 记录上的 query 区间。 |
| `prepare`，59 | 核对来源、读取序列、区分完整与局部观测，写出统一的 dBG 输入。 |

#### 输入校验和读取顺序

1. `models` 中该位点必须唯一且明确指定 `compact_dbg`；目录中也必须只找到一个位点，并且它不是 STR 结构模式。
2. 与点图 `parameters.tsv` 对照目录 SHA256、参考 `.fai` SHA256、参考大小与修改时间；检查当前 BAM 的路径集合是否与点图阶段一致。
3. 读取原锚定决策 `anchor_decisions.jsonl` 和完整分子表 `molecules.tsv`。分子键由样本、来源 namespace、read 名构成。
4. 对每个 BAM 运行 `samtools quickcheck` 和 `samtools view -h`，再由锚定模块解析记录。该版本将单个 BAM 的 SAM 文本一次性读入内存，并对 BAM 内容做 SHA256。

#### 完整与局部 OBS 的区别

完整 OBS 沿用 `molecules.tsv` 的 `qL/qR` 和精确序列，并要求当前某条匹配来源记录的 `SEQ[qL:qR]` 完全一致。其 `coverage=two_flank_anchored`。

未得到完整 OBS 时，程序可从比对投影提取局部片段，但排除已有明确分子冲突的情况。局部提取要求目标染色体、非 unmapped/secondary/QC-fail/duplicate，且 `MAPQ>=20`、`MAPQ!=255`；补充记录可以分别提供片段，但不会被拼接。`N` 操作断开区间，内部插入的 query 碱基保留；非 A/C/G/T 碱基再次断开精确片段。局部结果标为 `alignment_local_only`，边界标为 `CIGAR_projected_uncertain`。

同一分子的相同原始 query 区间与相同碱基被合并为一个片段，来源和平台取集合；同一分子的不同片段仍保留共同分子键，后续不会按多个独立分子计票。

#### 输出

- `input_paths.json`：REF、OBS、CAND 的精确序列、哈希、范围、分子身份及来源。
- `motif_annotations.json`：候选及参考的单元注释，用于路径图展示；拼接注释必须重建原序列。
- `input_audit.tsv`：每个 BAM 的解析统计与哈希。
- `molecule_audit.tsv`：每个分子被保留、排除或无法形成局部序列的原因。
- `classification.tsv`、`input_summary.json`：模型选择、数量统计和点图源文件哈希。

#### 边界

“当前 BAM 清单相同”不等于按内容哈希证明它与最早锚定时的全部 BAM 完全相同。程序验证完整切片并记录当前 BAM 哈希；报告应准确描述这些检查，不能扩大为整个输入历史的内容级认证。参考身份也采用大小、mtime、FAI 哈希及读取到的区间，而不是全 FASTA SHA256。

### 3. `vntrdbg.py`：精确有向 compact de Bruijn 图

#### 主要函数

| 函数（起始行） | 作用 |
|---|---|
| `words`，41 | 按序列原方向取有序 k-mer。 |
| `raw_graph`，46 | 构建原始节点、实际相邻边、有序路径及 OBS 分子集合。 |
| `supported_records`，77 | 根据同一样本内独立支持过滤，并把断开的路径拆成片段。 |
| `compact_graph`，122 | 压缩无分支链为 unitig，重建并验证每条明确路径。 |
| `write_gfa`，210 | 导出 S/L/P 和支持、来源标签。 |
| `color` / `topology`，223 / 230 | 设置显示颜色并调用 Graphviz。 |
| `paths_svg`，263 | 按实际保存路径展开节点遍历顺序与可用单元注释。 |
| `write_outputs`，316 | 写出每种视图的完整产物。 |
| `metrics`，337 | 计算节点边数量、候选边保留率、重复 k-mer 比例等。 |
| `build_all`，369 | 为全部 k 分别运行完整与过滤图流程。 |
| `main`，397 | 解析参数、拒绝覆盖非空输出目录、准备输入并构图。 |

#### 图语义

原始节点是精确、定向的 k-mer；只有输入路径中真实相邻的两个 k-mer 才建立边，不能只因为后缀/前缀相同就补边。序列反向互补不被自动合并。长为 L 的输入有 `max(0,L-k+1)` 个 k-mer 位置，重复位置保留在路径中；节点、边的支持则通过分子集合去重。

压缩仅连接出度和入度均允许的非分支链，且任何输入路径的起点、终点都会阻止越界压缩。unitig 序列由第一个 k-mer 加上后续 k-mer 的末碱基构成。相邻 unitig 的重叠为 `k-1`，所以重建路径时：

```text
重建序列 = 首节点全部序列 + 后续每个节点去掉前 k-1 bp 后的序列
```

代码要求此结果精确等于原路径序列。unitig 是图压缩单位，不是生物学重复单元；同一节点的重复访问也不能直接解释为 CN。

#### 支持图

默认 `--min-support=3`。对每个原始节点或边，分别计算各样本独立分子数，取最大值与门槛比较；不会把不同样本的低票数相加过门槛。这个规则生成的是位点合并视图，不保证每个样本都支持每条保留边。

REF 显式保留为零观测支持背景。其他路径经过低支持节点或边时拆段：例如原路径两端仍有支持，中间断开后只导出两个片段，不补一条跨越缺口的连接。片段保留 `parent_id`、`parent_start/end`，并标记 `filter_complete=False`。

#### 主要参数和输出

- 默认 k：17、21、31、33、49、65；要求 k 不重复且均不小于 2。
- 每个 `k<k>/` 同时生成全证据及 `.supported` 的 GFA、JSON、SVG、路径表和支持表。
- `paths-svg/` 保存每条明确路径的图；`summary.paths.svg` 仅包含 REF/CAND。
- `excluded_short_paths.tsv` 列出短于 k、不能进入该 k 图的路径；它们仍保存在冻结输入中。
- GFA `S` 的序列是精确 unitig；`L` 的重叠为 `k-1 M`；`P` 是有序输入路径。`ms` 为独立 OBS 分子数，`or` 为来源类型集合。

`candidate_edge_retention` 的分母是候选路径中的边出现次数，而不是候选数量。`candidate_repeated_kmer_fraction` 统计候选内属于重复 k-mer 的位置比例；它不直接衡量生物学错误率。节点数不超过 450 时用 Graphviz `dot`，否则用 `neato`；布局不会改变序列或支持。

### 4. `vntrdbg_analysis.py`：跨分支联合支持和 16 bp 单元注释

#### 输入与适用范围

此入口明确限制 `input_summary.json` 的位点为 `D1S80`，使用固定 `UNIT_BP=16`。主要输入是 dBG 的 `input_paths.json/input_summary.json` 及配套点图的 `motifs.tsv/tokens.tsv`。它核对点图文件哈希、路径数、路径 ID 唯一性及各序列哈希。默认 k 为 21、33，可指定其他不重复且 >=2 的值。

这是对冻结输入的附加分析，不重新调用候选、不改变原 dBG、不更新点图基因型。模块里的“事件”不是标准化后的生物学 SNV/indel。

#### 分支与事件定义

`branch_queries`（55）先以 REF/CAND 形成假设边，并加入一个样本中至少 3 个分子支持的 OBS 边。某候选路径经过具有多个出边或入边的位置即产生分支位置。连续位置及碱基范围重叠的 k+1-mer 签名合成一个分支块。

每条候选产生三类查询：

- `B`：单个分支块。
- `J`：相邻两个分支块以及它们之间的全部序列。
- `FULL`：整条候选序列。

`anchored_interval`（133）要求左右 k-bp 上下文在目标中各自唯一、顺序正确且不重叠；不会在重复命中中任取一个。`define_events`（142）用同一样本其他候选、相同可定位上下文内的精确序列构造竞争状态；相同序列的候选属于同一个状态。

#### 两套支持概念

`joint_support`（218）保留两个层次，不能混读：

1. `exact_support_molecules`：一个分子的某条连续 OBS 中出现过精确目标串；可能存在重复定位，所以是原始精确见证计数。
2. `support_molecules/conflict_molecules/uninformative_molecules`：结合唯一上下文和已定义状态后的可判断证据。目标状态为支持，显式竞争状态为冲突，其他差异、测序错误、重复定位或缺失覆盖为无信息。

`FULL` 必须是双侧锚定 OBS 的整条序列精确相等；局部片段即使序列相同也不能支持整条候选。一个分子的不同片段不连接，同一分子的重复命中只记一票；同一分子同时产生支持和冲突则整体为无信息。

```text
informative_molecules = support_molecules + conflict_molecules
support_fraction = support_molecules / informative_molecules
```

分母为零时输出 `NA`。例如支持 9、冲突 1、无信息 90 时，该字段为 0.9，但这不是 90% 分型准确率，也不是 90% 等位频率。候选原分配分子仅作元数据及辅助计数，不用于筛掉未分配的观测。

#### 16 bp 注释

`primitive_templates`（299）只选 usable、base、长度为 16 bp 且最小周期确为 16 bp 的模板。`primitive_units`（308）对每条 REF/OBS/CAND 调用核心固定相位动态规划分解器，保存主注释及一个次优注释；不是机械每 16 bp 切一刀，也不是穷举所有并列最优方案。

每个单元保留实际碱基、实际半开坐标、模板、方向、相位、编辑操作、歧义和来源路径。16 bp 模板允许匹配含 indel 的 15/17 bp 实际片段，`observed_unit_bp` 与编辑操作记录差异。主方案和备选方案分别执行无损拼接和连续坐标检查。

汇总中需要区分：

- `full_16bp_units`：该注释方案的完整模板单元数量。
- `partial_units` 和 `unknown_bp`：不完整或未解释部分。
- `length_equiv_cn=sequence_bp/16`：描述性长度比。
- `primitive_cn`：满足当前注释条件时才给值；片段、未知/部分或并列最优歧义均保持 `NA`。即使两个显示方案数量一样，也不消除未枚举方案的歧义。
- 终端模式开启后，所有 `primitive_cn` 一律 `NA`，另外输出 `annotated_full_units`。

#### 输出与函数索引

`template_for`（289）按方向和相位恢复实际注释模板；`terminal_evidence`（438）执行可选终端来源校验；`main`（491）组织整个流程。

产物包括 `joint_queries.tsv`、`events.tsv`、`joint_support.tsv`、`joint_molecules.tsv`、`molecule_events.tsv`、`primitive_units.tsv`、`unit_summary.tsv`、`motif_counts.tsv` 和 `manifest.json`。各表分别保存“查什么”“如何定义状态”“总体计数”“分子判断”“单条记录判断”“实际分解”“摘要”“家族计数”“来源和参数”。

### 5. `vntrterminal.py`：右端模板的独立证据

核心函数只有 `validate_terminals`（18）。它接受冻结路径、完整分子表、锚定决策、原始内部模板、真实右侧 16 bp 参考上下文以及原区间边界，返回候选模板表、逐分子证据和按样本分组的已接受模板。

模板可由 REF 或完整 OBS 末 16 bp 提名。CAND 既不提名也不投票，REF 只提名不投票。判定一个 OBS 为有效支持需同时满足：

1. 双侧锚定且原锚定状态为 `complete_support`；路径哈希正确。
2. 同一分子的完整记录在序列、query 坐标、原始方向坐标和来源记录上不冲突。
3. `molecules.tsv` 中完整 query 的原切片能精确恢复冻结 OBS。
4. 末端内侧 16 bp 距原始内部模板或其反向互补至多 1 次编辑。
5. query 中末端之外的 16 bp 与真实参考右上下文精确相等。
6. 末端模板为 16 bp 原始周期，且同一样本至少 3 个独立分子支持。

同一家族的不同循环相位不合并支持。不同平台小计可能重叠，不可直接相加。接受模板仅用于本样本完整区间的右端；局部 OBS 不启用，REF 注释可使用各样本接受模板的并集。注释阶段仍可能保留末端模板相对于原碱基的编辑，因此“模板经精确验证”不能理解为“以后所有末端注释均精确匹配”。

`terminal_evidence` 会额外核验原始分子表/决策/参数哈希、目录哈希、`boundary_shift=0`、目录声明为 confirmed、原参考区间与右上下文。这里的 confirmed 来自输入目录，不能扩展为程序新完成了实验边界确认。没有移动原边界，也没有新增全基因组唯一性结论。

### 6. `vntrdbg_report.py`：历史运行的预览和打包

`preview`（14）按 SVG `viewBox` 选择限制宽或高，调用 `rsvg-convert` 保持比例生成 PNG，并读取 PNG 宽高检查边界。

`main`（25）固定使用 `results/dbg/D1S80` 和历史验收目录。它先读取图验收、STR 回归结果并要求均为 PASS，然后生成各 k 的全证据图、支持图和一个完整候选预览，写 `RESULTS.md`、文件 SHA256 清单及 tar.gz。它是专用产物整理脚本，不是通用位点分析入口；也没有重新运行验收测试。

### 7. 检查发现和复现证据

#### 已复现：候选样本没有 OBS 时，联合证据导出崩溃（建议 P2）

- 位置：`vntrdbg_analysis.py` 注释后 553 行；baseline 529 行。
- 条件：全局有候选和 OBS，但候选属于样本 S，所有 OBS 属于样本 T。现有前置检查允许进入分析。
- 行为：`joint_support` 正确返回候选样本的 `tested_molecules=0` 和空 evidence；随后 `list(evidence[0])` 引发 `IndexError`，导致运行中断并留下部分输出目录。
- 合成复现包含一个 REF、一条 16 bp 的 S 样本候选、一条 T 样本 OBS，所有序列和输入文件哈希自洽。
- 本地证据：`reproduce_dbg_review.py`、`repro_dbg_no_sample_observation.json` 和 `repro_dbg_no_sample_observation/`。
- 建议后续修复：为空表使用显式字段定义，或者提前明确拒绝此输入；若允许零证据样本，应导出空分子表及正常的零支持摘要。
- 本次没有修改该行为；也没有证据表明当前冻结真实 D1S80 运行受到影响。

#### 参数呈现限制：三处“支持门槛”并不自动联动

1. `vntrdbg --min-support` 可调，默认 3。
2. `vntrdbg.py` 注释后 251 行的边颜色分界固定为 3。
3. `vntrdbg_analysis.branch_queries` 在主流程中采用默认 3，没有读取之前支持图门槛；`vntrdbg_report.py` 注释后 59 行的历史报告文字固定称“同样本 3 个分子”。

这不证明默认门槛为 3 的历史运行有错误，但若将门槛 1/2 的新结果交给同一报告脚本，报告文字会不准确；也不能假定联合查询已经跟随支持图门槛变化。建议报告及图例显式记录各自门槛。当前中文注释已说明颜色分界固定为 3。

#### 扩展规模时的资源限制

`vntrdbg_inputs.prepare` 的 `samtools view` 使用 `check_output` 一次保留整个 BAM 的 SAM 文本，BAM SHA256 使用 `read_bytes()` 一次读入文件。对现有按位点整理的小 BAM 未发现实际内存故障；扩展到更大的输入时应评估流式读取。此项是静态资源风险，不是已经测得的性能缺陷。

#### 已核实的正确边界

已有测试覆盖真实相邻边、路径端点保留、循环重复访问、同分子去重、REF/CAND 不投票、样本不混票、过滤不补桥、反向序列不擅自合并、随机序列重建、CIGAR N 拆分、完整候选支持限制、竞争状态、重复上下文、单元编辑重建、终端样本作用域和 CN 保留 NA 等。这里只描述测试覆盖的契约；全套测试的实际执行结果由总报告填写。

### 8. 本次注释验证

五个文件共增加 68 行独立中文 `#` 注释：`vntrdbg.py` 18 行、`vntrdbg_inputs.py` 12 行、`vntrdbg_analysis.py` 24 行、`vntrterminal.py` 9 行、`vntrdbg_report.py` 5 行。采用字节读取和写入保留原 LF，未改 docstring、字符串、条件、参数或算法。五个文件均通过 `ast.dump(ast.parse(...), include_attributes=False)` 与原文件相等检查；这证明 Python 语法结构没有变化，不替代完整运行环境下的回归验证。


## 10. read、STR 与结构展示


基线：服务器提交 `5be507688fe7d1f8da4c1786fd85d6df8dfccb0e`。检查对象为 `baseline/scripts/`，注释写入 `annotated/scripts/`。下文函数名对应当前实际实现；提及行号时均注明基线。未实现本次讨论的“保护结构的有损代表替换”。

### 一、先区分这六个文件的职责

| 文件 | 职责 | 核心输出/数据 |
|---|---|---|
| `vntrreads.py` | 收集全部 BAM 记录，逐条提取可见重复片段，生成读段级注释 | `rows → fragments → tokens`，读段 GFA/JSON/HTML、完整原始记录 |
| `vntrstr.py` | 给现有序列增加精确 STR 表达 | 精确重复块、字面碱基、主/备选分解、模型吻合状态 |
| `vntrstrgraph.py` | 将精确 STR 块构建成上下文相关的单元图 | 标准 GFA 路径、逐路径重复次数、拓扑与路径 SVG |
| `vntrstructure.py` | 种子/参考/HiFi 学习字典，给全部读段做结构分解并生成共享图 | 精确图、抽象模板图、字典、模式与审计数据 |
| `vntrsummary.py` | 将已注释片段归纳为结构骨架、次数亚型和局部模块 | 结构类型、亚型、模块、逐片段归属、参考结构报告 |
| `structure_summary_view.html` | 读取上述 JSON 并交互呈现 | 筛选、方向切换、图与读段联动、详情弹层 |

两条主要链路是：

```text
--read-view
  → vntrreads.main / run_locus
  → collect_reads → read_bam → extract_record
  → motif_templates → annotate
  → write_outputs

--structure-view
  → vntrstructure.main
  → seed_dictionary → add_reference_units
  → collect_reads
  → HiFi 分解/learn_variants，至多三轮
  → 用最终字典重新 annotate_rows（全部平台）
  → reference_report
  → emit → summarize + structural_summary + 页面/表格/GFA
```

如目录显式设置 `STR_STRUCTURE`，两个读段视图都可以额外调用 `vntrstr`，但这层注释不会重新提取边界或重新分型。`vntrstrgraph` 是由主程序按对应模式调用的独立精确 STR 图构建器。

### 二、vntrreads.py：全部读段证据，而非严格候选筛选

#### 1. 数据层次

`row` 对应一个 BAM 内的 QNAME；`fragment` 对应该 read 的一条可显示比对记录切片；`token` 是片段内的一段单元/残片/未知序列。补充与次级比对各自保留，不能把几个片段拼成未观察到的完整等位序列。

`read_bam` 调用 `samtools view -h <bam>`，没有传入区域过滤，因此会读取该 BAM 的全部记录。它核对目标染色体是否存在以及 BAM/参考长度是否一致。没有重复区序列的记录仍写入状态行，例如 `flank_only`、`sequence_unavailable`、`unmapped`。

展示行 ID 由 `[sample, platform, bam.name, QNAME]` 的摘要生成，作用是区分显示行；这不是后续结构模块的分子去重键。

#### 2. 区间提取

- `cigar_interval` 检查 CIGAR 是否能完整解析、消耗的 query 长度是否与 SEQ 一致。`M/=/X` 推进两种坐标，`I` 推进 query，`D/N` 推进参考，`S` 推进 query 但不自行认定为重复序列；`H` 只记录不可用硬剪切碱基。
- 重复区两端插入均保留。`N` 会产生 `reference_skip` 标记，后续覆盖状态变为 `split_reference_coverage`，不会误写为完整覆盖。
- `near_anchor` 首先寻找唯一精确命中，否则在 CIGAR 预期位置附近做近似匹配。这里没有做全基因组锚点唯一性验证。
- `extract_record` 使用两侧局部锚点改善切片边界，锚点匹配阈值为 0.72，且必须与 CIGAR 推定边界相差不超过 120 bp。锚点次序冲突时回退 CIGAR，并记录冲突，不直接删除记录。
- `coverage` 与 `anchor` 分开：前者描述参考投影覆盖，后者描述命中的两侧锚点。`spans_candidate_interval` 不是“严格确认两侧唯一锚定”。
- BAM 负链标志不触发对存储 SEQ 的二次反向互补；切片坐标始终对应存储序列。

#### 3. 注释算法

`motif_templates` 优先读取可用的验证种子文件和 BED 目录 motif；缺少这些模板时才从参考发现。复合单元会补充其基本周期。

`annotate` 使用精确短种子确定候选匹配起点，再用动态规划为片段选择解释。HiFi 候选编辑容忍量按单元长度的 12% 计算，其他平台为 22%，再经上下限约束。它们是匹配候选范围，不是质量概率。

当前代价：普通未知碱基每个 2；特定弱末端位置为 0.8；完整单元为 `0.5 + 1.5 × 编辑距离`；允许的精确末端残片为 0.8。相邻未知块会合并显示，但保留所有碱基。每个 token 有实际切片、模板、编辑距离、方向、相位和歧义标志。

`annotation_results` 的串行与多进程版本都保持输入顺序，使输出注释与片段一一对应。`write_outputs` 再次检查 token 拼接结果及每个 token 对应的切片。

#### 4. 输出与限制

- `records.jsonl`：原始记录状态、完整存储 SEQ/QUAL，方便核查剪切与坐标。
- `<locus>.reads.gfa`：每条片段独立成链。节点 ID 含片段及局部序号，`RC:i:1` 表示该独立节点的单条记录属性，不是共享图的跨分子支持。
- `read_nodes.tsv`：逐 token 的坐标、真实碱基和注释。
- `reads.tsv`：所有展示行，包括没有可提取重复序列的行。
- HTML 每页 20 行，静态 SVG 每页 20 行；分页不等于抽样。
- `parameters.json`、`status.json`：参数、源码摘要及数量核验。

这种视图的目标是检查所有可见证据。低 MAPQ、QC 标记和仅单侧锚定仍可显示，不能直接把显示行数称为合格完整分子数。

### 三、vntrstr.py：精确 STR 注释

`parse_structure` 只接受 `(2..6 bp ACGT motif)*` 与 ACGT 字面片段的组合。单核苷酸基本周期被排除；`validate_config` 还会先检查原始目录 motif 长度，避免长 VNTR 经基本周期约简后进入该试验性 STR 模式。

`annotate` 对配置中的精确 motif 集合做动态规划，不自动循环移位或容忍错配。连续至少两次才形成 `repeat` 块；单次单元、突变、中断、N 和不完整尾部保存成 `literal` 块。目标按字典序最小化 `[字面碱基数, 重复块数]`，保存主方案以及至多一个不同的备选方案。同分才标记切分歧义；备选集合不穷举。

重要区别：切分时使用模型里的 motif 集合，并不强制所有重复块按照 `STR_STRUCTURE` 中的顺序排列；完整语法的顺序和固定字面片段通过额外的正则全匹配计算 `model_match`。因此，`model_match=False` 仍可得到可重建的精确表达；这符合保留不符合模型的实际序列这一职责。

`reconstruct` 将 `motif × count` 和字面碱基拼接回输入。主方案和备选都做断言验证。`reverse_blocks` 同时反转区块顺序、坐标和每一段的碱基方向，避免只反转 motif 而留下中断方向错误。

例如配置 `(CAG)*`，输入 `CAG×5 + CAA + CAG×7 + CA`，表达为 `(CAG)5–CAA–(CAG)7–CA`。保留下来的 CAA 与 CA 不被换成 CAG，也不贡献 CAG 次数。

主要输出：`str_paths.json` 保留全部注释，`str_paths.tsv` 是摘要，`str-svg/` 每页最多 40 条。结构框宽度不与 bp 成比例。这里的 count 是本条序列中精确游程的次数，不是新推断的基因型或局部片段的整个位点 CN。

### 四、vntrstrgraph.py：上下文相关的 STR 单元图

`build_graph` 接受已经带精确 STR 注释的显式记录，以及左右参考侧翼：

1. 验证路径 ID 唯一且格式安全；对每条注释重新拼接并检查序列哈希。
2. 将每个区块编码为 `(块类型, motif或字面碱基)`，构建不含次数的前缀树。
3. 逆序计算包含 token、是否终止以及后继状态的签名，合并相同的后缀状态。
4. 为每条路径把重复次数展开成多次访问，增加左右参考节点，再验证拼接结果。
5. 对每个节点和边按传入分子 ID 的集合计算支持；一个分子重复经过节点十次仍只贡献一份支持。

该设计允许前后多个 CTG 区段分别保留上下文节点，也允许相同结构但次数不同的路径共享节点。次数保存于 `walk` 和 `runs`，不把一个固定次数放在共享节点上。`compact_walk` 中的访问次数 n 与自环经过 n−1 次不能混淆。

左右侧翼是附加参考上下文，不证明输入分子真正覆盖到它们。图本身也可能容许显式路径之外的其他遍历或次数；只有保存的 P 路径是该输出明确记录的序列。

`write_gfa` 写展开的标准 P 记录，并附加非标准但 ASCII 的 `rp` 游程标签；JSON 保存全部元数据。`write_topology` 使用 Graphviz 展示精确 S/L 拓扑；`write_paths_svg` 逐路径绘制实际节点 ID 和次数。两种图服务于不同问题，拓扑自环不能单独提供某条序列的 CN。

### 五、vntrstructure.py：结构字典、近似注释和两套图

#### 1. 字典如何产生

- `seed_dictionary` 要求 BED 明确提供非空 ACGT 种子。
- `normalize_dictionary` 为每个模板计算基本周期/家族信息，补充缺少的基本周期，并将能由较短单元拼成的模板标记为复合层。后续单元分解只使用基础层。
- `add_reference_units` 遍历种子的起始相位，优先选择近似命中多、精确命中多的相位，加入参考中的单元变体；又独立扫描精确周期，减少完全依赖种子的偏差。参考来源本身不增加分子支持。
- `learn_variants` 只使用 HiFi、主比对、MAPQ≥20 且无 duplicate/QC-fail 标记的片段。候选来自局部精确串联、已匹配单元的变体和两侧有单元的未知中断，并要求该分子上存在精确重复上下文。
- 候选按相位等价家族分组。实际准入条件是跨样本去重总支持数至少 `minimum`，且样本内最高比例至少 `min_fraction`；这不是“同一样本至少 minimum 条”。默认分别为 3、0.10。
- 默认 `--max-variants=256`，超限报错，不悄悄截断。最多学习三轮；每轮新增后重新注释 HiFi，最后用最终字典注释所有平台。

#### 2. 单元分解

`segment` 测试模板正向及反向互补，使用短种子限制起点，`TwoPaths` 记录主方案和一个备选。固定模板相位与相对模板方向分开记录，不能把反向模板命中自动解释为染色体倒位。

候选编辑距离上限：单碱基为 0；2–6 bp 为 1；更长单元为 HiFi 长度的 32% 向下取整、其他平台 38%，最少 1。允许宽度偏移另有 8%/16% 及最多 8 bp 的限制。

最终代价为：未知每碱基 1；单元 `0.1 + min(3, 0.6×模板长度)×编辑距离 + 2×长度差绝对值`；满足长度要求的精确末端残片 0.2。高编辑上限只是候选范围，最终是否采用由代价决定。它不是用户讨论的90%相似度聚合，也不意味着输出序列被改写。

主/备选 token 都保存 `exact_bases`、实际区间、模板、方向、相位、编辑距离。未知可合并为显示块，但实际碱基不删除。

#### 3. 图、模式和支持计数

`summarize` 生成三个不同层次：

- 精确节点 `E...` 按实际切片碱基共享，逐路径可重建真实序列。
- 结构节点按 motif ID 与方向共享；unknown 和 partial 使用实际碱基摘要区分，不能把所有未知串当成同一序列。
- 家族显示节点把 unknown/partial 汇总成 U/P 桶，只用于概览。

`support_key` 使用 `[sample, QNAME]`；同一分子跨 BAM、平台标签和补充比对不重复增加支持。`occurrences` 则累计实际经过次数，可能包含多条记录的重复观察。精细模式按覆盖状态及有序结构节点序列分组，既不是按真实全序列精确分组，也不是确认的等位型。

`run_blocks` 写连续次数；`compound_blocks` 提供独立的复合结构层。基本单元、复合块和目录块可能重叠，不能把三层计数相加。

`write_gfa(..., abstract=False)` 写实际序列；抽象模式写 `S=*` 并通过 `TS` 保存模板。抽象图不能用于重建观测碱基。

#### 4. 主要输出

`emit` 生成精确/结构 GFA、字典及候选表、结构类型、局部模块、逐片段归属、节点/边/模式表、备选分解、参考报告、完整 JSON 和交互 HTML。静态 read SVG 每页 40 行；所有读段保留。`overview.svg` 是完整结构类型摘要，`family-transitions.svg` 才是家族邻接概览。

### 六、vntrsummary.py：完整类型与局部模块是两种证据

`blocks` 仅压缩相邻相同完整单元；`structure_text(..., skeleton=True)` 用 ×n 表示骨架，忽略次数但不消除中断的次序。

`structural_summary` 的第一遍给片段生成骨架 ID 与次数亚型 ID。以下任一条件使片段不参与完整类型：覆盖不全、含 unknown/partial、反向候选、注释歧义、空序列。同一分子的合格主比对若出现多个骨架/次数组合，就作为冲突分子，归为未分配。

完整类型的强证据要求：HiFi、合格主比对、覆盖候选区间、双侧锚点且两侧局部锚点匹配分数至少0.9。某类型要标记 `supported`，必须至少一个样本同时达到默认3分子及该样本合格HiFi分母的10%。该0.9是锚点匹配分数，不是整条等位序列正确率。

`local_candidates` 从同一片段中寻找两类局部结构：

- `tandem`：基本 token 周期2–12，至少两轮，选当前起点最长跨度，并将循环相位规范化用于模块ID；保留实例的实际相位。
- `shared_block`：2或3个相邻完整游程，单个分子可以只观察到一次。

局部模块不跨 U/P 或反向候选，也排除所选 token 的切分歧义。模块状态使用全批次合格HiFi分子总数达到默认3，而非完整类型的样本内比例规则。`exact_support` 和 `approximate_support` 可以重叠：同一分子可能在不同位置同时含精确和近似实例。

代表来自实际观测片段：亚型内部优先强边界、再优先编辑代价较小者；类型代表来自支持最多的亚型。各段次数的最小/最大范围只是描述，不能随意组合成新路径。

`reference_report` 给出核心序列及上下游窗口的精确周期扫描、坐标、来源摘要和边界复核提示。扫描未发现跨界重复不等于确认边界。

### 七、structure_summary_view.html：只呈现，不重新分析

模板包含 `__TITLE__`、`__GRAPH__`、`__DATA__` 三个占位符。Python 将完整数据嵌入 JSON 脚本，转义 `< > &` 后由页面解析。前端建立 motif、节点、read、fragment 和 assignment 索引。

- `matchingAssignments` 统一应用样本、平台、证据、结构类型、搜索和模块筛选，让结构、模块、精细模式及 read 同步变化。
- `stateTag` 保留全量分析的类别状态；`support` 按当前筛选内的分子集合计算数量。筛选后支持变少不会重新改变原全量状态。
- `renderReference`、`renderDictionary` 展示参考、种子来源、实际序列和不同支持口径。
- `renderGraphs`、`renderClasses` 每页8类；类型使用实际代表片段，次数亚型单独展开。
- `renderModules` 每页10模块，实例表默认前30个，完整实例仍在导出文件；局部支持不等于完整路径支持。
- `renderPatterns` 每页30种细模式；`renderReads` 每页20行，保留未分配片段和无片段状态行。
- 基因负链显示会反转区块次序并反向互补实际序列；存储坐标仍是原BAM坐标。部分局部模块表达明确固定为参考方向。
- 统一bp刻度允许比较片段长度；独立铺满模式会让每条片段占满宽度，不能据屏幕宽度比较长度。
- 点击类别/模块/字典联动筛选；点击单元展示真实碱基、坐标、全量节点支持和经过次数。

### 八、检查结果与后续关注点

#### 已确认的模块差异与限制

1. **分子命名空间没有统一到旧结构链路。** 基线 `vntrstructure.py:29`、`vntrsummary.py:106` 使用 `[sample,QNAME]`；旧read行身份虽然区分BAM/平台，但后续分子支持仍不用namespace。主点图的 `[sample,namespace,read_name]` 规则不能直接套用到这里。若同一样本的两个独立文库恰有相同QNAME，这里会保守合成一个分子，可能低估支持；结构不同还会标为冲突。已有测试主动验证跨BAM同名只计一次，说明这是现行契约。新报告应明确范围，后续多文库支持可单独设计，不在注释任务中改行为。

2. **旧说明存在参数漂移。** `STRUCTURE_VIEW.md` 写候选模板最多12，当前CLI默认256；旧文档的评分说明也不足以表达当前短单元使用 `min(3,0.6×长度)` 的实现。本次新解读按源码记录实际默认值，不修改历史说明。

3. **同名“支持门槛”对应不同分母。** 字典学习使用跨样本总支持+样本内最高比例；完整结构类型使用同一样本支持与比例；局部模块使用全批次合格分子总数。不能把三个字段都解释成同一样本三条分子，详解与新增注释已分别说明。

4. **STR精确分解的长纯重复性能风险已做小规模复现。** `vntrstr.annotate` 在每个位置枚举多个可能结束点，纯重复序列中候选数随长度呈二次增长，并逐位置排序。当前本地Python3.12对300/600/1200/2400bp纯CAG分别约0.009/0.043/0.187/0.945秒，表达均正确；原始结果见 `str_scaling_probe.json`。这不是分型错误，也不是服务器性能基准，但说明全量长STR read可视化的耗时不能从短序列测试直接外推。长序列可在独立优化任务中用等价的游程/动态规划优化验证。

5. **已复现结构分解备选 token 歧义标记不一致（P2）。** 基线 `vntrstructure.py` 的 `segment` 在约第153–175行创建 `alternate_set` 并通过同一个 `convert` 转换主、备选方案。两次转换都用备选集合比较，导致备选方案实际上与自身比较，备选 token 的 `annotation_ambiguity` 为假，即使全局明确存在同分歧义。最小调用为 `segment(('ACGT', seed_dictionary(['ACAT', 'ACCT']), 'HiFi'), return_details=True)`：全局 `annotation_ambiguity=True`、`score_gap=0.0`，主 token 的标记为真，备选 token 的标记为假。当前该输入的主模板为正向 `ACAT`，备选为同一模板的反向互补 `ATGT`；两者均能重建输入碱基。

   影响范围是结构分解的备选 token 元数据，包括写入 `annotation_alternatives.json` 的备选方案。全局歧义仍保留，当前页面也显示全局歧义，不能据此声称整条路径被判定为确定基因型。建议未来修复时让转换函数显式接收“另一方案”的签名集合，并加入主/备选对称标记测试。本轮源码仅添加注释，保留该缺陷供独立修复。

   可复现脚本为 `structure_ambiguity_probe.py`，结果为 `structure_ambiguity_probe.json`；脚本读取 baseline，断言问题仍可重现，并检查两种分解均精确重建输入。当前审阅没有穷尽所有输入，也不构成生物学准确率证明；全量测试由主审汇总。

#### 注释改动验证

5个Python文件共新增66行独立中文注释，HTML新增10行独立JS中文注释。未改已有docstring、字符串、表达式或逻辑；各文件保持LF。Python逐文件AST完全一致；按行diff均只有注释插入；HTML去除新增注释后与基线字节完全一致。结果保存于 `views_comment_validation.json`。

重点阅读过的相关测试：`test_all_reads.py`、`test_str.py`、`test_str_graph.py`、`test_structure.py`、`test_structure_summary.py`，并检查了 `test_phase1.py` 与 `test_revision.py` 的相关覆盖范围。它们涵盖真实切片保留、部分覆盖、负链、未知序列、重复支持去重、精确STR重建、路径次数、结构冲突及中断保留。


## 11. 本轮验证与使用边界

- 15 个源码文件仅新增 **216 行独立中文注释**：14 个 Python 模块及 1 个 HTML/JavaScript 模板。
- Python AST（忽略行号）完全一致；原有行的字节和顺序保留，LF 换行未改变。HTML 模板移除新增注释后与原版字节一致，提取出的 JavaScript 通过 `node --check`。
- Slurm 作业 **167017** 在 `cn01` 完成，状态 `COMPLETED 0:0`；Python 3.13.9。
- **186 项现有单元测试全部通过**。
- 合成 BAM 的主点图、全量 reads、结构图三项集成测试全部通过。
- 16 bp VNTR 注释前后回归通过：**17 个结果文件逐字节相同**；`parameters.tsv` 按设计记录不同源码哈希和输出路径，不参与字节相等断言。
- 测试前后源码 SHA256 校验通过。真实 BAM 和既有冻结科学结果未重跑。
- 初次作业 167016 已通过单元与三项集成测试，但我编写的验证包装脚本漏传回归所需 `--baseline`，因此退出 2；补齐参数后在独立 v2 快照完整通过。这不是被检程序单元测试失败。
- 本地 Windows 缺少 edlib 所需的 C++ 编译器，因此没有以本地依赖安装失败冒充程序失败，完整验证使用服务器已有依赖。
- 原版源码保存在本次 `baseline` 快照并随交付提供；源码未提交 Git，现有分支提交号不变。

验证范围限于上述测试和最小复现。现有测试通过并不抵消本报告 R1–R3 已复现的缺陷；这些触发条件不在原测试覆盖内。

完整检查条目和后续修复建议见 [代码检查报告](CODE_REVIEW_20261029.zh-CN.md)。

## 12. 源码函数索引

下列行号对应本次添加中文注释后的源码。HTML 阅读版提供带行号的源码浏览，模块解读中的原版行号另有说明。

### vntranchors.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `dependencies` | 函数 | [27](../scripts/vntranchors.py#L27) |
| `canonical` | 函数 | [36](../scripts/vntranchors.py#L36) |
| `nonoverlapping` | 函数 | [41](../scripts/vntranchors.py#L41) |
| `load_equivalences` | 函数 | [50](../scripts/vntranchors.py#L50) |
| `normalized_position` | 函数 | [77](../scripts/vntranchors.py#L77) |
| `reference_profiles` | 函数 | [89](../scripts/vntranchors.py#L89) |
| `endpoint_costs` | 函数 | [164](../scripts/vntranchors.py#L164) |
| `AnchorOverflow` | 类 | [194](../scripts/vntranchors.py#L194) |
| `flank_candidates` | 函数 | [199](../scripts/vntranchors.py#L199) |
| `Record` | 类 | [254](../scripts/vntranchors.py#L254) |
| `Record.original_base` | 函数 | [273](../scripts/vntranchors.py#L273) |
| `Record.original_interval` | 函数 | [277](../scripts/vntranchors.py#L277) |
| `Record.original_sequence` | 函数 | [281](../scripts/vntranchors.py#L281) |
| `Record.original_start` | 函数 | [285](../scripts/vntranchors.py#L285) |
| `mapped_bases` | 函数 | [290](../scripts/vntranchors.py#L290) |
| `parse_records` | 函数 | [305](../scripts/vntranchors.py#L305) |
| `consistent_sequences` | 函数 | [356](../scripts/vntranchors.py#L356) |
| `localization` | 函数 | [374](../scripts/vntranchors.py#L374) |
| `competing_mapping` | 函数 | [400](../scripts/vntranchors.py#L400) |
| `choose_pair` | 函数 | [435](../scripts/vntranchors.py#L435) |
| `extract_records` | 函数 | [454](../scripts/vntranchors.py#L454) |
| `extract_parsed` | 函数 | [461](../scripts/vntranchors.py#L461) |
| `extract_bam` | 函数 | [586](../scripts/vntranchors.py#L586) |
| `extract_bams` | 函数 | [597](../scripts/vntranchors.py#L597) |

### vntrcore.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `revcomp` | 函数 | [23](../scripts/vntrcore.py#L23) |
| `compact` | 函数 | [27](../scripts/vntrcore.py#L27) |
| `Read` | 类 | [33](../scripts/vntrcore.py#L33) |
| `Read.molecule_key` | 函数 | [53](../scripts/vntrcore.py#L53) |
| `Read.observed_platforms` | 函数 | [57](../scripts/vntrcore.py#L57) |
| `Token` | 类 | [63](../scripts/vntrcore.py#L63) |
| `Token.__iter__` | 函数 | [76](../scripts/vntrcore.py#L76) |
| `Allele` | 类 | [82](../scripts/vntrcore.py#L82) |
| `deduplicate_reads` | 函数 | [107](../scripts/vntrcore.py#L107) |
| `edit_distance` | 函数 | [138](../scripts/vntrcore.py#L138) |
| `similarity` | 函数 | [170](../scripts/vntrcore.py#L170) |
| `AlignmentUnresolved` | 类 | [174](../scripts/vntrcore.py#L174) |
| `align_to_seed` | 函数 | [179](../scripts/vntrcore.py#L179) |
| `alignment_counts` | 函数 | [227](../scripts/vntrcore.py#L227) |
| `make_consensus` | 函数 | [247](../scripts/vntrcore.py#L247) |
| `_candidate_groups` | 函数 | [289](../scripts/vntrcore.py#L289) |
| `split_stable_variants` | 函数 | [310](../scripts/vntrcore.py#L310) |
| `generate_candidates` | 函数 | [356](../scripts/vntrcore.py#L356) |
| `assign_reads` | 函数 | [392](../scripts/vntrcore.py#L392) |
| `call_candidates` | 函数 | [434](../scripts/vntrcore.py#L434) |
| `cluster_reads` | 函数 | [514](../scripts/vntrcore.py#L514) |
| `canonical_motif` | 函数 | [519](../scripts/vntrcore.py#L519) |
| `motif_orientation` | 函数 | [523](../scripts/vntrcore.py#L523) |
| `discover_motifs` | 函数 | [532](../scripts/vntrcore.py#L532) |
| `motif_dictionary` | 函数 | [575](../scripts/vntrcore.py#L575) |
| `cigar_edit` | 函数 | [619](../scripts/vntrcore.py#L619) |
| `decompose` | 函数 | [631](../scripts/vntrcore.py#L631) |
| `decompose.convert` | 函数 | [672](../scripts/vntrcore.py#L672) |
| `token_counts` | 函数 | [696](../scripts/vntrcore.py#L696) |
| `build_graph` | 函数 | [704](../scripts/vntrcore.py#L704) |
| `graph_evidence` | 函数 | [730](../scripts/vntrcore.py#L730) |
| `dbg_evidence` | 函数 | [754](../scripts/vntrcore.py#L754) |

### vntrdbg.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `digest` | 函数 | [20](../scripts/vntrdbg.py#L20) |
| `save_json` | 函数 | [24](../scripts/vntrdbg.py#L24) |
| `tsv` | 函数 | [28](../scripts/vntrdbg.py#L28) |
| `words` | 函数 | [41](../scripts/vntrdbg.py#L41) |
| `raw_graph` | 函数 | [46](../scripts/vntrdbg.py#L46) |
| `supported_records` | 函数 | [77](../scripts/vntrdbg.py#L77) |
| `supported_records.enough` | 函数 | [83](../scripts/vntrdbg.py#L83) |
| `compact_graph` | 函数 | [122](../scripts/vntrdbg.py#L122) |
| `compact_graph.evidence` | 函数 | [196](../scripts/vntrdbg.py#L196) |
| `write_gfa` | 函数 | [210](../scripts/vntrdbg.py#L210) |
| `color` | 函数 | [223](../scripts/vntrdbg.py#L223) |
| `topology` | 函数 | [230](../scripts/vntrdbg.py#L230) |
| `paths_svg` | 函数 | [263](../scripts/vntrdbg.py#L263) |
| `write_outputs` | 函数 | [316](../scripts/vntrdbg.py#L316) |
| `metrics` | 函数 | [337](../scripts/vntrdbg.py#L337) |
| `build_all` | 函数 | [369](../scripts/vntrdbg.py#L369) |
| `main` | 函数 | [397](../scripts/vntrdbg.py#L397) |

### vntrdbg_analysis.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `read_tsv` | 函数 | [25](../scripts/vntrdbg_analysis.py#L25) |
| `write_tsv` | 函数 | [31](../scripts/vntrdbg_analysis.py#L31) |
| `occurrences` | 函数 | [42](../scripts/vntrdbg_analysis.py#L42) |
| `branch_queries` | 函数 | [55](../scripts/vntrdbg_analysis.py#L55) |
| `anchored_interval` | 函数 | [133](../scripts/vntrdbg_analysis.py#L133) |
| `define_events` | 函数 | [142](../scripts/vntrdbg_analysis.py#L142) |
| `assess_observation` | 函数 | [185](../scripts/vntrdbg_analysis.py#L185) |
| `joint_support` | 函数 | [218](../scripts/vntrdbg_analysis.py#L218) |
| `template_for` | 函数 | [289](../scripts/vntrdbg_analysis.py#L289) |
| `primitive_templates` | 函数 | [299](../scripts/vntrdbg_analysis.py#L299) |
| `primitive_units` | 函数 | [308](../scripts/vntrdbg_analysis.py#L308) |
| `terminal_evidence` | 函数 | [438](../scripts/vntrdbg_analysis.py#L438) |
| `main` | 函数 | [491](../scripts/vntrdbg_analysis.py#L491) |

### vntrdbg_inputs.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `rows` | 函数 | [20](../scripts/vntrdbg_inputs.py#L20) |
| `intervals` | 函数 | [26](../scripts/vntrdbg_inputs.py#L26) |
| `prepare` | 函数 | [59](../scripts/vntrdbg_inputs.py#L59) |

### vntrdbg_report.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `preview` | 函数 | [14](../scripts/vntrdbg_report.py#L14) |
| `main` | 函数 | [25](../scripts/vntrdbg_report.py#L25) |

### vntrgraph.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `Locus` | 类 | [33](../scripts/vntrgraph.py#L33) |
| `IndexedFastaReference` | 类 | [44](../scripts/vntrgraph.py#L44) |
| `IndexedFastaReference.__init__` | 函数 | [47](../scripts/vntrgraph.py#L47) |
| `IndexedFastaReference.__enter__` | 函数 | [58](../scripts/vntrgraph.py#L58) |
| `IndexedFastaReference.__exit__` | 函数 | [61](../scripts/vntrgraph.py#L61) |
| `IndexedFastaReference.contig_length` | 函数 | [64](../scripts/vntrgraph.py#L64) |
| `IndexedFastaReference.fetch` | 函数 | [67](../scripts/vntrgraph.py#L67) |
| `read_catalog` | 函数 | [87](../scripts/vntrgraph.py#L87) |
| `discover_bams` | 函数 | [122](../scripts/vntrgraph.py#L122) |
| `reference_to_query` | 函数 | [139](../scripts/vntrgraph.py#L139) |
| `flank_hit` | 函数 | [158](../scripts/vntrgraph.py#L158) |
| `extract_sam_records` | 函数 | [207](../scripts/vntrgraph.py#L207) |
| `write_tsv` | 函数 | [218](../scripts/vntrgraph.py#L218) |
| `source_hashes` | 函数 | [226](../scripts/vntrgraph.py#L226) |
| `write_gfa` | 函数 | [231](../scripts/vntrgraph.py#L231) |
| `compress_tokens` | 函数 | [249](../scripts/vntrgraph.py#L249) |
| `write_svg` | 函数 | [262](../scripts/vntrgraph.py#L262) |
| `process_locus` | 函数 | [316](../scripts/vntrgraph.py#L316) |
| `main` | 函数 | [447](../scripts/vntrgraph.py#L447) |

### vntrmotifs.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `primitive_unit` | 函数 | [12](../scripts/vntrmotifs.py#L12) |
| `forms` | 函数 | [21](../scripts/vntrmotifs.py#L21) |
| `canonical` | 函数 | [26](../scripts/vntrmotifs.py#L26) |
| `definition` | 函数 | [30](../scripts/vntrmotifs.py#L30) |
| `exact_runs` | 函数 | [38](../scripts/vntrmotifs.py#L38) |
| `exact_occurrences` | 函数 | [63](../scripts/vntrmotifs.py#L63) |
| `context_support` | 函数 | [74](../scripts/vntrmotifs.py#L74) |
| `PathState` | 类 | [102](../scripts/vntrmotifs.py#L102) |
| `TwoPaths` | 类 | [109](../scripts/vntrmotifs.py#L109) |
| `TwoPaths.__init__` | 函数 | [111](../scripts/vntrmotifs.py#L111) |
| `TwoPaths.offer` | 函数 | [117](../scripts/vntrmotifs.py#L117) |
| `TwoPaths.trace` | 函数 | [128](../scripts/vntrmotifs.py#L128) |
| `TwoPaths.audit` | 函数 | [139](../scripts/vntrmotifs.py#L139) |
| `base_metrics` | 函数 | [150](../scripts/vntrmotifs.py#L150) |
| `compound_annotations` | 函数 | [170](../scripts/vntrmotifs.py#L170) |
| `observed_compounds` | 函数 | [198](../scripts/vntrmotifs.py#L198) |

### vntrreads.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `cigar_interval` | 函数 | [37](../scripts/vntrreads.py#L37) |
| `near_anchor` | 函数 | [80](../scripts/vntrreads.py#L80) |
| `extract_record` | 函数 | [95](../scripts/vntrreads.py#L95) |
| `read_bam` | 函数 | [145](../scripts/vntrreads.py#L145) |
| `collect_reads` | 函数 | [202](../scripts/vntrreads.py#L202) |
| `annotation_results` | 函数 | [213](../scripts/vntrreads.py#L213) |
| `motif_templates` | 函数 | [223](../scripts/vntrreads.py#L223) |
| `annotate` | 函数 | [248](../scripts/vntrreads.py#L248) |
| `annotate.update` | 函数 | [276](../scripts/vntrreads.py#L276) |
| `write_table` | 函数 | [325](../scripts/vntrreads.py#L325) |
| `svg_read` | 函数 | [330](../scripts/vntrreads.py#L330) |
| `svg_read.text` | 函数 | [334](../scripts/vntrreads.py#L334) |
| `write_outputs` | 函数 | [409](../scripts/vntrreads.py#L409) |
| `run_locus` | 函数 | [471](../scripts/vntrreads.py#L471) |
| `main` | 函数 | [509](../scripts/vntrreads.py#L509) |

### vntrstr.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `parse_structure` | 函数 | [21](../scripts/vntrstr.py#L21) |
| `validate_config` | 函数 | [43](../scripts/vntrstr.py#L43) |
| `reconstruct` | 函数 | [53](../scripts/vntrstr.py#L53) |
| `format_blocks` | 函数 | [59](../scripts/vntrstr.py#L59) |
| `reverse_blocks` | 函数 | [65](../scripts/vntrstr.py#L65) |
| `reverse_blocks.rc` | 函数 | [67](../scripts/vntrstr.py#L67) |
| `_blocks` | 函数 | [80](../scripts/vntrstr.py#L80) |
| `annotate` | 函数 | [94](../scripts/vntrstr.py#L94) |
| `path_record` | 函数 | [148](../scripts/vntrstr.py#L148) |
| `annotate_rows` | 函数 | [153](../scripts/vntrstr.py#L153) |
| `svg_blocks` | 函数 | [168](../scripts/vntrstr.py#L168) |
| `write_paths_svg` | 函数 | [194](../scripts/vntrstr.py#L194) |
| `write_outputs` | 函数 | [226](../scripts/vntrstr.py#L226) |

### vntrstrgraph.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `compact_walk` | 函数 | [21](../scripts/vntrstrgraph.py#L21) |
| `_token` | 函数 | [28](../scripts/vntrstrgraph.py#L28) |
| `build_graph` | 函数 | [34](../scripts/vntrstrgraph.py#L34) |
| `write_gfa` | 函数 | [140](../scripts/vntrstrgraph.py#L140) |
| `_lines` | 函数 | [156](../scripts/vntrstrgraph.py#L156) |
| `write_topology` | 函数 | [161](../scripts/vntrstrgraph.py#L161) |
| `write_paths_svg` | 函数 | [179](../scripts/vntrstrgraph.py#L179) |
| `write_outputs` | 函数 | [236](../scripts/vntrstrgraph.py#L236) |

### vntrstructure.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `support_key` | 函数 | [30](../scripts/vntrstructure.py#L30) |
| `seed_dictionary` | 函数 | [36](../scripts/vntrstructure.py#L36) |
| `normalize_dictionary` | 函数 | [49](../scripts/vntrstructure.py#L49) |
| `add_reference_units` | 函数 | [70](../scripts/vntrstructure.py#L70) |
| `exact_tandem_units` | 函数 | [106](../scripts/vntrstructure.py#L106) |
| `segment` | 函数 | [112](../scripts/vntrstructure.py#L112) |
| `segment.update` | 函数 | [139](../scripts/vntrstructure.py#L139) |
| `segment.convert` | 函数 | [161](../scripts/vntrstructure.py#L161) |
| `segment_details` | 函数 | [187](../scripts/vntrstructure.py#L187) |
| `annotate_rows` | 函数 | [192](../scripts/vntrstructure.py#L192) |
| `learn_variants` | 函数 | [201](../scripts/vntrstructure.py#L201) |
| `run_blocks` | 函数 | [261](../scripts/vntrstructure.py#L261) |
| `compound_blocks` | 函数 | [274](../scripts/vntrstructure.py#L274) |
| `summarize` | 函数 | [280](../scripts/vntrstructure.py#L280) |
| `summarize.node` | 函数 | [283](../scripts/vntrstructure.py#L283) |
| `summarize.add` | 函数 | [286](../scripts/vntrstructure.py#L286) |
| `summarize.edge` | 函数 | [289](../scripts/vntrstructure.py#L289) |
| `summarize.finish` | 函数 | [326](../scripts/vntrstructure.py#L326) |
| `color` | 函数 | [337](../scripts/vntrstructure.py#L337) |
| `graph_svg` | 函数 | [348](../scripts/vntrstructure.py#L348) |
| `rows_svg` | 函数 | [372](../scripts/vntrstructure.py#L372) |
| `write_gfa` | 函数 | [398](../scripts/vntrstructure.py#L398) |
| `emit` | 函数 | [413](../scripts/vntrstructure.py#L413) |
| `main` | 函数 | [504](../scripts/vntrstructure.py#L504) |

### vntrsummary.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `identity` | 函数 | [17](../scripts/vntrsummary.py#L17) |
| `blocks` | 函数 | [22](../scripts/vntrsummary.py#L22) |
| `structure_text` | 函数 | [35](../scripts/vntrsummary.py#L35) |
| `reference_report` | 函数 | [41](../scripts/vntrsummary.py#L41) |
| `is_primary` | 函数 | [65](../scripts/vntrsummary.py#L65) |
| `strong_boundary` | 函数 | [70](../scripts/vntrsummary.py#L70) |
| `local_candidates` | 函数 | [75](../scripts/vntrsummary.py#L75) |
| `structural_summary` | 函数 | [110](../scripts/vntrsummary.py#L110) |
| `structural_summary.counts` | 函数 | [166](../scripts/vntrsummary.py#L166) |
| `summary_svg` | 函数 | [213](../scripts/vntrsummary.py#L213) |

### vntrterminal.py

| 定义 | 类型 | 注释后行号 |
|---|---|---:|
| `validate_terminals` | 函数 | [18](../scripts/vntrterminal.py#L18) |
