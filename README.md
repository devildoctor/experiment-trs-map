# VNTR 点位图程序

本程序独立读取 BAM/BAI、参考 FASTA/FAI 和位点表，生成有独立分子证据的候选序列及 GFA/SVG。根据《VNTR 方案源码核查与修订建议》修订；不调用 TRGT、LongTR、vamos、danbing-tk 或 PGGB，也不声称复现它们的分型性能。

2026-10-19 主点位图锚点更新详见 [ANCHOR_UPDATE.md](ANCHOR_UPDATE.md)：多长度模板分级、整份点位 BAM 读取、分子级补充定位验证，以及可选的显式 ALT 同点位映射。展示模式的行数不能代替主图完整支持。

## 输入和边界

BAM 按 bamdata_by_locus/点位/样本/{HiFi,ONT,ONTUL}/文件名.位点.bam 整理；data_by_locus 和 bamdata 是兼容旧脚本的目录链接。要求 BAI 和匹配的参考。提供的 BAM 是区间左右各 500 bp 的完整可用输入，主点位图顺序读取全部记录，不再次添加参考区间过滤。BED 坐标为 0-based 半开区间。第 4 列支持 ID、MOTIFS，以及成对的 REPEAT_START、REPEAT_END 和 BOUNDARY=confirmed。第 2、3 列保留搜索窗口，明确重复区间必须位于窗口内。

未明确设置 BOUNDARY=confirmed 时使用候选窗口截取观测序列，标记 candidate_boundaries。当前主图目录包含 10 个位点（原有 8 个，加 ATXN2 和 D1S80）；目录中继承的 confirmed 标注并非本程序重新验证的生物学边界。首尾真实碱基不会为凑整数拷贝而修剪。

## 提取和分子证据

BAM 存储 SEQ 与同记录 CIGAR 配套，FLAG 0x10 不再翻转 SEQ。只在 SAM 文本 POS 解析时减 1。原始 query 区间另存，考虑硬剪切与方向。

全参考检查每侧 150/250/350/500 bp 内的全部 31-mer；至少两个不重叠的唯一片段才具备完整支持所需的参考区分能力。read 模板采用编辑距离匹配，内侧切点固定于参考定义，插入保留实际碱基。主/补充 CIGAR 验证匹配位置，低 MAPQ 和次级比对保留审计但不能独立提供可信定位。补充记录不拼接、不填充 N 或参考序列；只有一条记录保存完整双锚及中间实际序列时才能提供完整观测。

分子键是 [sample, namespace, read_name]。namespace 优先来自可选 sources.tsv，其次来自 BAM 的 LB；来源缺失时使用空 namespace，同一样本内相同 read ID 跨文件、跨平台保守去重。RG、LB、HP 和原始记录索引仍保留。同一分子的完整观测不一致时记录 conflicting_sequences 并排除明确支持。

sources.tsv 列名为 bam、namespace；相对 BAM 路径相对该表所在目录解析。同一原始分子的 ONT/ONTUL 数据必须使用同一 namespace，独立文库可显式区分。平台名和文件名不自动代表独立分子。

## 候选和分配

各平台均可提出候选。重复支持的同长度序列差异与簇内稳定残差模式可以拆成不同候选。共识使用实际读段 medoid 和每分子一票的修复；HiFi 可优先提供骨架。它不是质量加权 POA，不随机抽样，不纠正回参考。

距离为完整序列的精确 Levenshtein 距离，采用任意精度位向量，不退化为长度差。回溯扩大比对带宽，并与精确距离核对；超过 600 万 trace 单元返回明确 unresolved，保留候选和原因。

全部原始读段共同评分，输出替换、普通插入/缺失、同聚物 indel 和最优与次优分差。HiFi 每次编辑代价 1；ONT/ONTUL 同聚物 indel 代价 0.5，其余 1。代价在确定性的最短编辑比对上计算，未经过校准，不是概率或全局最优的质量加权比对。默认最小分差 0.5；相近候选保持 ambiguous，无法解释的读段单列。

候选允许超过 2 条。默认按倍性 2 枚举有支持的候选组合，输出描述性代价；超过 12 条候选时明确记录资源限制。超过倍性的候选保留 ploidy_conflict，不强取两个最大簇。仅观察一条序列不报告纯合。标签仅在本位点有效；保留 HP 信息不等于已实施跨位点相位。

## motif 和图

motif 发现与字典分解分开。参考和有支持的候选提出基础单元与复合块，新单元需独立分子证据。发现模块是有限窗口的探索基线，不保证发现全部生物学单元。

token 保存真实切片、起止位置、家族、相位、方向、编辑操作、full/partial/unknown 和注释歧义。基础单元与复合块分别注释，完整、部分单元和未知碱基分别计数，不以总长除单一周期。

默认 VNTR GFA 只报告显式保存的完整路径。不同路径上下文和位置使用不同节点；只有完全相同的 token 路径才共享整条内部节点链。这样不会通过共享中间 motif 拼出未观察到的重组路径。节点是互不重叠切片，连接为 0M；两端各 25 bp 是明确标注的参考上下文。显式启用 STR 的单元图另见 STR_GRAPH.md，其自环拓扑须结合明确 P 路径解释。

SVG 使用短单元 ID 和完整碱基图例，连续相同单元以 ×N 压缩。×N 表示观测序列显示次数；生物学拷贝数能否报告取决于边界和证据。局部有向 k-mer 图仅核查证据，不剪枝或自由生成等位。

始终验证 token 拼接等于共识、GFA 路径去除参考上下文后等于共识、SVG 压缩展开等于共识。支持数可追溯到分子集合；平台分项可能重叠，不得相加代替独立分子总数。

## 输出

| 文件 | 内容 |
| --- | --- |
| 位点.svg / 位点.gfa / 位点.json | 图、明确路径、节点与边的分子支持 |
| status.json / parameters.tsv | 状态、错误、边界、参数、源码与输入校验信息 |
| molecules.tsv / dedup.tsv / bam_audit.tsv | 分子、SEQ、锚点、原始记录、去重及拒绝原因 |
| candidates.tsv / audit.tsv | 全部候选 / 通过支持门槛的序列、平台冲突与残差 |
| assignments.tsv / genotypes.tsv | 所有分子的评分与分配、候选组合及不确定状态 |
| motifs.tsv / tokens.tsv | 字典与证据、逐个真实切片的注释 |
| 样本_C编号.dbg.json | 有向 k-mer 邻接、分子支持与有序路径 |

无完整跨越证据时生成说明原因的 SVG 和状态，不伪造 GFA。输入错误导致非零退出；生物学证据不足标为 unresolved。已有结果目录非空时拒绝覆盖，需使用新的运行目录。

## 验证与 Slurm

Python 3.11+ 和 samtools；主点位图另需 requirements-anchors.txt 中固定版本的 edlib 与 pyahocorasick。集群依赖部署于 .deps/vntr-anchors，由 Slurm 设置 PYTHONPATH。集群 Python 是 /home/data/user/stzhu/miniconda3/bin/python3；samtools 是 /home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools。

```bash
python3 -m unittest discover -s tests -v
python3 tests/integration_synthetic.py
```

合成 BAM 测试包含正/负链、主/补充比对去重、ONT/ONTUL 分子去重、仅由 ONTUL 跨越的扩增，并验证图的碱基重构及平台支持数。

先把 scripts、tests、catalogs、slurm 和 README.md 保存到 runs-result/ 下新的独立快照目录，并为这些文件生成 SHA256SUMS，然后提交。主图目录为 catalogs/grch38_pointgraph_10_loci_v1.bed；read 和 structure 展示使用 catalogs/grch38_structure_seeds_il1rn_per3_abca7_v1.bed；汇总视图使用 catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed。

```bash
RUN_DIR="$PWD/runs-result/新快照目录"
bash "$RUN_DIR/slurm/submit.sh" "$RUN_DIR"
```

jobs.tsv 保存作业号。ATXN2 试运行成功后，其余 9 个位点开始。日志位于快照目录的 logs/，结果位于 results/。作业执行前核验源码校验和，根目录后续改动不影响排队任务。每个位点使用 1 节点、2 CPU、16 GB 和 6 小时上限。

当前是可审计的保守基线，尚未完成真实 BAM 的性能及分型准确性验证。候选边界、平台误差模型和复杂比对歧义仍需结合真实结果核查。

## 2–6 bp STR 表达试点

ATXN2 可使用 `catalogs/grch38_str_expression_atxn2_v1.bed` 显式启用精确
STR 表达，例如 `(CTG)8–TTG–(CTG)4`。固定中断、残片与原始证据均保留；
>6 bp VNTR 和原目录默认流程不启用。规则、输出和 Slurm 入口见
[STR_EXPRESSION.md](STR_EXPRESSION.md)。

STR 点图现在输出重复单元节点、自环和标准 GFA P；静态逐路径图显示
`×n`。主图、旧图及证据文件的成套命名、Graphviz 依赖和独立验收入口见
[STR_GRAPH.md](STR_GRAPH.md)。新入口为 `slurm/vntrstrgraph.sbatch`。
