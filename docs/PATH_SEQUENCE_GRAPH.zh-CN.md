# VNTR 序列结构图：GFA S/L/P

本模块把完整观测区间序列组织为共享节点、分支、汇合与跳跃边，保留每条原始序列的有序 P 路径。独立入口 `scripts/vntrpathgraph_cli.py`；不改变已有位点图、单元划分或分型流程。

## 图的含义

- S 节点存实际碱基；无分叉链可无损合成较长 S。节点不必等于 motif，也不按固定 16 bp 机械切片。
- L 的 `0M` 表示相邻片段，无重叠，不表示未知序列缺口。直接越过骨架片段的边描述所选比对下的跳过；不自动赋予生物学删除事件。
- 每条 P 对应一种完整观测序列；相同字面序列合并分子数，`members.tsv` 保留来源。重复单元不会全局并成一个自环节点，长度与顺序保留在 P 中。
- 这是以实际观测为骨架的 center-star 全局编辑距离比对图，使用 edlib 的一个最优比对。不是 PGGB 算法；未枚举重复区并列最优对应，拓扑不证明唯一同源关系。
- 只有 P 记录代表纳入的观测。图上其他可走通的边组合不一定被实际观测。

## 观测与计数

沿用输入 BAM 的参考方向，从单条 primary 记录提取目录区间，两侧各取 100 bp 参考锚点。完整覆盖、MAPQ ≥20 且不为255、两侧局部锚点匹配度 HiFi ≥0.90 / ONT与ONTUL ≥0.80 才纳入；排除 unmapped、secondary、supplementary、QC fail、duplicate、参考 N 跳跃和锚点次序冲突。阈值是第一版可审计筛选规则，不是测序错误校正或准确度估计。

分子键沿用 `[sample, namespace, read_name]`，相同分子跨 BAM/平台重复记录只计一次。同分子出现不同合格完整序列时整体排除并记录冲突。`excluded.tsv` 的数量是片段/无片段行数量，不是去重分子数；各排除原因可重叠。HiFi 与其他平台小计也可能重叠。

`--structures-json` 可直接导入已有结构输出的 rows，筛选条件相同；该模式不重新检查原 BAM。真实试运行采用 BAM 输入，并保存 `records.jsonl` 以回查 SAM SEQ。

这些条件只描述当前 BAM 中的完整区间观测，不证明全基因组唯一定位或真实等位型。输入目录的 `BOUNDARY=confirmed` 是沿用的声明。不会补回上游已删除的 secondary 记录，不输出基因型/CN。

## 无损压缩与有损结构

默认 `--loss-budget 0`，所有差异保留：SNP、CAG/CAA 等中断、插入删除、长度及单元次序。基础图先保留单碱基节点，再只合并没有分叉/汇合、没有 P 起止屏障的链；逐条验证拼接结果。

可选有损模式必须同时提供 `--merge-rules rules.json --loss-budget 0.01` 等明确预算。第一版仅允许等长 A/C/G/T 替换：源/目标都是已观测完整序列，逐位置许可、禁止替换链或环、逐 path 限额。目标不是自动构造共识。未获得具体保护规则前不自动舍弃真实位点的 SNP 或中断。

```json
{
  "schema_version": 1,
  "locus": "D1S80",
  "rules": [{
    "source_sha256": "从 paths.tsv 读取源序列哈希",
    "target_sha256": "从 paths.tsv 读取代表序列哈希",
    "allowed_positions": [50],
    "reason": "此位置允许用于展示压缩的具体依据"
  }]
}
```

位置为重复区实际观测序列的 0-based 下标，预算为被替换碱基数/源序列长度；这是表示损失，不是测序错误率。有损后共用结构路线的原始 P 仍分别保留；`lb:i`、`losses.tsv` 与原始 FASTA 说明变化。规则制定者必须保护具有生物学意义的中断和块次序；程序不代替这种判断。

## 文件与静态图

| 文件 | 用途 |
|---|---|
| `*.base.gfa` | 原始碱基级图，全部合格 P |
| `*.exact.gfa` | 无损链压缩图，全部合格 P |
| `*.structure.gfa` | 结构图；默认与 exact 同拓扑、同序列 |
| `paths.original.fa` / `paths.represented.fa` | 原始序列及结构图表示的序列 |
| `paths.tsv` / `members.tsv` | 路径长度、结构组、计数、遍历及分子来源 |
| `nodes.tsv` / `edges.tsv` | 实际碱基、比对列、路径成员与分子支持 |
| `losses.tsv` / `excluded.tsv` | 人工替换与筛选排除原因 |
| `*.preview.gfa/dot/svg/png` | 默认前10条路径重建的静态子图 |
| `summary.json` / `validation.json` | 参数、代码哈希、重建核验与限制 |
| `verification.json` | 额外独立脚本核验：GFA→FASTA→SAM切片和分子计数 |

私有 GFA tags：`mc:i` 为独立分子键数量，`pc:i` 为原始 P 数，`vc:i` 为 P 访问数；L 的 `ps:Z` 列出经过路径；P 的 `lb:i` 为人工替换碱基数。不能用节点访问数直接推断分子数或 CN。

静态预览按 HiFi 分子数、总分子数、长度和哈希排序取前10条（可 `--preview-paths 5`），并使用相同规则重新压缩所选子图。**预览节点编号仅对应 preview.gfa**，不是全图编号。图中盒子给出节点号、长度和短序列/周期表达；长非周期串显示首尾省略，完整序列在 S 中。边上的短编号与顶部 P 图例一一对应，`1-3` 表示路径1、2、3；灰边是全部展示路径共享。默认显示全部小差异；没有使用 hover 或交互。图宽不按碱基比例，图中不直接打印每条长序列的全部碱基。

超过20个节点时，静态图分为由上到下的连续面板，优先在跨面板边较少的位置分行。跨面板边两端以 `to 节点 (panel ...)` / `from 节点 (panel ...)` 虚线接续框精确配对；接续框只用于排版，不是新增 GFA 节点。每个真实节点仍只画一次，每条边都在图中保留。另提供 `*.preview.linear.svg` 的完整连续线性布局。面板拼接使用 [Graphviz gvpack](https://graphviz.org/pdf/gvpack.1.pdf)；未安装 gvpack 时回退到线性预览。

## 运行

在已分配计算节点上运行，使用新的输出目录：

```bash
python scripts/vntrpathgraph_cli.py \
  --bam-dir /path/to/bamdata_by_locus \
  --reference /path/to/GRCh38.fa \
  --catalog catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed \
  --locus ATXN2 --out /path/to/new-run/results --samtools /path/to/samtools
python tests/verify_pathgraph_results.py /path/to/new-run/results/ATXN2
```

把位点改为 D1S80 可运行另一位点。`slurm/pathgraph_pilot.sbatch` 在冻结源码目录执行全部单元测试、合成 BAM 集成、两个真实位点和独立重建检查。依赖沿用 `requirements-anchors.txt`、samtools 和 Graphviz；未安装 Graphviz 时 GFA 仍生成、预览状态为不可用。试运行要求实际 SVG/PNG 成功生成。

验证通过只证明保存的序列与路径、计数一致，不构成测序正确性、边界、同源定位或完整基因型的独立验证。
