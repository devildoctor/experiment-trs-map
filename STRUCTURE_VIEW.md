# 共享重复结构与全部 read 证据

以下参数已按当前源码核对；开头 r3 问题保留为历史背景。

## r3 已确认的问题

- 每次出现都使用 read 特有节点 ID，全部三个位点的图中分支节点数、自环数都是 0；RC 被写死为 1。这是独立链集合，不是共享结构图。
- PER3 只有一个 54 bp 注释模板。r3 中未解释碱基占 67.9%，HiFi 也达 68.69%。不能将该结果解释为不存在重复。
- UNKNOWN 每碱基罚 2、模板编辑每次罚 1.5，鼓励模板吸收相邻碱基。首条 PER3 HiFi read 在提取片段 [55,109) 有精确 54 bp 模板，却被注释成 [55,116) 的 61 bp 单元。
- ABCA7 的 21 类模板用 6 色表示，没有连续重复压缩、共享节点支持数或明确的复合排列展示。

## 改进依据

- TRGT 的 BED 使用 `ID`、`MOTIFS`、`STRUC`。本程序使用同样的种子字段，`STRUC=<TR>` 作为说明，不声称实现 TRGT 分型或完整结构语法。
  https://github.com/PacificBiosciences/trgt/blob/main/docs/repeat_files.md
- TRGT 的 waterfall 图按 read 原有顺序显示重复片段。本程序借鉴该视图，不将不同起点的部分 read 假装对齐到同一参考位置。
  https://github.com/PacificBiosciences/trgt/blob/main/docs/trvz-plots.md
- MotifScope 支持已知 motif 引导与单元连续次数表达。本程序借鉴表达方式，未复制或调用该项目的算法。
  https://github.com/holstegelab/MotifScope

## 种子目录

`catalogs/grch38_structure_seeds_il1rn_per3_abca7_v1.bed` 仅包含 IL1RN、PER3、ABCA7，分别给定 1、2、2 个候选种子。
IL1RN 为 r2/r3 已使用的 86 bp 模板；PER3 为原 54 bp 模板及多条 r3 HiFi 读段观测到的另一 54 bp 单元；ABCA7 为 r3 模板中的 25/26 bp 单元。
这些是数据支持的候选种子，不是已校准的参考结构或病理定义。候选边界保持不变，不将原 BAM ±500 bp 简单当作 SEQ 剪切长度。

## 处理流程与关键参数

1. 读取全部 BAM 记录，沿用单侧锚点/CIGAR 提取，无 MAPQ、支持数和双锚数量的 read 删除门槛。所有只覆盖侧翼或 SEQ 缺失的记录也保留状态。
2. 用 BED 种子和参考序列补充初始字典，再迭代注释 HiFi、发现候选模板。只从合格 HiFi 主比对（MAPQ≥20，排除 duplicate/QC-fail）学习；候选需精确邻接重复上下文、默认跨样本合计至少3个独立分子，并在至少一个样本达到10%比例；这里不要求单一样本达到3条。`--variant-support`、`--min-sample-fraction` 可调整。默认 `--max-variants=256`，超过上限报错；`--discovery-rounds` 默认3，只允许1–3轮，无新增即提前结束。这些门槛控制字典学习和结构支持展示，不删除 read。
3. 全部 HiFi/ONT/ONTUL 用扩充字典重新注释。正反模板分别测试；UNKNOWN 保留实际碱基；末端少于 2–8 bp（由 motif 长度决定）的短匹配不硬称部分重复。
4. 评分为 UNKNOWN=1/碱基，完整单元 `0.1 + min(3, 0.6×模板长度)×编辑距离 + 2×长度差`，精确末端残片代价0.2。长度1模板不允许编辑，长度2–6允许至多1次编辑；更长模板的候选编辑上限为 HiFi 32%、其他平台38%（向下取整且至少1）。长度偏移另有限制，最终由总代价选路。这些值是搜索与评分参数，不是准确率或有损替换比例。保留主方案及一个次优方案，主、备选 token 的同分歧义相互比较。
5. 精确序列层按完全相同的实际碱基合并，GFA 路径必须重构原片段。结构层按已注释 motif/候选变体/方向合并；未知序列仅按精确相等合并，不用一个 UNKNOWN 生物节点冒充所有未知结构。
6. 支持数统一按 `[sample, namespace, QNAME]` 计数，主图、read-view、结构图与summary一致。namespace优先取 `--sources`，否则取记录 RG 对应的 LB，缺失时为空。一个分子内的多次重复、补充比对、跨BAM或平台标签副本不增加支持；经过次数另列。不同样本或不同namespace的同名读段分别计数。旧行缺少namespace时兼容为空，不能凭历史文件恢复原先未保存的LB。
7. 连续同类单元显示 `M1×n`；观测到重复的 2–6 单元组合显示复合排列候选；内部未知片段和反向 motif 匹配分别标记。反向 motif 匹配不是已确认倒位，不自动宣称缺失、插入等变异。

## 输出与查看

- 首先打开 `位点.structure.html`：BED 字典、全量家族概览图、全部 read 的同尺度条带图、家族/变体切换、连续次数、全部节点和支持数、结构模式分组。
- `位点.overview.svg`：结构类型与局部模块摘要。`位点.family-transitions.svg` 保存家族层共享图；U/P 仅是显示桶，不能解释为相同未知序列。
- `svg/index.html`：每页 40 read，保留全部分页；长路径完整碱基见 HTML 与 TSV，SVG 概要长标签带完整 title。
- `位点.reads.gfa`：精确合并节点，每条 read 路径可重构真实碱基。
- `位点.structure.gfa`：抽象节点 `S=*`、模板存储于 TS 标签，不能用于重构原始碱基。
- `read_nodes.tsv`：两个图之间的逐节点映射、实际坐标、方向、编辑距离与原始碱基。
- `structure_nodes.tsv`、`structure_edges.tsv`：read 支持数、主比对 MAPQ≥20 支持数、经过次数、分样本支持。
- `patterns.tsv`：相同节点顺序且覆盖状态相同的观察模式；不是确定的等位分型。
- `status.json`、`parameters.json`、`verification.json`：数量、未知比例、源码与 BED 哈希、独立重构及相对于 r3 的 read/碱基保留检查。

## 视觉检查约定

这是科研图程序的专用输出，不是通用业务报告。主比较为 read 之间的重复单元顺序、次数与中断；采用 motif 路径、共享图和按实际 bp 的条带图。家族使用蓝/橙等固定颜色，变体保留文字 ID，未知使用灰色，部分单元虚线。所有 read 保留，分页不是抽样，固定全局 bp 比例尺。

## 验收与限制

每个 Slurm 作业先执行全部单元测试和合成 BAM 集成测试，再处理真实 BAM，最后独立验证精确 GFA 重构、支持数去重、与 r3 全部 read 清单和提取碱基完全一致。
候选窗口的边界问题仍保留。局部模板注释可能歧义，生物结构判断仍需参考校准及多 read 支持。共享图能够组合出从未被真实 read 观察到的路径，因此只允许把实际保存的 read 路径作为证据。不能为了得到预想图形而制造分支、复合重复或倒位。
