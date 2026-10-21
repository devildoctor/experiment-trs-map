# STR 单元图与压缩路径显示

仅对已启用 `STR_STRUCTURE`、目录单元为 2–6 bp 的位点生效。ATXN2
使用现有 `grch38_str_expression_atxn2_v1.bed`；>6 bp VNTR 继续原流程。
不变更坐标、候选发现、支持门槛、分子分配或分型状态。

## 图语义

`scripts/vntrstrgraph.py` 使用现有精确块注释生成图：

1. 去掉重复块计数后，按 `(类型, 碱基)` 构建前缀树。
2. 仅合并 token、后继状态及终止状态完全相同的后缀状态。
   两个位置均为 CTG 不足以合并；同一路径被中断隔开的 CTG 保持不同节点。
3. 重复节点的 `S` 序列只有一个 2–6 bp 单元，重复访问对应 `0M` 自环。
   固定中断、N、尾部残片保留完整碱基。
4. 每条 `P` 逐次写出节点，拼接后等于候选序列加左右各 25 bp 参考上下文。
   `REFERENCE` 是参考路径；样本路径仍是候选，不能当作确认的二倍体分型。
5. 可选用户标签 `rp:Z:L+,N1+*8,...` 保存可读的连续访问次数。
   它不是标准 `P` 第三列的替代语法，通用查看器可能忽略该标签。
6. 静态路径图将连续访问画成一个节点和一个自环，旁标 `×8`。
   此处 8 是该路径中的单元次数；自环实际走了 7 次，不再把 7 画到图上。
   不同路径可在同一节点有不同次数，因此总图不标固定拷贝数。

同一分子对一个节点或一条边只计一次支持。循环次数不会放大支持数。
拓扑允许一些没有显式 P 的计数或组合；报告序列时必须使用保存的 P。
图保持参考正链 CTG；反向互补 CAG 表达仍保存在 `str_paths.json`。

## 文件规范

| 文件 | 内容 |
| --- | --- |
| `ATXN2.gfa` | 新主图，GFA 1.0 S/L/P；完整路径及可选 rp 标签 |
| `ATXN2.json` | 与主 GFA 配套的节点、边、路径、次数、序列哈希及分子证据 |
| `ATXN2.svg` / `ATXN2.dot` | 全部节点及边；Graphviz 布局，同 ID 同碱基 |
| `ATXN2.paths.svg` | 所有显式路径，每行一条，显示该路径的 ×n |
| `paths-svg/ATXN2.<path_id>.svg` | 一条路径；如 `ATXN2.LCL5_C13.svg` |
| `ATXN2.detail.gfa` / `ATXN2.detail.json` / `ATXN2.detail.svg` | 原点图及其配套证据，成套保留 |
| `str_paths.tsv` / `str_paths.json` | 参考及全部候选的表达，包括未通过支持门槛的候选 |
| `graph_validation.json` | 独立验收脚本的结果 |

主图和逐路径 SVG 只显示节点 ID、实际碱基、连接和路径次数。颜色：蓝色
重复单元，金色字面碱基，灰色参考两翼。框宽不表示基因组长度。
支持数、歧义、平台差异和候选状态可查询 JSON/TSV。点图分支不新生成 HTML；
既有 read/structure 模式仍按其原有规则输出。

## Slurm 与验证

渲染依赖 Graphviz 的 `dot`。测试和运行都在 Slurm：

```bash
sbatch --export=ALL,RUN_DIR="$RUN_DIR",BASELINE_DIR="$BASELINE_DIR",BEFORE_RESULT="$BEFORE_RESULT",OUT_DIR="$OUT_DIR" \
  "$RUN_DIR/slurm/vntrstrgraph.sbatch"
```

先创建新的 UTC 时间戳运行目录 `<UTC>_str-gfa_ATXN2`，保存 baseline 和
after 源码及 SHA256SUMS；结果写入该运行的 `results/pointgraph/ATXN2/`，
日志放 `logs/`，作业号记入 `jobs.tsv`。已有非空点图目录仍拒绝覆盖。

验收包括原有测试、STR 图语义测试、合成 BAM 闭环、16 bp VNTR 修改前后
逐字节比较、实际 ATXN2 全路径重构、独立解析 GFA/rp、SVG 节点/边/次数
逐项对照及原始候选/分配结果一致性。退出 0 不替代人工图形检查。

格式依据：[GFA 1 官方规范](https://gfa-spec.github.io/GFA-spec/GFA1.html)。
本次没有运行 PGGB 或 ExpansionHunter。
