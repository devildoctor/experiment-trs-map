# 重复单元识别第一阶段（2026-10-19）

本次同时保留基础单元和复合单元，先修公共规则，再对齐主点位图与结构图。结果仍是观察序列的注释，不是已校准的基因分型。

## 处理逻辑

1. 保留原始序列。所有 token 都记录原序列的 `[start,end)` 切片，主方案和备选方案均能逐碱基还原。
2. 找候选：主图保留原有近似周期扫描，补充局部精确串联扫描；结构图的局部发现独立于既有模板，避免错误种子把真实重复吞成“近似匹配”。默认局部新候选长度为 2–200 bp；单碱基目录模板可注释。
3. 将精确幂拆为基础单元，例如 `ACGTACGT` 的基础单元是 `ACGT`。目录原单元保留为复合层。同一字典内能由更短模板完全铺满的单位也列为复合层；这一层级依赖字典，不意味着唯一的生物学定义。
4. 在原始 read 中检验候选的精确命中位置。相邻同族单元，或两侧均有紧邻单元的中间变体，才计为重复上下文证据。`AABAA` 中一个 B 可以由三条独立分子共同支持；三条 read 各自孤立出现一次 B 则不够。支持数按分子去重。
5. 目录模板允许参与注释，但 `usable` 与 `validated` 分开；没有观察证据时不再因为来自目录就自动验证通过。主图使用已提取的分子；结构图发现/验证使用 primary、MAPQ≥20、非 duplicate/QC-fail 的 HiFi 分子，两者证据范围保留各自原有策略。
6. 按固定模板相位做动态规划，主图与结构图保留总代价最低的两个路径。检查正反模板方向，但不在每一个单元自由轮换相位。`phase` 是相对于 canonical family 的偏移；结构图的 `orientation` 仍相对于目录模板，另外保存 `canonical_orientation`，避免把命名方向误作倒位。
7. 完整单元、编辑后的完整单元、端部残片、未知碱基分开计数。端部 1 bp 的弱匹配保持 UNKNOWN；完整率不包括 partial。A×4 按四个 1 bp 完整单元计数。
8. 复合层从主方案的完整基础 token 派生，禁止跨 UNKNOWN 或 partial。自动报告混合 token 模式；纯同类的 8/12/... bp 谐波不无限新增，只有目录明确指定的复合单位另行统计。不同复合定义可以重叠，次数不能相加。

例如 `ACGT×30`、目录同时给出 `ACGTACGT` 和 `ACGTACGTACGT` 时：基础层是 30；8 bp 目录块是 15；12 bp 目录块是 10。这是三种分别报告的描述，不能合计成 55 个单位。

## 输出变更

- `motifs.tsv`：增加层级、primitive sequence、相位规则、证据状态、分子数和支持位置。结构图的 `final_read_support` 是注释命中的 read 数，`evidence_read_support` 才是符合精确重复上下文的 HiFi 支持数。
- 主图 `tokens.tsv` 仅包含基础层；复合信息单独写入 `compound_blocks.tsv`，避免与基础计数相加。旧代码若读取 compound token，应改读新表。
- `annotation_alternatives.json`：主/备选代价、分差、同分歧义、完整/partial/unknown 碱基统计以及备选的实际序列切片。
- 主图 SVG 显示主/备选完整单元数和分差；结构 HTML 显示证据状态、完整率及可展开的备选路径、复合块。原始 GFA 保存的主路径依然能精确还原序列。
- read-view 修复基础幂、1 bp 单元、弱端部匹配和完整率，并保留平台特定评分；双路径展示本阶段限于主点位图和结构图。

## 运行与验收

```bash
export PYTHONPATH="$PWD/.deps/vntr-anchors"
/home/data/user/stzhu/miniconda3/bin/python3 -m unittest discover -s tests -p 'test_*.py' -v
/home/data/user/stzhu/miniconda3/bin/python3 tests/run_phase1.py \
  --out runs-result/<新目录>/acceptance \
  --samtools /home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
```

新增 `test_phase1.py` 共 14 项验收。`run_phase1.py` 保留 4 组示例的 SVG/GFA/HTML/JSON、验收表，并调用三个已有端到端程序，保留它们的合成 FASTA/SAM/BAM 和输出。主/备选重建、分子去重、单锚定/部分覆盖、不跨未知区的复合计数均有断言。

## 第一阶段范围

局部精确扫描依然不能保证找到含高错误率、无连续精确拷贝的新单元；原有主图候选数量上限仍为 32（目录更大时按目录大小）。备选是当前候选字典和评分下的第二条路径，不是所有可能解释。分数差不是概率，仅同分被自动标为歧义。不同入口保留原有平台和编辑代价，噪声 read 的结果仍可能不同。真实数据的阈值校准、统一概率模型和更全面的局部近似发现留待后续阶段。
