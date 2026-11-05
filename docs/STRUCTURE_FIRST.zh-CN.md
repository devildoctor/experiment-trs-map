# 结构整理 → 有损归并 → 按结构建图

新增 `vntrpathgraph_cli.py --graph-mode structure`，沿用原有完整观测筛选、分子身份和冻结划分算法。旧的 sequence 模式保持原有参数含义。

## 使用

在已分配的计算节点运行：

```bash
python scripts/vntrpathgraph_cli.py \
  --bam-dir /path/to/bamdata_by_locus --reference /path/to/GRCh38.fa \
  --catalog catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed \
  --locus ATXN2 --out /path/to/new-output --samtools /path/to/samtools \
  --graph-mode structure --locus-model catalogs/ATXN2.structure.v1.json \
  --structure-policy catalogs/ATXN2.structure-policy.v1.json --max-loss-bp 2
```

D1S80 使用同名位点模型和策略。`--max-loss-bp` 默认 0，首版支持 0、1、2；是整条 path 的总编辑碱基数。结构模式不使用旧的比例参数 `--loss-budget`。`--structures-json` 可导入同一模型哈希下的既有结构注释；BAM 模式对每条去重完整序列调用现有 `seed_dp_v2`，使用固定的 HiFi 代价配置以避免平台混合改变同一序列的解释。

## 结构和不确定性

- 所有冻结 motif 身份、顺序和次数均受保护；连续段和复合模式直接显示在节点。
- 重要单元类型不会因只差 1 bp 而归并，例如 ATXN2 的 CTG/TTG。
- 未知和残片保留原始碱基，不能变成跳过边。结构模式记录种子候选 DAG 的全部最优转移，通过路径计数确定所有最优划分共同经过的单元，只有其间未定区段收成未定块；不会把两个保存的解释误当作全部可能。没有该局部证据的旧输入仍整条标为未定。划分解释和歧义区间保存在 `path_features.jsonl`。
- 单元连续段对骨架做类型比对。两端相同完整段按固定比较约定共享；中间存在多种最优对应时保留整个比较块，不强选精确删除位置。
- 前后不同位置的相同单元不能全局合并成一个自环。分支、汇合、path 起止和未定区边界处停止链压缩。
- 共同节点的表示碱基完全相同。节点表达可展开为有序类型串；真实碱基仍在 GFA S 中。
- 相位、边界声明和对齐位置均来自固定模型/比较约定，不构成生物学唯一性证明。单元数不报告为基因型或已验证 CN。

## 有损规则

只有完整、未报告并列/冲突的相同受保护类型串进入归并。按既有 HiFi 分子数、总分子数、长度和哈希排序选择实际观测的固定代表，每条成员直接与最终代表比较，不做传递归并。

允许的位置在策略文件中明确给出：当前试运行只开放冻结模板共同具有的内部碱基列，保护亚型区分列、首尾单元、专用端部。修改过的单元在代表序列中必须精确匹配模板。重复出现的编辑亚型保留，避免抹掉交替/成块规律。这个策略是显示取舍，不证明获准隐藏的变化在生物学上不重要。

在对应单元内部计算编辑脚本。若任何最优脚本触及禁止位置则拒绝。S 表示替换、I 表示向表示序列添加、D 表示从原始序列移除：

```
loss_bp = S + I + D
表示长度 = 原始长度 + I - D
```

原始与表示坐标分别保存。GFA 重建表示序列，表示序列和账本恢复原始序列。代表原始序列的精确支持与归并路线支持分别输出。

D1S80 模型固定沿用旧运行的 9 个内部基本模板，不从本轮样本学习。RT 末端模板只在旧证据支持的 LCL5–LCL8 样本、完整区间的精确右端匹配时使用，不加入内部词典；这表示与冻结端部模板相容，没有重新宣称独立末端支持。模板来源及哈希保存在模型中。

## 输出和验证

`*.exact.gfa` 是零损失结构图；`*.structure.gfa` 是按预算得到的表示图，保留全部原始 P。`*.preview.gfa/svg/png` 显示前 10 条原始 P，节点编号沿用全图。图中的路径编号对应实际 P，其他可拼接的路线不增加为观测。

`paths.tsv`、`path_features.jsonl` 保存结构、原始/表示长度和注释状态。`merge_decisions.tsv` 保存候选决策；`losses.tsv/jsonl` 保存双坐标编辑脚本；`structure_groups.tsv` 分开列出路线分子数和代表序列精确分子数。原始 FASTA、source_members 和 SAM 提取记录不丢弃。

GFA 私有标签 `lb` 保持替换数含义；新增 `le` 总编辑数、`li` 添加数、`ld` 移除数、`ol/rl` 原始/表示长度。未知片段和对齐歧义有独立标签/说明，不以节点访问次数推断分子数。

`slurm/structure_first_pilot.sbatch` 只运行 pathgraph、冻结划分的相关测试、新增的 11 个核心用例和 ATXN2/D1S80 真实 BAM，不运行全仓库测试。每个位点由独立脚本重新读取 GFA、核对存储的 SAM 切片、原始/表示序列、残差恢复、分子身份及静态 SVG 的节点/边，并与上一版完整观测集合比较。
