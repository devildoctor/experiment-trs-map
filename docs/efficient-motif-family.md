# 有效字典、HMM 与可展开家族图

入口 `scripts/run_efficient_family.py` 按位点类型分流。VNTR 运行完整字典 HMM 注释、q=0.1 有效字典选择、有效字典 HMM 重注释、逐 P 3% 有损归并及展示家族选择。STR 沿用冻结完整字典、`seed_dp` 和原有无损结构图；不进入 HMM/SKIP、精简字典或新家族流程。边界和读段筛选沿用现有输入。

模型可显式声明 `repeat_type: STR` 或 `VNTR`。为保持旧冻结字典文件及其导入注释哈希，现有模型不改写：内部模板均为 1–6 bp 时走 STR，均大于 6 bp 时走 VNTR；混合长度模型须显式声明。末端模板不参与分类。ATXN2 因而走 STR，D1S80 走 VNTR。

STR 使用模型旁或内置目录中的 `<locus>.structure-policy.v1.json`，保留 CTG/TTG 等独立单元及中断类型，旧版最大改写为 0 bp。`--q` 和 `--loss-fraction` 仅作用于 VNTR。即使主结构图入口对 STR 传入 HMM policy，也转入旧版 STR 输出；专门的 HMM 重建入口只接受 VNTR。

## 字典选择

`vntreffmotifs.py` 按固定版本 [vamos 的 effMotifILP.py](https://github.com/ChaissonLab/vamos/blob/3dd3b429922399e9f14e9bfebc0ea3056b18fd0a/buildTRs/lib/effMotifILP.py) 的预算、频率约束和目标函数选择代表。当前两个位点只有 2 或 9 个内部模板，使用子集枚举求解，无需增加 OR-Tools 依赖。

统计口径是各登记 P 的完整字典单元出现次数，包含参考路径；不是分子加权频率或群体等位频率。模板按出现次数降序输入预算计算。每种模板统一映射到出现次数不低于它的代表；优化 `delta * 代表数 + 按次数加权的模板编辑距离`，且后者不超过 delta。q 控制模板替换预算，和实际逐 P 碱基改写比例独立。

完整字典及其注释保存在输出中。有效字典作为新的版本传入 HMM，重注释仍逐碱基保留输入字符串。RT 等独立末端模板不参加内部字典压缩；已有中断类型只要保留在有效字典中，仍沿用现有有损归并保护规则。

## 家族与改写

`vntrfamilies.py` 搜索家族划分，每个成员到其代表的直接编辑距离最多 1 bp，避免相似关系沿链条扩大。家族共享限定于相同比较位置、方向、角色的单部分内部节点，各候选节点不能出现在同一登记 P 中。原类型、实际字符串、有序变体、weak 和 skip 参数保存在每次访问中。

3% 是每条 P 的实际修改上限，不要求耗尽。`plan_edits` 的可选展示构建器使每次候选修改按参数化家族图的实际 N+E 接受或拒绝。完成改写后再次选择家族，并与改写前固定家族规则比较，保留更小的图。该流程不声称对字典、分解、改写和家族联合求得全局最优。

## 运行及输出

```bash
python scripts/run_efficient_family.py \
  --input /path/to/source_members.json \
  --model catalogs/D1S80.structure.v1.json \
  --out /path/to/new_results --q 0.1 --loss-fraction 0.03
```

ATXN2 的旧观察路径输入若没有独立参考 P，需要另传 `--reference /path/to/GRCh38.fa`。程序调用已有的参考路径分配函数，补入参考并保存新旧 P 编号映射，不删除任何观察记录。

- `baseline_3pct`：完整字典、weak/skip 参数化、3% 改写对照。
- `expanded_3pct`：有效字典及本轮选定的改写，家族展开状态。
- `family_3pct`：相同字符串和账本的家族图。
- `comparison.tsv`：各阶段节点、边及实际修改量。
- `dictionary_selection.json`、`family_selection.json`：实际选定的字典和家族。
- `<locus>.families.html`：离线查看家族/展开图，按 P 高亮并查看各节点的单元原文。

STR 输出在 `str_legacy/`，其中 `index.html`、`<locus>.preview.svg`、`<locus>.structure.gfa` 和 `paths.tsv` 沿用旧格式。根目录 `processing.json` 记录实际使用的流程。

家族与展开状态的 FASTA、GFA 和改写账本必须逐字节一致。图上的任意拼接路线不是新增观测；登记路径及其访问参数用于序列还原。
