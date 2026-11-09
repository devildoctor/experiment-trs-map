# B2：固定核心与左右 250 bp 侧翼

2026-11-09 起，点位结构图、结构分解和 all-read 入口统一自动使用 B2。
无需指定边界模型，也没有 B0、B1、多尺度扩边或 100 bp 诊断选项。
旧实验结果继续保存在 runs-result/20261108_flank_v2，历史源码保存在 Git。

核心区间来自 catalog；左右各取 250 bp 参考模板，要求每侧至少 175 个匹配碱基。
参考模板完全消费，read 两端自由；错配代价 2、长度 k 的 gap 代价 5+k。
CIGAR 只提供搜索位置，搜索 padding 为 120；未覆盖模板时不外推。
左右候选按总代价、匹配数、与 CIGAR 距离及坐标排序，保留全部最低代价同分候选。
保留歧义标记；规范代表不是独立的真实边界证明。

P 入口必须通过 B2 完整核心证据及原有质量/分子冲突检查。
部分 reads 仍保留原始记录、CIGAR 显示片段和失败原因，不作为完整 P。
原始核心字符串用于 P 折叠；固定词典、结构策略与原有损失预算不变。
命令输出 boundary_model.json、extraction_audit.jsonl 和全部原始 records.jsonl。
导入 structures.json 时也必须提供 reference/catalog，并核对 B2 模型身份。
旧版结构 JSON 需从现有 BAM 重新提取，不能绕过 B2 证据准入。

完整双位点复核入口为 scripts/run_b2_pipeline.py；Slurm 模板为 slurm/b2_pipeline.sbatch。
主流程：全部 BAM 记录 → B2 → 按核心原始字符串汇总 P → 冻结词典与统一 HiFi 评分分解 → 结构图。
额外审计运行逐 read 分解及导入建图；该旧入口保留平台相关分解评分，因此节点划分可能不同。
主图保存于 graphs/，逐 read 分解保存于 structures/，导入交叉核对保存于 import_audit/。
核对旧 B2 每条记录的切点、候选集合和核心序列；核对所有 P 的分子成员以及主图完整 GFA。
导入审计要求原始路径序列和成员一致，并独立通过精确重建，不要求不同分解评分的节点划分一致。
默认测试流程使用 max-loss-bp=0，所有路径必须精确重建。
