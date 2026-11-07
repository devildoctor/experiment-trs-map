# 左右侧翼 v2：显式运行与复算

本分支新增完整侧翼定位和冻结外围边界发现。无 `--flank-boundary-model` 时保持旧提取流程。
核心 P 仍是原始 read 核心子串的精确折叠；外围上下文另存，不改变核心字典，不代表基因型。

## 环境

新模式使用 Python 3.13、NumPy 2.2.6、Numba 0.62.1、llvmlite 0.45.1，版本见
`requirements-flanks.txt`。安装到独立环境或独立目标目录，保留原有 edlib、samtools 等依赖。
本次服务器依赖位于 `runs-result/20261108_flank_v2/deps`；不修改系统 Python。

```bash
PROJECT='/home/data/user/stzhu/重复区域分析/vntr-point-graph'
RUN="$PROJECT/runs-result/20261108_flank_v2"
cd "$PROJECT"
export PYTHONPATH="$RUN/deps:$PROJECT/.deps/vntr-anchors"
export NUMBA_CACHE_DIR="$RUN/numba-cache-manual"
PYTHON=/home/data/user/stzhu/miniconda3/bin/python3
SAMTOOLS=/home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools
```

以下计算命令应放入 Slurm 脚本提交到计算节点。输出目录必须是新目录。

## 使用已冻结模型

以 ATXN2 W1 为例；本轮所有 W 窗口得到相同边界，W1 仅是较简配置，尚未获得默认切换验收。

```bash
"$PYTHON" scripts/vntrpathgraph_cli.py \
  --bam-dir "$PROJECT/bamdata_by_locus" \
  --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --catalog catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed \
  --locus ATXN2 --out "$RUN/manual-new-run" --samtools "$SAMTOOLS" \
  --flank-boundary-model "$RUN/final_replay/models/ATXN2/W1.json" \
  --graph-mode structure --locus-model catalogs/ATXN2.structure.v1.json \
  --structure-policy catalogs/ATXN2.structure-policy.v1.json --max-loss-bp 0
```

`scripts/vntrstructure.py` 接受相同的 `--flank-boundary-model` 参数。
从其输出的 `*.structures.json` 导入点位图时，也必须提供同一个边界模型、参考和目录；
模型哈希或片段提取版本不同会拒绝导入。旧结果不带新模型，继续按旧语义处理。
不将旧缓存静默标记为新结果。

## 重新发现外围

发现输入是保存全部 BAM 原记录的 `records.jsonl` 及其同目录 `summary.json`，不是已保留的 P。
只使用 LCL5–LCL8 的 HiFi 主比对、MAPQ 50–254；同分子重复去重、冲突排除。
完整读取给定的位点 BAM，在固定取数窗口内计算覆盖、变异和归一化。

```bash
"$PYTHON" scripts/vntrflankbounds.py \
  --records "$RUN/comparison/B0/ATXN2/records.jsonl" \
  --reference "$PROJECT/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  --catalog catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed \
  --locus ATXN2 --out "$RUN/rediscovered-ATXN2"
```

固定模型参数来自 vclust 的已锁定提交，仅用于上述 PacBio HiFi 发现集；来源与原许可见
`third_party/vclust/`。本项目独立实现 read 定位，不运行 TRGT/vclust 程序。

## 本轮运行组织

- 旧源代码：`runs-result/20261108_flank_v2/baseline_source`，基线提交记录在 `baseline_commit.txt`。
- 第一轮：`comparison/B0,B1,B2,W1,...,W8`；两个位点各 11 组。
- 最终代码重放：`final_replay/`；B0 链接到同一只读基线，其他 10 组在新目录完整重算。
- 冻结模型、逐位置发现审计：各轮 `models/`。
- 共同分母、敏感性、逐分子留一：`audit/`。
- 详细补充表、隐藏方法身份的病例材料：`report/`。
- 所有日志：`logs/`。早期源快照和结果保留，不覆盖。

`run_flank_comparison.py` 运行已预先复现 B0 的十个新配置；
`audit_flank_results.py --sensitivity` 运行固定目标、13 种搜索域设置的对照；
`audit_flank_leave_one.py` 对固定发现集逐分子留一；
`report_flank_evidence.py` 从已有结果生成明细，并核对最终重放一致性。
相应 Slurm 脚本放在 `slurm/flank_*.sbatch`。整轮重算请使用新的 run 目录，并调整脚本中的 RUN。

## 固定语义和限制

主锚 250 bp，错配代价 2、gap 为 5+k；完整参考锚、query 两端免费；每侧至少 175 个匹配。
只在原 CIGAR 提示附近 ±120 bp 完整枚举端点，超过记录的 DP 上限记未完成。
保留最低总代价的全部有序切点组合，再按匹配数、原提示距离、query 坐标确定规范代表。
规范代表不消除原始歧义，也不证明全基因组唯一定位。

核心和外围分别定位。`core_only` 可在其余质量条件合格时进入 P；
`analysis_only`、核心未完成或规范代表包含冲突不能进入完整核心 P。
歧义核心可形成带标记的观测 P，但不能升级为结构摘要的 strong 证据。
未知/残片、原 CIGAR 覆盖与原始序列保留。新锚不会再经过旧 0.9/0.8 评分门槛。

所有主比较关闭有损合并；每条 P 必须逐字重建。新增分子、更多 P、边界留一稳定，均不能
单独证明误截减少。Q17 还要求独立证据支持的旧误截修复、无新增确认误截及同目标敏感性不恶化。
未满足时保留旧默认；不合并到原分支、不替换旧结果。
