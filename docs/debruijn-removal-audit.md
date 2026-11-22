# de Bruijn 图删除后的使用与核验（2026-11-17）

当前分支：`feat/efficient-motif-family-3pct-20261111`。
删除提交：`fa05fa03af3734d1f5163247c6209ac604afb039`；对照版本：`9c5739f18adf31add56b9fb36752a5eeb527ad89`。

## 当前入口

| 功能 | 入口 / 说明 |
|---|---|
| 候选点位图 | `python scripts/vntrgraph.py run ...`，见 [README](../README.md) |
| read / structure 视图 | `scripts/vntrgraph.py --read-view ...` / `--structure-view ...` |
| 完整观测路径图 | `scripts/vntrpathgraph_cli.py`，见 [路径图](PATH_SEQUENCE_GRAPH.zh-CN.md) |
| B2 全流程 | `scripts/run_b2_pipeline.py`，见 [B2](flank-b2.md) |
| VNTR HMM/SKIP、有效字典与家族图 | `scripts/run_efficient_family.py`，见 [分流与字典](efficient-motif-family.md) |
| 纯序列 16 bp 注释库 | `scripts/vntrunits.py:primitive_units`；无独立 CLI |
| 终端证据库 | `scripts/vntrterminal.py:validate_terminals`；无独立 CLI |

STR 保留原有 seed_dp 结构流程；HMM/SKIP 与有效字典用于 VNTR。
路径结构图中的 compact_graph 和 graph_support 仍用于路径压缩及分子支持，
不是已删除的 k-mer 图模块，不应按名称继续删除。

## 删除带来的预期变化

- 不再提供 `--dbg-view`、多 k 图、联合分支分析及其 Slurm/验收脚本。
- 不再生成 `*.dbg.json`；候选/audit 表不再含 `dbg_k`、`graph_support`。
- 不再依据局部 k-mer 连边覆盖率生成 `low_local_graph_support`。
  如果它曾是某候选唯一的不确定性原因，该候选的 uncertain 及依赖该标记的展示/计数可能变化。
  这不表示新增了生物学分型证据。不能声称删除前后所有元数据逐字节相同。
- 图构建以显式保存的路径为准；图上任意游走仍不代表已观测等位序列。
- 历史分析 CLI 的来源哈希和参考一致性检查不会因保留 Python 库而自动执行。

## 残留引用审计

当前 scripts 中没有导入已删除模块；现行 CLI 和 Slurm 中没有旧图入口。
测试中的旧输出名称是用于断言它们不再出现的负向检查。
README 的删除说明是兼容性记录，不是可用入口。

2026-10-29 的代码审阅、修复记录、旧代码解读和
`VALIDATION_20261030.json` 是历史快照，允许保留旧模块、命令和哈希。
旧代码解读整体不作为当前调用链说明；历史文档的源码链接已固定到删除前提交，
不会再跳转到当前分支中不存在的文件。历史 JSON 的原始哈希不改写。
终端单元文档已移除当前用法中的旧 CLI/Slurm 指令，改为库接口说明。

## 服务器验证

在服务器现有 CPython 3.13、B2 NumPy/Numba 依赖与锚定依赖环境中运行：

- 全量 unittest：293 项通过。
- 其中 HMM/SKIP 8 项、weak/SKIP 参数化 5 项、有效字典/家族图 7 项、
  路径图 18 项、结构图 11 项、STR/VNTR 分流 4 项通过。
- `tests/integration_synthetic.py` 通过：合成输入经 samtools 生成 BAM，
  点图 3 条路径逐碱基重建、分子去重/平台支持、SVG 和状态均通过；
  新增断言确认无旧辅助图文件、字段或标记。
- `tests/integration_pathgraph.py` 通过：B2 提取、重复 BAM 去重、
  参考加两条观测路径、5 个独立分子、GFA/TXT 重建及 PNG 渲染通过。
- 与删除前源码比较：HMM/SKIP、路径图和结构图模块未修改；
  vntrcore 保留的函数/类（除 Allele 删除字段外）AST 一致。
  拆出的 template_for、primitive_templates、primitive_units 三个函数 AST 一致。

本次验证包括真实 BAM 文件格式的合成集成测试，未重跑全部真实样本数据；
测试通过不等于已经证明所有真实输入的结果完全不变。
