# vntr-point-graph

**2026-11-17：已移除 de Bruijn 图功能。** 不再提供 `--dbg-view`、多 k 图与联合分支分析入口，主点图不再输出 `*.dbg.json`、`dbg_k`、`graph_support` 或 `low_local_graph_support`。现有分析结果未删除。纯序列单元注释保留于 `scripts/vntrunits.py`，终端证据库保留于 `scripts/vntrterminal.py`。下方历史分支及旧文档记录仅描述对应历史版本。

**观察路径的 B2 取段。** 在 BAM 观察路径、结构分解和 all-read 入口，按 catalog 核心和左右 250 bp 模板取段。普通候选点图使用自己的锚定、候选聚类及分配流程。B0/B1/W1–W8 和 100 bp 诊断已移除，历史结果保留。完整流程见 [B2 使用与核对规则](docs/flank-b2.md)。B2 最初实现于历史分支 `feat/b2-only-20261109`；下方早期分支说明为历史记录。

面向长读段 VNTR / STR 的可审计点位图分析工具。输入按位点整理的 BAM、匹配参考及位点目录，输出明确保存的序列路径、分子证据、重复单元注释和 GFA/SVG。

主要功能：

1. 依据两侧锚点提取实际读段序列，按独立分子去重，保留定位与候选冲突。
2. 生成候选序列点图，并提供可选 STR 单元图和 VNTR 路径结构图。
3. 支持重复单元注释、VNTR HMM/SKIP 分解和精简字典；保留独立的 D1S80 16 bp 单元注释库。

程序独立实现，不依赖 danbing-tk、TRGT、LongTR、vamos 或 PGGB。当前属于研究实现；图路径、局部支持和单元注释不能自动视为已验证的基因型。

新增完整观测序列的 **GFA S/L/P 结构图**：`scripts/vntrpathgraph_cli.py`，输出共享节点、分支、跳跃边与静态 SVG/PNG，默认无损；可通过显式逐位替换规则生成带损失台账的结构图。详见[设计、运行与限制](docs/PATH_SEQUENCE_GRAPH.zh-CN.md)。本功能分支为 `feat/path-structure-graph-20261105`。

新增 **结构优先模式** `--graph-mode structure`：接入冻结单元划分，以重复段/模块建图，局部保留划分歧义，支持每条 path 累计 0–2 bp 的可逆显示归并。见[使用与规则](docs/STRUCTURE_FIRST.zh-CN.md)和[ATXN2/D1S80 计算节点验证](docs/STRUCTURE_FIRST_PILOT_20261106.zh-CN.md)。分支：`feat/structure-first-20261106`。

## Versions and branches

**当前分支为 `codex/structure-supernodes-20261118`**。当前功能和删除后的兼容性变化见 [使用与核验说明](docs/debruijn-removal-audit.md)。下表是历史开发分支。

| 分支 | 用途 |
|---|---|
| [review-fixes](https://github.com/devildoctor/experiment-trs-map/tree/fix/vntr-review-20261029) | 当前修复、回归测试和中文代码解读 |
| [vntr-point-graph](https://github.com/devildoctor/experiment-trs-map/tree/vntr-point-graph) | 基础点图、STR图及多k值VNTR图 |
| [joint-support-unit-decomposition](https://github.com/devildoctor/experiment-trs-map/tree/feature/d1s80-joint-support-unit-decomposition) | 初版联合支持与候选单元分解 |
| [joint-units-complete](https://github.com/devildoctor/experiment-trs-map/tree/feature/d1s80-joint-units-complete) | 完整事件证据及所有路径的单元注释 |
| [terminal-unit-evidence](https://github.com/devildoctor/experiment-trs-map/tree/feature/d1s80-terminal-unit-evidence) | 经分子验证的右端模板，保留修复前版本 |

各分支保存独立开发阶段；旧分支不具备后续分支的全部功能。仓库只提供程序、测试和目录配置，不提供真实 BAM、参考 FASTA 或逐分子分析结果。

## Installation on Linux

需要 Python 3.11+、samtools，以及图布局使用的 Graphviz (`dot` / `neato`)。锚定模块依赖固定版本的 edlib 和 pyahocorasick。

```bash
git clone --branch codex/structure-supernodes-20261118 https://github.com/devildoctor/experiment-trs-map.git vntr-point-graph
cd vntr-point-graph
conda create -n vntr-point-graph -c conda-forge -c bioconda python=3.13 samtools graphviz pip
conda activate vntr-point-graph
python -m pip install -r requirements-anchors.txt -r requirements-flanks.txt
python scripts/vntrgraph.py --help
```

samtools 和参考构建版本应记录在运行记录中。集群脚本中的账号、分区、Python和数据路径是原部署配置，换环境后必须先调整。

## 从真实输入到图

以下命令从仓库根目录执行。将环境变量 `REF`、`BAMS`、`MEMBERS` 指向带 FAI 的 GRCh38 FASTA、按位点整理的 BAM 目录及已有成员 JSON。输出目录应为新目录。

```bash
CAT=catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed
python scripts/vntrpathgraph_cli.py --bam-dir "$BAMS" --reference "$REF" --catalog "$CAT" --locus D1S80 --graph-mode structure --locus-model catalogs/D1S80.structure.v1.json --supernode-policy catalogs/D1S80.supernodes.v1.json --q 0.2 --out out_bam
python scripts/run_efficient_family.py --input "$MEMBERS" --reference "$REF" --model catalogs/D1S80.structure.v1.json --supernode-policy catalogs/D1S80.supernodes.v1.json --q 0.2 --out out_members
python scripts/vntrsupernode_export.py out_bam/D1S80 --source-root .
python scripts/vntrpathgraph_cli.py --bam-dir "$BAMS" --reference "$REF" --catalog "$CAT" --locus ATXN2 --graph-mode structure --locus-model catalogs/ATXN2.structure.v1.json --structure-policy catalogs/ATXN2.structure-policy.v1.json --out out_str
python scripts/vntrgraph.py run --bam-dir "$BAMS" --reference "$REF" --catalog "$CAT" --locus D1S80 --out out_candidates
```

D1S80 打开 `D1S80.supernodes.html`。ATXN2 的 `ATXN2.observed-str.html`、同名 GFA 与 SVG 使用相同的观察路径。候选入口的候选及其支持数另有统计口径。省略超节点的 `--q` 使用全字典。B2 v2.2 在切点审计及汇总中区分完整空核心和空/非空并列。全切点重合的观察按分子键登记，GFA 的观察 P 保存非空原串。

超节点包中的 `*.path_sequences.txt` 按 GFA 路径名列出原串。成员导入记录可独立核实的位点字段；BAM 输入清单登记每份 BAM、BAI 的 SHA-256、参考核心、边界模型及目录。来源文件哈希可在对应 Git 检出上复核。历史输入与旧交付包按其记录的版本读取。

## Input preparation

| 输入 | 格式与要求 |
|---|---|
| 位点读段 | BAM及匹配BAI；与参考使用相同染色体命名和组装版本 |
| 参考 | FASTA及FAI，保留实际采用的ALT/decoy信息 |
| 位点目录 | BED，0-based半开坐标；明确搜索窗口与重复区边界 |
| 可选分子来源 | TSV列名 `bam, namespace`；相对路径相对此TSV解析 |
| 可选参考等价关系 | 明确验证的参考区间等价表，不能用来默认忽略ALT |

推荐目录：

```text
bamdata_by_locus/<locus>/<sample>/<HiFi|ONT|ONTUL>/<name>.<locus>.bam
```

BED第4列使用 `ID=...;MOTIFS=...;REPEAT_START=...;REPEAT_END=...;BOUNDARY=confirmed`。
第2、3列是搜索窗口，重复区必须位于窗口内。`confirmed` 是输入目录的声明，不是程序新提供的实验验证。不能为凑整数单元而裁剪端部。

主图、read-view、structure/summary 的分子键统一为 `[sample, namespace, read_name]`。`namespace` 优先取 `--sources` 的显式配置，否则取该记录 RG 对应的 LB，缺失时为空；旧结果缺 namespace 时也按空值兼容。BAM 路径和平台不自动成为 namespace。同一分子跨BAM、跨平台重复出现只计一次；平台小计可能重叠。主/补充记录不会拼接成未观测的完整序列。

## Usage

以下命令应在已分配的计算节点执行；每次使用新的输出目录。

### Scenario 1: anchored point graph

```bash
REFERENCE=/path/to/reference.fa
BAMS=/path/to/bamdata_by_locus
CATALOG=catalogs/grch38_pointgraph_10_loci_v1.bed
RUN=/path/to/new-run
python scripts/vntrgraph.py run --bam-dir "$BAMS" --reference "$REFERENCE" \
  --catalog "$CATALOG" --locus D1S80 --samtools "$(command -v samtools)" \
  --out "$RUN/pointgraph"
```

D1S80严格模式可能因主染色体和ALT上的锚点重复而得不到完整候选。只有验证参考与等价表匹配后，才显式增加
`--reference-equivalence catalogs/grch38_d1s80_alt_equivalence_v1.tsv`；结果应继续保留等价映射假设。

## Output files

| 文件/目录 | 解释 |
|---|---|
| `<locus>.gfa/.svg/.json` | 点图及明确保存的路径 |
| `molecules.tsv`, `anchor_decisions.jsonl` | 原始分子切片、来源、锚定判断与拒绝原因 |
| `candidates.tsv`, `assignments.tsv`, `genotypes.tsv` | 候选、分配、候选组合与不确定性 |
| `motifs.tsv`, `tokens.tsv` | 字典与原始切片注释 |
| `manifest.json`, `parameters.tsv` | 方法版本、参数、输入来源与哈希 |

原始碱基和分子身份可能存在于TSV/JSON/GFA中；这些分析输出不是源码发布内容。

点图切片之间使用0M连接。GFA节点、边数或自环次数不能直接当作生物学CN。

## Tests and reproducibility

```bash
python -m unittest discover -s tests -v
python tests/integration_synthetic.py
```

真实数据运行应保存源码提交、输入哈希、参考身份、参数、调度日志和独立验收结果。Slurm模板位于 `slurm/`，要求先准备其指定的源码快照及校验清单，不能直接在登录节点运行分析。

## Interpretation and limitations

- REF是参考上下文，OBS是实际观测，CAND是推断候选；REF/CAND不贡献观测分子票。
- 局部片段可以支持所覆盖区间，不能据此证明整条候选或完整CN。
- 完整单元、partial与unknown分别报告。长度除以16仅是描述性长度比。
- 描述性编辑代价、条件支持比例均不是基因型概率或等位频率。
- 注释歧义、倍性冲突和位置等价假设必须保留；不能强取两个最大簇作为确定基因型。
- 当前测试证明实现与产物契约，不构成人群分型准确率或临床性能验证。

## Documentation

- [2026-10-29 历史修复与验证](docs/CODE_FIXES_20261029.zh-CN.md)
- [历史中文代码解读](docs/CODE_WALKTHROUGH.zh-CN.md)、[原始检查记录](docs/CODE_REVIEW_20261029.zh-CN.md)

- [完整方法与历史说明](METHODS.md)
- [锚定与来源规则](ANCHOR_UPDATE.md)
- [STR表达](STR_EXPRESSION.md)与[STR单元图](STR_GRAPH.md)
- [读段视图](READ_VIEW.md)、[结构视图](STRUCTURE_VIEW.md)

README组织方式参考 [ChaissonLab/danbing-tk](https://github.com/ChaissonLab/danbing-tk) 的功能、依赖、分场景运行和文件说明。算法与实现均按本仓库实际代码说明，不代表复现danbing-tk的性能。问题反馈请附源码提交、命令、参数和已去除敏感数据的错误日志，提交至仓库Issues。
