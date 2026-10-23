# vntr-point-graph

面向长读段 VNTR / STR 的可审计点位图分析工具。输入按位点整理的 BAM、匹配参考及位点目录，输出明确保存的序列路径、分子证据、重复单元注释和 GFA/SVG。

主要功能：

1. 依据两侧锚点提取实际读段序列，按独立分子去重，保留定位与候选冲突。
2. 生成候选序列点图，并提供可选 STR 单元图与 VNTR compact de Bruijn 图。
3. 在冻结的 D1S80 图输入上分析连续跨分支支持与16 bp基本单元分解。

程序独立实现，不依赖 danbing-tk、TRGT、LongTR、vamos 或 PGGB。当前属于研究实现；图路径、局部支持和单元注释不能自动视为已验证的基因型。

## Versions and branches

**当前分支：`feature/d1s80-joint-units-complete`** — 完整证据版：显式支持/冲突/无信息状态，覆盖全部参考、观测与候选路径。

| 分支 | 用途 |
|---|---|
| [vntr-point-graph](https://github.com/devildoctor/experiment-trs-map/tree/vntr-point-graph) | 基础点图、STR图及多k值VNTR图 |
| [joint-support-unit-decomposition](https://github.com/devildoctor/experiment-trs-map/tree/feature/d1s80-joint-support-unit-decomposition) | 初版联合支持与候选单元分解 |
| [joint-units-complete](https://github.com/devildoctor/experiment-trs-map/tree/feature/d1s80-joint-units-complete) | 完整事件证据及所有路径的单元注释 |
| [terminal-unit-evidence](https://github.com/devildoctor/experiment-trs-map/tree/feature/d1s80-terminal-unit-evidence) | 经分子验证的右端模板，建议从此分支开始 |

各分支保存独立开发阶段；旧分支不具备后续分支的全部功能。仓库只提供程序、测试和目录配置，不提供真实 BAM、参考 FASTA 或逐分子分析结果。

## Installation on Linux

需要 Python 3.11+、samtools，以及图布局使用的 Graphviz (`dot` / `neato`)。锚定模块依赖固定版本的 edlib 和 pyahocorasick。

```bash
git clone --branch feature/d1s80-joint-units-complete https://github.com/devildoctor/experiment-trs-map.git vntr-point-graph
cd vntr-point-graph
conda create -n vntr-point-graph -c conda-forge -c bioconda python=3.11 samtools graphviz pip
conda activate vntr-point-graph
python -m pip install -r requirements-anchors.txt
python scripts/vntrgraph.py --help
python scripts/vntrdbg_analysis.py --help
```

samtools 和参考构建版本应记录在运行记录中。集群脚本中的账号、分区、Python和数据路径是原部署配置，换环境后必须先调整。

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

分子键为 `[sample, namespace, read_name]`。同一分子跨BAM、跨平台重复出现只计一次；平台小计可能重叠。主/补充记录不会拼接成未观测的完整序列。

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

### Scenario 2: compact VNTR graph and k comparison

先完成同一批输入的点图提取，再运行：

```bash
python scripts/vntrgraph.py --dbg-view --evidence-dir "$RUN/pointgraph/D1S80" \
  --bam-dir "$BAMS" --reference "$REFERENCE" --catalog "$CATALOG" \
  --models catalogs/grch38_graph_models_v1.tsv --locus D1S80 \
  --k 17 21 31 33 49 65 --samtools "$(command -v samtools)" --out "$RUN/dbg/D1S80"
python tests/validate_dbg.py "$RUN/dbg/D1S80"
```

边来自真实连续k+1-mer。支持视图要求同一样本至少3个独立分子；被过滤的路径拆段，不补桥。k是图上下文长度，不是重复单元长度。完整语义和输出见 [VNTR_DBG.md](VNTR_DBG.md)。

### Scenario 3: joint branch evidence and repeat units

```bash
python scripts/vntrdbg_analysis.py --dbg-dir "$RUN/dbg/D1S80" \
  --pointgraph-dir "$RUN/pointgraph/D1S80" --out "$RUN/joint-units" --k 21 33
```

仅接受可核对来源哈希的冻结D1S80输入。连续跨分支支持要求一个观测路径覆盖完整有序区间，不能把不同片段的票相加当作联合支持。完整候选匹配要求双侧锚定及整条序列精确一致。详见 [VNTR_DBG_ANALYSIS.md](VNTR_DBG_ANALYSIS.md)。



## Output files

| 文件/目录 | 解释 |
|---|---|
| `<locus>.gfa/.svg/.json` | 点图及明确保存的路径 |
| `molecules.tsv`, `anchor_decisions.jsonl` | 原始分子切片、来源、锚定判断与拒绝原因 |
| `candidates.tsv`, `assignments.tsv`, `genotypes.tsv` | 候选、分配、候选组合与不确定性 |
| `motifs.tsv`, `tokens.tsv` | 字典与原始切片注释 |
| `input_paths.json`, `input_summary.json` | dBG冻结输入和校验信息 |
| `k*/` | 各k的完整证据图、支持图及逐路径SVG |
| `joint_queries.tsv`, `joint_support.tsv`, `joint_molecules.tsv` | 有序区间、去重支持与逐分子证据 |
| `primitive_units.tsv`, `unit_summary.tsv` | 基本单元的主/备选分解及完整/部分/未知计数 |
| `manifest.json`, `parameters.tsv` | 方法版本、参数、输入来源与哈希 |

完整证据版另有 `events.tsv`, `molecule_events.tsv`, `motif_counts.tsv`；初版不提供全部这些文件。
原始碱基和分子身份可能存在于TSV/JSON/GFA中；这些分析输出不是源码发布内容。

点图切片之间使用0M连接；compact dBG使用k−1重叠，重建路径时必须去除重叠。GFA节点、边数或自环次数不能直接当作生物学CN。

## Tests and reproducibility

```bash
python -m unittest discover -s tests -v
python tests/integration_synthetic.py
```

真实数据运行应保存源码提交、输入哈希、参考身份、参数、调度日志和独立验收结果。Slurm模板位于 `slurm/`，要求先准备其指定的源码快照及校验清单，不能直接在登录节点运行分析。

此分支保留对应开发阶段的测试与验证入口。终端模板的186项测试和最新D1S80结果属于后续终端功能分支，不应归属于本历史版本。

## Interpretation and limitations

- REF是参考上下文，OBS是实际观测，CAND是推断候选；REF/CAND不贡献观测分子票。
- 局部片段可以支持所覆盖区间，不能据此证明整条候选或完整CN。
- 完整单元、partial与unknown分别报告。长度除以16仅是描述性长度比。
- 描述性编辑代价、条件支持比例均不是基因型概率或等位频率。
- 注释歧义、倍性冲突和位置等价假设必须保留；不能强取两个最大簇作为确定基因型。
- 当前测试证明实现与产物契约，不构成人群分型准确率或临床性能验证。

## Documentation

- [完整方法与历史说明](METHODS.md)
- [锚定与来源规则](ANCHOR_UPDATE.md)
- [STR表达](STR_EXPRESSION.md)与[STR单元图](STR_GRAPH.md)
- [读段视图](READ_VIEW.md)、[结构视图](STRUCTURE_VIEW.md)
- [VNTR图](VNTR_DBG.md)、[联合支持与分解](VNTR_DBG_ANALYSIS.md)

README组织方式参考 [ChaissonLab/danbing-tk](https://github.com/ChaissonLab/danbing-tk) 的功能、依赖、分场景运行和文件说明。算法与实现均按本仓库实际代码说明，不代表复现danbing-tk的性能。问题反馈请附源码提交、命令、参数和已去除敏感数据的错误日志，提交至仓库Issues。
