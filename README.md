# vntr-point-graph · VNTR / STR 点位图

在一个串联重复位点，哪些单元连续出现了多少次，中断在哪里，哪些分子实际走过同一条路径？根目录的 `vntr-point-graph` 从长读段 BAM 或已登记的路径成员建图，保留核心原串、单元顺序和分子支持。网页上的观察路径可以沿 GFA 走法逐碱基核对。

`main` 是仓库首页和当前运行入口。原有的候选读段筛选工具保留在 [experiment-trs-map/](experiment-trs-map/README.md)，可以独立使用。

## 选择入口

| 我想得到 | 从仓库根目录使用 | 运行后先看 |
|---|---|---|
| D1S80 的 VNTR 重复单元、路径和分子支持 | `scripts/vntrpathgraph_cli.py`，结构模式加超节点策略 | `D1S80/D1S80.supernodes.html`，默认对比五条观察路径 |
| ATXN2 的 STR 单元次数和中断 | 同一入口，使用 ATXN2 冻结模型和 STR 结构策略 | `ATXN2/index.html`，也可打开 `ATXN2.observed-str.svg` 或同名 GFA |
| 用已有成员重建 VNTR 图，比较全字典与精简字典 | `scripts/run_efficient_family.py` 加 `--supernode-policy` | 输出目录内的 `D1S80.supernodes.html`、`summary.json` 和改写账本 |
| 按锚点建立候选序列点图 | `scripts/vntrgraph.py run` | 位点目录中的 `candidates.tsv`、`assignments.tsv`、`anchor_decisions.jsonl` |
| 从单个 BAM 筛选候选 reads，输出 BAM、FASTQ 和审计表 | `experiment-trs-map/scripts/extract_candidate_reads.py` | 筛选工具的 `*.manifest.tsv`、`*.audit.tsv`。[用法与参数](experiment-trs-map/README.md) |

## 安装与输入

下列示例在 Linux 上从仓库根目录运行。项目提供锚点比对和 B2 数值计算的固定依赖文件，已使用 CPython 3.13 验证。samtools 用于读写 BAM 和建立索引，Graphviz 用于静态布局。离线 HTML 可直接用浏览器打开。

```bash
git clone https://github.com/devildoctor/experiment-trs-map.git
cd experiment-trs-map
conda create -n vntr-point-graph -c conda-forge -c bioconda python=3.13 samtools graphviz pip
conda activate vntr-point-graph
python -m pip install -r requirements-anchors.txt -r requirements-flanks.txt
python scripts/vntrpathgraph_cli.py --help
```

- **参考序列**使用与 BAM 比对时相同的 FASTA 和配套 `.fai`。示例目录针对 GRCh38。
- **位点 BAM** 按 `<BAMS>/<位点>/<样本>/<HiFi|ONT|ONTUL>/<前缀>.<位点>.bam` 存放，每份 BAM 需要匹配的 BAI。
- **位点目录**用 0 起点、右端不含的 BED。第 2、3 列定义搜索窗口，`REPEAT_START` 和 `REPEAT_END` 可指定其中的核心。比较位点时保留 `BOUNDARY` 的原始声明。[示例 BED](catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed)
- **分子身份**按样本、命名空间和读段名去重。默认从读组的文库标识取得命名空间，也可用 `--sources` 提供带 `bam`、`namespace` 两列的 TSV。跨平台小计可能重叠。

## 从 BAM 到观察路径图

先把 `REF` 与 `BAMS` 改为实际路径，并为 `OUT` 选择新的运行目录。两条命令使用同一份目录表，各自写入位点子目录。

```bash
REF=/path/to/GRCh38.fa
BAMS=/path/to/bamdata_by_locus
CAT=catalogs/grch38_structure_summary_atxn2_d1s80_abca7_v1.bed
OUT=run-point-graphs

python scripts/vntrpathgraph_cli.py --bam-dir "$BAMS" --reference "$REF" \
  --catalog "$CAT" --locus D1S80 --graph-mode structure \
  --locus-model catalogs/D1S80.structure.v1.json \
  --supernode-policy catalogs/D1S80.supernodes.v1.json --q 0.2 --out "$OUT"

python scripts/vntrpathgraph_cli.py --bam-dir "$BAMS" --reference "$REF" \
  --catalog "$CAT" --locus ATXN2 --graph-mode structure \
  --locus-model catalogs/ATXN2.structure.v1.json \
  --structure-policy catalogs/ATXN2.structure-policy.v1.json --out "$OUT"
```

D1S80 示例以 HMM 注释并选择有效字典，省略 `--q` 即使用全字典。`q` 控制字典选择，随附超节点策略中的 3% 是每条观察路径允许的表示层碱基改写上限。ATXN2 示例沿用冻结字典和 `seed_dp`。图上分别标出该路径的连续重复次数和支持它的分子键数。

### 从同一批成员重新建图

下面直接使用上一步生成的成员文件。成员入口的 `--out` 就是最终包目录，可对两个入口分别做完整性检查。

```bash
python scripts/run_efficient_family.py --input "$OUT/D1S80/source_members.json" \
  --reference "$REF" --model catalogs/D1S80.structure.v1.json \
  --supernode-policy catalogs/D1S80.supernodes.v1.json --q 0.2 \
  --out "$OUT/D1S80-from-members"

python scripts/vntrsupernode_export.py "$OUT/D1S80" --source-root .
python scripts/vntrsupernode_export.py "$OUT/D1S80-from-members" --source-root .
```

成员与模型有明确冲突时，程序在创建结果目录前报错。旧成员缺少位点名、组装等来源字段时，`input_identity.json` 逐项列出能核实的证据和缺项。

### 生成候选序列或筛选读段

已有按位点整理的 BAM 时，第一条命令按锚点建立 D1S80 候选点图。只有单个已排序并建索引的 BAM、需要先选出读段时，可以使用第二条命令，示例目录指定 DMPK。

```bash
python scripts/vntrgraph.py run --bam-dir "$BAMS" --reference "$REF" \
  --catalog "$CAT" --locus D1S80 --out "$OUT/candidates"

python -B experiment-trs-map/scripts/extract_candidate_reads.py \
  --bam /path/to/sample.bam --repeats experiment-trs-map/catalogs/example.repeats.bed \
  --repeat-id DMPK --output-dir "$OUT/reads-DMPK"
```

## 从结果追溯到碱基

| 输出 | 能核对什么 |
|---|---|
| `<位点>.path_sequences.txt`、`source_members.json` | 按 GFA 的 P 名查角色、原始长度、完整核心序列和分子成员。P0001 是参考，观察路径从 P0002 起。 |
| `extraction_audit.jsonl`、`boundary_model.json` | BAM 入口的逐记录切点、并列和使用的边界模型。B2 v2.2 对完整空核心及空/非空并列单列证据。 |
| `<位点>.original.gfa`、`<位点>.represented.gfa`、`losses.jsonl` | VNTR 超节点包的原始走法、表示走法和逐碱基改写。按账本可以反向还原原串。 |
| `path_visits.jsonl`、`display_graph.json` | VNTR 结构块的有序访问、次数、位置可比性和对应的原始 GFA 切片。 |
| `input_manifest.json`、`provenance.json`、`validation.json` | VNTR 超节点包的输入、模型与源码版本、文件哈希以及重建检查。 |
| `ATXN2.observed-str.gfa`、`ATXN2.observed-str.svg` | STR 观察核心的同一组登记走法，以连续次数和冻结中断单元展示。 |

同一条 P 可由多个分子键支持。两条相接边各自有支持时，仍须查看有多少条登记 P 连续经过两边。节点数、自环次数和描述性长度比均需结合具体路径及边界证据解读。报告生物学拷贝数或基因型还需独立验证。

## 测试与其他工作流

```bash
export SAMTOOLS="$(command -v samtools)"
python -m unittest discover -s tests
python tests/integration_synthetic.py
python -B -m unittest discover -s experiment-trs-map/tests
```

前两项分别检查点位图的接口和合成 BAM 到图的流程，第三项检查旧筛选工具。前端的投影与序列比较可用 Node.js 运行 `node tests/test_supernode_projection.js` 和 `node tests/test_supernode_compare.js`。

- [B2 核心提取](docs/flank-b2.md) · [完整观察路径的定义与 GFA](docs/PATH_SEQUENCE_GRAPH.zh-CN.md) · [结构优先的分解规则](docs/STRUCTURE_FIRST.zh-CN.md)
- [锚点和来源](ANCHOR_UPDATE.md) · [读段视图](READ_VIEW.md) · [STR 单元图](STR_GRAPH.md) · [结构视图](STRUCTURE_VIEW.md)
- [VNTR 有效字典与家族图](docs/efficient-motif-family.md) · [HMM 参数化](docs/hmm-parameter-display.md) · [Slurm 作业模板](slurm/)。集群运行前请按本地环境调整分区、数据和依赖路径。

## 分支与历史版本

下表列出 2026-11-22 合入 `main` 的 GitHub 分支快照。链接固定到当时的提交，便于复现旧接口。现行根目录程序以经过验收的 `77d6bfd` 为基础。多 k 值 de Bruijn 图与联合分支证据可在相应历史版本中查看。[接口演变说明](docs/debruijn-removal-audit.md)

| 开发阶段 | 分支快照 |
|---|---|
| 原 `main` 的候选读段工具 | [main · 70dbcdb](https://github.com/devildoctor/experiment-trs-map/tree/70dbcdb451a8e075ac298a7c6565d3c03d05fd02)，现位于 [experiment-trs-map/](experiment-trs-map/README.md) |
| 点图、D1S80 联合证据与末端模板 | [vntr-point-graph](https://github.com/devildoctor/experiment-trs-map/tree/47d50079c255656a00f1e0caa5cb795e74b75c7b) · [feature/d1s80-joint-support-unit-decomposition](https://github.com/devildoctor/experiment-trs-map/tree/bd9542c1e213063dcfda31f2efc34dd5761d828e) · [feature/d1s80-joint-units-complete](https://github.com/devildoctor/experiment-trs-map/tree/5be0d5cfe9b707c1a64aac8e4e520ff27094329f) · [feature/d1s80-terminal-unit-evidence](https://github.com/devildoctor/experiment-trs-map/tree/5be507688fe7d1f8da4c1786fd85d6df8dfccb0e) |
| 证据语义、冻结字典与结构图 | [fix/vntr-review-20261029](https://github.com/devildoctor/experiment-trs-map/tree/bb3bb8eb0731283b4b3b99ba6d129ae6ae6bff2a) · [feat/frozen-partition-20261105](https://github.com/devildoctor/experiment-trs-map/tree/491ae44f732498a99f30f6d18b8ab5aeec455210) · [feat/path-structure-graph-20261105](https://github.com/devildoctor/experiment-trs-map/tree/e9bcadec55bd5b57a99e1611e3fbeddc2d7cd957) · [feat/structure-first-20261106](https://github.com/devildoctor/experiment-trs-map/tree/9c48c2b85a6516e6f6350c8dc464662f3fec1c38) |
| HMM、有效字典与观察路径超节点 | [feat/hmm-parameter-display-20261110](https://github.com/devildoctor/experiment-trs-map/tree/69df097fa485e438d21a52bd15a0c14c6952e09b) · [feat/efficient-motif-family-3pct-20261111](https://github.com/devildoctor/experiment-trs-map/tree/cd3adc0f7a6fb6b95345a5b19b7bd1745d9c3f18) · [codex/structure-supernodes-20261118](https://github.com/devildoctor/experiment-trs-map/tree/77d6bfdd7d374cfc5fbfc83274667d5c1ef8665e) |
