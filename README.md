# VNTR 点位图生成程序

本分支独立实现 Word 方案中的 VNTR 点位图流程。输入为 `bamdata/样本/平台/*.位点.bam`、匹配的 `reference/*.fa` 与 `.fai`、以及 `catalogs/*.bed`。程序不会调用旧的重复区域提取程序。内部坐标采用 0-based、半开区间。

## 文件说明

| 文件 | 作用 |
| --- | --- |
| `scripts/vntrgraph.py` | BAM/FASTA 读取、双侧翼定界、聚类共识、局部 k-mer 图检查、motif 分解、GFA/SVG 输出与命令行 |
| `catalogs/VNTR_candidates.bed` | 8 个目标位点的 GRCh38 候选窗口；无预设 motif，须从序列发现 |
| `tests/test_vntrgraph.py` | 坐标、定界、CIGAR、等位路径可逆性及未知序列回归测试 |
| `slurm/vntrgraph.sbatch` | 在创新港计算节点上运行单个位点或全部位点 |
| `.gitignore` | 排除大体积 BAM、参考序列和运行产物 |

## 方法与判定

1. 核对 BED、FAI 与 BAM 的染色体长度，扫描全套参考 FASTA，要求左右侧翼各至少一个 31-mer 在双链上唯一。
2. 按位点读取 BAM 主比对和补充比对，用 CIGAR 只限定搜索范围；在 read 内重新匹配左右参考侧翼，拒绝缺失、歧义或逆序锚点。同一 read 只计一次。
3. 按样本和平台分别聚类。每簇至少 3 条双锚 read，以多数证据形成共识；只有全长序列高度一致才跨平台合并。读段数代表观测支持，不表示染色体拷贝数。
4. 建立 read 真实邻接形成的局部 k-mer 图，检查共识是否有读段支持。再从参考和共识寻找候选周期，用动态规划得到精确碱基 token；未解释区段标记为 `UNKNOWN`。
5. GFA 的每个节点存精确碱基，路径重复列出重复节点，能逐碱基还原完整序列。SVG 将连续相同 token 压缩为 `×N`，并显示平台支持。低证据、未知区段或复杂多等位标为待定，不报告整数拷贝数。

BED 边界是候选搜索窗口，并不自动成为生物学上的精确重复边界。对于复杂区域，须检查 `audit.tsv`、`bam_audit.tsv` 和原始 BAM 后再解释图。无可靠双锚 read 时只输出说明原因的 SVG。

## 环境

- Linux Python 3.11 或更新；无需 Python 第三方包。服务器的系统 `python3` 为 3.6，应使用 `/home/data/user/stzhu/miniconda3/bin/python3`。
- `samtools` 1.20 或更新；BAM 必须按坐标排序且有 BAI。
- 与 BAM 比对所用版本一致的 GRCh38 FASTA 和 `.fai`。
- 全基因组锚点唯一性扫描会读取完整 FASTA，建议在计算节点运行。

## 运行

在仓库根目录执行：

```bash
/home/data/user/stzhu/miniconda3/bin/python3 scripts/vntrgraph.py run \
  --bam-dir bamdata \
  --reference reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --catalog catalogs/VNTR_candidates.bed \
  --out results/vntrgraph \
  --samtools /home/data/user/stzhu/miniconda3/envs/samtools_env/bin/samtools \
  --locus ABCA7
```

去掉 `--locus ABCA7` 即处理 BED 中全部 8 个位点。结果在 `results/vntrgraph/位点/`，包括 `位点.gfa`、`位点.svg`、`audit.tsv`、`bam_audit.tsv` 和 `parameters.tsv`。没有可信等位路径时不会生成 GFA。

## 详细测试步骤

1. 检查输入：`find bamdata -name '*.bam' | wc -l` 应显示 68；`test -f reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai`；运行 `samtools quickcheck -v bamdata/LCL8/HiFi/*ABCA7.bam`，应无错误输出。
2. 运行单元测试：`/home/data/user/stzhu/miniconda3/bin/python3 -m unittest discover -s tests -v`。必须全部 `ok`。
3. 运行合成 BAM 闭环测试：`/home/data/user/stzhu/miniconda3/bin/python3 tests/integration_synthetic.py`。应输出 `synthetic integration OK`。
4. 提交单个位点：`sbatch --export=ALL,LOCUS=ABCA7 slurm/vntrgraph.sbatch`。用 `squeue -u "$USER"` 查状态，作业完成后查看 `vntrgraph.<jobid>.out` 和 `.err`。
5. 核对 `results/vntrgraph/ABCA7/bam_audit.tsv` 的每个 BAM 入选 read 数、拒绝原因；检查 `audit.tsv` 中支持 read 身份、平台支持、回贴相似度、k-mer 图支持和待定标记。
6. 若存在 `ABCA7.gfa`，逐条读取 P 记录，按 S 节点序列连接。该序列必须等于左侧 25 bp + 对应共识重复序列 + 右侧 25 bp；程序构图时也执行同一断言。检查每条 L 边都在至少一条 P 路径中出现。
7. 浏览 `ABCA7.svg`，确认不同样本分行、精确 motif 变体分色、连续同序列单元显示 `×N`、待定路径不出现整数拷贝数。SVG 与 GFA 中等位名称和平台支持应相同。
8. 单个位点通过后，执行 `sbatch slurm/vntrgraph.sbatch` 处理 8 个位点。对不出 GFA 的位点，检查 SVG 错误说明，不手工填入猜测拷贝数。

## 当前判定边界

程序仅在跨越 read 充分且序列可解释时显示确定拷贝数。多数表决共识与 read 支持的局部图用于保守筛选；复杂重排、超长扩增和多个邻接子区可能保留为待定路径，需要进一步人工核查。本实现不把候选 BED 或读段覆盖峰视为确定边界。
