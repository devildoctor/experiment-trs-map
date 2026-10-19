# 全部原始 read 节点图 3.0

入口为 `scripts/vntrgraph.py --read-view`，实现于 `scripts/vntrreads.py`。
原候选共识入口仍可使用。此次指定任务仅运行 IL1RN、PER3、ABCA7。

每个 BAM 内的每个 QNAME 均保留一行。无 MAPQ、重复标记、QC 标记、最小支持数或双锚门槛。
单侧锚点允许定位；锚点均不可靠时，用 CIGAR 与位点区间相交得到实际观测片段。
左右各 500 bp 是用户提供 BAM 的选取上下文，不意味着 BAM 中的 SEQ 已裁剪到窗口；因此按 CIGAR 和 read 坐标提取，不能对 SEQ 简单去掉首尾 500 bp。
不进行全基因组唯一性前置拦截。双侧局部匹配仅作为截取证据，不声称全基因组唯一。

同一 QNAME 的主、补充和次级比对在同一行中分别保留为片段，不拼接重叠 SEQ。
跨 BAM 重复 QNAME 仍展示在各自文件下，不将行数解释为独立分子支持。
负链 BAM 存储 SEQ 与本记录 CIGAR 配套，不再次反向互补。
保留重复区内部以及边界上的插入；软剪切未知部分不推断为 VNTR；硬剪切缺失碱基明确不可用。
仅覆盖侧翼、SEQ 缺失等 read 仍有状态行，完整比对及存储 SEQ 写入 records.jsonl。

先复用已验证 motif 字典作为注释模板，没有字典时从参考发现。每条 read 单独分解实际碱基。
注释使用短种子限定候选位置的动态规划；种子未命中和不能解释的序列保留为 UNKNOWN。
HiFi/ONT 的编辑上限为单元长度的 12%/22%，限制 1 至 10；这是探索性参数，未校准为概率。
原始 read 错误不纠正成共识，不把观测单元数称作确定基因型或完整等位拷贝数。
每个节点为真实非重叠切片；GFA、SVG、HTML 统一节点 ID，并验证节点拼接恢复提取的序列。

每个位点输出：

- `位点.reads.html`：包含全部 read 数据的本地浏览页面，按样本、平台、read ID 筛选，分页只控制显示，无抽样。
- `svg/index.html`：全部静态 SVG 分页索引，每页 20 条 read；长节点路径自动换行，无节点截断。
- `位点.reads.gfa`：所有有效比对片段的节点与路径，0M 连接，无虚构侧翼节点。
- `位点.reads.json`、`reads.tsv`、`read_nodes.tsv`：全量路径、read 清单、真实碱基和节点经过记录。
- `records.jsonl`：全部原始比对信息、BAM 存储 SEQ 和质量，包括无法截取的记录。
- `status.json`：逐 BAM 记录数、read 行数、覆盖与锚点类别、重构检查计数。
- `parameters.json`、`motifs.tsv`：参数、源码哈希与 motif 来源。

每次运行使用新快照目录，旧结果不覆盖。Slurm 调度脚本 `slurm/submit_reads.sh` 固定提交三个点位。
