# STR 精确表达层（ATXN2 试点）

## 启用与范围

在独立 BED 目录中添加 `STR_STRUCTURE=(CTG)*`。试点目录为
`catalogs/grch38_str_expression_atxn2_v1.bed`，坐标保持
`chr12:[111598950,111599019)`，搜索区间保持 `[111598900,111599070)`。
原有三个 BED 和所有 >6 bp VNTR 点位配置保持不变。

仅接受 `(2–6 bp ACGT 单元)*` 与固定 ACGT 字面序列，例如
`(CAG)*TTGACC(CTG)*`。`STRUC=<TR>` 仍保持原义，不自动启用。
没有 STR_STRUCTURE 时不产生 STR 文件，不改变既有注释、分组和推断。
显式启用但目录含 1 bp、>6 bp 或缺少 motif 时拒绝配置，避免把 VNTR
内部偶然出现的短重复当成整个位点的 STR。配置语法错误会终止运行。

模型借鉴 ExpansionHunter 的结构描述形式。本模块只做既有序列的精确
表示，没有运行 ExpansionHunter，也没有实现其分型或置信区间估计。

## 计数规则

- 用配置中的 motif 原始相位进行精确匹配，不接受近似编辑。
- motif 若本身是精确复合单元（例如 CAGCAG）则以 CAG 计数，同时保存
  catalog_motif 与 repeat_factor。不会先缩短 >6 bp 目录单元再纳入试点。
- 连续至少两次才形成重复块；孤立一个单元、N、错误、中断、端部残片
  均保留为字面碱基。例如 `(CAG)5–CAA–(CAG)7–CA`。
- 不同位置的重复块分别计数。按字面碱基数、重复块数依次最小化；
  相同代价时按坐标和 motif 排序，保留一个备选并标记同分歧义。
  不是穷举所有解释，代价不是概率。
- 模型只提供允许压缩的 motif 和吻合检查；不强迫实际片段符合模型。
  `model_match=false` 可以表示真实中断、片段截断或测序误差。
- 表达式能重建完整输入序列。候选序列的次数只描述候选；
  observed_fragment 的次数只描述实际截取片段；参考单独标记。
- 存储方向是参考正链。反向展示时整个块序列反向互补，固定中断也一起
  反向互补。不会把参考 CTG 直接改名为 CAG 而忽略顺序。

## 输出与图形

三个模式共用 scripts/vntrstr.py，新增 str_paths.tsv 和 str_paths.json；
read/structure 模式保留 str-svg/ 分页 SVG。JSON 保存准确碱基、片段内起止、单位、次数、备选、
模型吻合、序列 SHA256 和来源。TSV 便于查询。

主图 ATXN2.svg 显示 STR 单元节点与边，ATXN2.paths.svg 按路径显示次数。
主 GFA/JSON 同步为单元图；原 GFA/JSON/SVG 成套另存 ATXN2.detail.*。
完整命名和图语义见 STR_GRAPH.md。候选表、分子分配和阈值均保留。
read HTML/SVG 用实际 motif 与字面碱基展示，完整节点仍在 HTML 详情。
structure HTML 的参考、代表路径、亚型和读段使用相同表达层；
原 M 单元注释可展开。类和模块仍按原方法生成。
蓝色框是精确重复，金色框是字面碱基；框宽不代表 bp 长度。

部分覆盖、锚点、MAPQ、读段 flag 和候选 uncertainty 均保留。
本次表达改进不解决 ATXN2 的 ploidy_conflict，也不自动挑出两个等位。

## Slurm

按 README 的约定创建包含 scripts/tests/catalogs/slurm 的版本快照并生成
SHA256SUMS。设置 RUN_DIR 后提交 `slurm/vntrstr.sbatch`，仅运行 ATXN2，
输出为 `results/{pointgraph,readview,structure}/ATXN2/`。
测试、参考读取、BAM 解码和图形渲染均应在计算节点进行。

参考：
- https://github.com/Illumina/ExpansionHunter/blob/master/docs/04_VariantCatalogFiles.md
- 官方 hg38 ATXN2 为 `(GCT)*`、`chr12:[111598949,111599018)`，与本项目
  相差一个碱基。试点不自动改边界；表达保持本项目参考相位并核验重构。
