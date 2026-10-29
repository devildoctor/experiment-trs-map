# VNTR 审查问题修复与验证

更新日期：2026-10-30。修复以 `5be507688fe7d1f8da4c1786fd85d6df8dfccb0e` 及其已添加的中文注释为基础。GitHub 分支：[fix/vntr-review-20261029](https://github.com/devildoctor/experiment-trs-map/tree/fix/vntr-review-20261029)。

本次解决原检查报告的 R1–R6，并保留中文注释和详细功能解读。改动集中于输出一致性、输入路径、分子身份及阈值呈现。此前讨论的有损代表单元替换仍属于待设计功能；本次近似注释继续保存实际碱基。

## 修复内容

| 编号 | 触发条件 | 修复后行为 | 主要入口 |
|---|---|---|---|
| R1 | 两个同分分解方案，备选 token 原来与自身比较 | 主、备选各与另一方案比较；共同 token 不标歧义，非同分备选不误标 | `vntrcore.decompose`、`vntrstructure.segment` |
| R2 | BED 的 ID 为 `.`/`..`，或位点目录通过符号链接指向根外 | 解析时拒绝特殊 ID；主图/read/structure 写入前解析路径并要求严格位于输出根内 | `vntrgraph.validate_locus_id`、`locus_output_dir` |
| R3 | 全局有 OBS，但某个候选所属样本没有 OBS | 正常导出空分子表的固定表头；摘要 tested/support/conflict=0，fraction=NA，状态 `no_observed_molecules` | `vntrdbg_analysis.main`、`joint_support` |
| R4 | 两个独立文库复用 QNAME；或同一分子重复保存到不同 BAM/平台 | 全链路采用 `[sample, namespace, QNAME]`；独立文库分开计数，副本不重复投票 | `load_source_namespaces`、`read_bam`、`row_molecule_key`、结构与summary支持统计 |
| R5 | 建图或联合分支分析采用非默认门槛 | 颜色按实际阈值和样本内最大支持判断；报告读取参数；联合分析有独立 CLI 门槛，独立验收使用同一记录值 | `vntrdbg`、`vntrdbg_report`、`vntrdbg_analysis`及验收脚本 |
| R6 | 旧结构文档参数停留在早期版本 | 同步字典上限256、最多3轮、实际评分公式及不同支持口径；历史解读明确标注版本 | README、READ_VIEW、STRUCTURE_VIEW、VNTR_DBG、VNTR_DBG_ANALYSIS |

### 分子身份如何传递

1. `load_source_namespaces(path)` 统一解析 `bam/namespace` TSV。相对 BAM 路径以 TSV 所在目录为基准，同一规范路径重复出现直接报错。
2. 单个 BAM 有显式映射时采用该值，空字符串同样是有效覆盖；没有映射时取记录 RG 对应的 `@RG LB`，LB 缺失时使用空字符串。平台和 BAM 文件名不充当文库身份。
3. read-view 在每个 BAM 内按 `(namespace, QNAME)` 建行，主/补充/次级比对仍各自保留真实片段。显示行不会因跨 BAM 去重而删除。
4. structure 图的节点、边、模式、模板学习以及 summary 的类、模块、冲突判断都按三段分子键去重。
5. TSV/JSON 输出保存 namespace 与分子键；页面显示 namespace，区分 read 行数与独立分子数。原始片段顺序、实际碱基和覆盖状态保留。

旧行未保存 namespace 时按空字符串兼容。旧结果中已丢失的 LB 无法凭空恢复，需用原始 BAM 重跑才能区分历史同名文库。显示行 ID 现在包含 namespace，不能直接用新旧 ID 相等判断是否丢读段；历史验收改为比对原 BAM/QNAME 下的记录数与片段多重集合。

### 三种支持门槛分别记录

| 参数 | 所属阶段 | 默认值 | 记录位置 |
|---|---|---:|---|
| `--min-support` | compact dBG 支持视图 | 3 | 图 JSON、graph_semantics、metrics、k_comparison |
| `--branch-min-support` | 联合证据的分支定义 | 3 | `manifest.branch_min_support` 与 `branch_support_scope` |
| `--terminal-min-support` | 可选右端模板独立分子证据 | 3 | 原终端证据参数与审计 |

联合分支门槛不会自动继承支持图门槛。比较不同运行时应显式记录这三个阶段的实际配置。Slurm 联合分析脚本新增 `BRANCH_MIN_SUPPORT` 环境变量，默认3；终端历史验收脚本仍使用原固定默认值以对照旧结果。

图的 `support_scope` 和联合 `branch_support_scope` 为 `max_unique_molecules_within_one_sample`。`max_sample_support` 表示任一样本内的最大独立分子数；颜色据此判断是否达标。线宽及总分子标签仍描述跨样本去重总数，不能用总数替代样本内门槛。参考上下文可能在支持不足时仍被保留。

报告逐 k 读取实际参数。旧图若只在 semantics 文本中记录阈值，仍可读取；完全没有记录时才明确显示“推定历史默认值3”，不会把推定值冒充已记录参数。联合分析的独立验收同样标记旧默认值来源。

### 输出兼容性

- 新 TSV 身份与支持字段追加至表尾，原有列序号保持不变；仍建议按表头读取。
- 正常非空 `joint_molecules.tsv`、`molecule_events.tsv` 字段顺序保持不变。某候选样本没有 OBS 时新增状态；全局完全缺少 OBS 仍按既有输入规则拒绝。
- R1 只纠正备选 token 的歧义标记；主路径、真实序列、分解代价和备选方案选择保持。历史快照测试经过逐字段核对后更新，原哈希保留；先单独检查修正值，再仅规范化该字段核对旧摘要。
- 旧科学结果与冻结输入不覆盖。新运行使用新的输出目录，当前修复不会自动重算既有结果。

## 验证结果

- Slurm 作业 **167021**，计算节点 `cn01`，Python **3.13.9**，最终 `COMPLETED 0:0`。
- **221 项单元测试全部通过，无跳过**；包括 Linux 符号链接路径检查。
- 4项合成 BAM 集成测试通过：主点图、全量 read-view、structure，以及新增 namespace 测试（两个文库复用QNAME、跨BAM副本、相对来源路径覆盖、两种视图）。
- 16 bp VNTR 前后回归：17个结果文件逐字节一致；`parameters.tsv` 按设计记录不同源码哈希与运行路径，不要求字节一致。
- 同一组 R1/R2/R3/R4 最小探针确认：旧版能重现问题，修复版达到预期。R5 使用不同门槛的 CLI、图渲染参数、报告及独立验收回归。
- 测试前后源文件 SHA256 一致；HTML 模板 JavaScript 语法检查通过。
- 首次作业167020执行221例，唯一失败是历史分解快照包含旧错误歧义字段。逐字段核对后保留旧生物字段断言、补充修复字段断言，再在新快照运行完整验证。失败日志与修复说明保留，未将其计为通过。

验证限于实现契约和合成输入；没有重新分析真实 BAM，也没有据此声明分型准确率、疾病风险或 CN 已被验证。内存用量和长纯重复分解性能仍是已记录的扩展限制，本次没有做未经测量的性能重构。

### 复现命令

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
python tests/integration_synthetic.py
python tests/integration_all_reads.py "$SAMTOOLS"
python tests/integration_structure.py "$SAMTOOLS"
python tests/integration_review_namespaces.py "$SAMTOOLS"
python tests/regression_vntr_graph.py --baseline /path/to/preserved-source
```

依赖与计算节点使用方式见 [README](../README.md)。原始详细解读见 [CODE_WALKTHROUGH](CODE_WALKTHROUGH.zh-CN.md)，其中历史行号需以对应版本阅读；下面给出当前源码定义位置。

## 当前源码函数索引

本索引从当前 AST 生成，包含顶层函数、类及其内部定义。

### vntranchors.py

| 定义 | 当前行号 |
|---|---:|
| `dependencies` | [27](../scripts/vntranchors.py#L27) |
| `canonical` | [36](../scripts/vntranchors.py#L36) |
| `nonoverlapping` | [41](../scripts/vntranchors.py#L41) |
| `load_equivalences` | [50](../scripts/vntranchors.py#L50) |
| `normalized_position` | [77](../scripts/vntranchors.py#L77) |
| `reference_profiles` | [89](../scripts/vntranchors.py#L89) |
| `endpoint_costs` | [164](../scripts/vntranchors.py#L164) |
| `AnchorOverflow` | [194](../scripts/vntranchors.py#L194) |
| `flank_candidates` | [199](../scripts/vntranchors.py#L199) |
| `Record` | [254](../scripts/vntranchors.py#L254) |
| `Record.original_base` | [273](../scripts/vntranchors.py#L273) |
| `Record.original_interval` | [277](../scripts/vntranchors.py#L277) |
| `Record.original_sequence` | [281](../scripts/vntranchors.py#L281) |
| `Record.original_start` | [285](../scripts/vntranchors.py#L285) |
| `mapped_bases` | [290](../scripts/vntranchors.py#L290) |
| `parse_records` | [305](../scripts/vntranchors.py#L305) |
| `consistent_sequences` | [356](../scripts/vntranchors.py#L356) |
| `localization` | [374](../scripts/vntranchors.py#L374) |
| `competing_mapping` | [400](../scripts/vntranchors.py#L400) |
| `choose_pair` | [435](../scripts/vntranchors.py#L435) |
| `extract_records` | [454](../scripts/vntranchors.py#L454) |
| `extract_parsed` | [461](../scripts/vntranchors.py#L461) |
| `extract_bam` | [586](../scripts/vntranchors.py#L586) |
| `extract_bams` | [597](../scripts/vntranchors.py#L597) |

### vntrcore.py

| 定义 | 当前行号 |
|---|---:|
| `revcomp` | [23](../scripts/vntrcore.py#L23) |
| `compact` | [27](../scripts/vntrcore.py#L27) |
| `Read` | [33](../scripts/vntrcore.py#L33) |
| `Read.molecule_key` | [53](../scripts/vntrcore.py#L53) |
| `Read.observed_platforms` | [57](../scripts/vntrcore.py#L57) |
| `Token` | [63](../scripts/vntrcore.py#L63) |
| `Token.__iter__` | [76](../scripts/vntrcore.py#L76) |
| `Allele` | [82](../scripts/vntrcore.py#L82) |
| `deduplicate_reads` | [107](../scripts/vntrcore.py#L107) |
| `edit_distance` | [138](../scripts/vntrcore.py#L138) |
| `similarity` | [170](../scripts/vntrcore.py#L170) |
| `AlignmentUnresolved` | [174](../scripts/vntrcore.py#L174) |
| `align_to_seed` | [179](../scripts/vntrcore.py#L179) |
| `alignment_counts` | [227](../scripts/vntrcore.py#L227) |
| `make_consensus` | [247](../scripts/vntrcore.py#L247) |
| `_candidate_groups` | [289](../scripts/vntrcore.py#L289) |
| `split_stable_variants` | [310](../scripts/vntrcore.py#L310) |
| `generate_candidates` | [356](../scripts/vntrcore.py#L356) |
| `assign_reads` | [392](../scripts/vntrcore.py#L392) |
| `call_candidates` | [434](../scripts/vntrcore.py#L434) |
| `cluster_reads` | [514](../scripts/vntrcore.py#L514) |
| `canonical_motif` | [519](../scripts/vntrcore.py#L519) |
| `motif_orientation` | [523](../scripts/vntrcore.py#L523) |
| `discover_motifs` | [532](../scripts/vntrcore.py#L532) |
| `motif_dictionary` | [575](../scripts/vntrcore.py#L575) |
| `cigar_edit` | [619](../scripts/vntrcore.py#L619) |
| `decompose` | [631](../scripts/vntrcore.py#L631) |
| `decompose.convert` | [671](../scripts/vntrcore.py#L671) |
| `token_counts` | [697](../scripts/vntrcore.py#L697) |
| `build_graph` | [705](../scripts/vntrcore.py#L705) |
| `graph_evidence` | [731](../scripts/vntrcore.py#L731) |
| `dbg_evidence` | [755](../scripts/vntrcore.py#L755) |

### vntrdbg.py

| 定义 | 当前行号 |
|---|---:|
| `digest` | [20](../scripts/vntrdbg.py#L20) |
| `save_json` | [24](../scripts/vntrdbg.py#L24) |
| `tsv` | [28](../scripts/vntrdbg.py#L28) |
| `words` | [41](../scripts/vntrdbg.py#L41) |
| `raw_graph` | [46](../scripts/vntrdbg.py#L46) |
| `supported_records` | [77](../scripts/vntrdbg.py#L77) |
| `supported_records.enough` | [83](../scripts/vntrdbg.py#L83) |
| `compact_graph` | [122](../scripts/vntrdbg.py#L122) |
| `compact_graph.evidence` | [196](../scripts/vntrdbg.py#L196) |
| `write_gfa` | [212](../scripts/vntrdbg.py#L212) |
| `color` | [225](../scripts/vntrdbg.py#L225) |
| `topology` | [232](../scripts/vntrdbg.py#L232) |
| `paths_svg` | [267](../scripts/vntrdbg.py#L267) |
| `write_outputs` | [320](../scripts/vntrdbg.py#L320) |
| `metrics` | [341](../scripts/vntrdbg.py#L341) |
| `build_all` | [373](../scripts/vntrdbg.py#L373) |
| `main` | [401](../scripts/vntrdbg.py#L401) |

### vntrdbg_analysis.py

| 定义 | 当前行号 |
|---|---:|
| `read_tsv` | [29](../scripts/vntrdbg_analysis.py#L29) |
| `write_tsv` | [35](../scripts/vntrdbg_analysis.py#L35) |
| `occurrences` | [46](../scripts/vntrdbg_analysis.py#L46) |
| `branch_queries` | [59](../scripts/vntrdbg_analysis.py#L59) |
| `anchored_interval` | [137](../scripts/vntrdbg_analysis.py#L137) |
| `define_events` | [146](../scripts/vntrdbg_analysis.py#L146) |
| `assess_observation` | [189](../scripts/vntrdbg_analysis.py#L189) |
| `joint_support` | [222](../scripts/vntrdbg_analysis.py#L222) |
| `template_for` | [294](../scripts/vntrdbg_analysis.py#L294) |
| `primitive_templates` | [304](../scripts/vntrdbg_analysis.py#L304) |
| `primitive_units` | [313](../scripts/vntrdbg_analysis.py#L313) |
| `terminal_evidence` | [443](../scripts/vntrdbg_analysis.py#L443) |
| `main` | [496](../scripts/vntrdbg_analysis.py#L496) |

### vntrdbg_inputs.py

| 定义 | 当前行号 |
|---|---:|
| `rows` | [20](../scripts/vntrdbg_inputs.py#L20) |
| `intervals` | [26](../scripts/vntrdbg_inputs.py#L26) |
| `prepare` | [59](../scripts/vntrdbg_inputs.py#L59) |

### vntrdbg_report.py

| 定义 | 当前行号 |
|---|---:|
| `support_threshold` | [15](../scripts/vntrdbg_report.py#L15) |
| `preview` | [29](../scripts/vntrdbg_report.py#L29) |
| `main` | [40](../scripts/vntrdbg_report.py#L40) |

### vntrgraph.py

| 定义 | 当前行号 |
|---|---:|
| `Locus` | [33](../scripts/vntrgraph.py#L33) |
| `IndexedFastaReference` | [44](../scripts/vntrgraph.py#L44) |
| `IndexedFastaReference.__init__` | [47](../scripts/vntrgraph.py#L47) |
| `IndexedFastaReference.__enter__` | [58](../scripts/vntrgraph.py#L58) |
| `IndexedFastaReference.__exit__` | [61](../scripts/vntrgraph.py#L61) |
| `IndexedFastaReference.contig_length` | [64](../scripts/vntrgraph.py#L64) |
| `IndexedFastaReference.fetch` | [67](../scripts/vntrgraph.py#L67) |
| `validate_locus_id` | [86](../scripts/vntrgraph.py#L86) |
| `locus_output_dir` | [92](../scripts/vntrgraph.py#L92) |
| `load_source_namespaces` | [102](../scripts/vntrgraph.py#L102) |
| `read_catalog` | [118](../scripts/vntrgraph.py#L118) |
| `discover_bams` | [154](../scripts/vntrgraph.py#L154) |
| `reference_to_query` | [171](../scripts/vntrgraph.py#L171) |
| `flank_hit` | [190](../scripts/vntrgraph.py#L190) |
| `extract_sam_records` | [239](../scripts/vntrgraph.py#L239) |
| `write_tsv` | [250](../scripts/vntrgraph.py#L250) |
| `source_hashes` | [258](../scripts/vntrgraph.py#L258) |
| `write_gfa` | [263](../scripts/vntrgraph.py#L263) |
| `compress_tokens` | [281](../scripts/vntrgraph.py#L281) |
| `write_svg` | [294](../scripts/vntrgraph.py#L294) |
| `process_locus` | [348](../scripts/vntrgraph.py#L348) |
| `main` | [479](../scripts/vntrgraph.py#L479) |

### vntrmotifs.py

| 定义 | 当前行号 |
|---|---:|
| `primitive_unit` | [12](../scripts/vntrmotifs.py#L12) |
| `forms` | [21](../scripts/vntrmotifs.py#L21) |
| `canonical` | [26](../scripts/vntrmotifs.py#L26) |
| `definition` | [30](../scripts/vntrmotifs.py#L30) |
| `exact_runs` | [38](../scripts/vntrmotifs.py#L38) |
| `exact_occurrences` | [63](../scripts/vntrmotifs.py#L63) |
| `context_support` | [74](../scripts/vntrmotifs.py#L74) |
| `PathState` | [102](../scripts/vntrmotifs.py#L102) |
| `TwoPaths` | [109](../scripts/vntrmotifs.py#L109) |
| `TwoPaths.__init__` | [111](../scripts/vntrmotifs.py#L111) |
| `TwoPaths.offer` | [117](../scripts/vntrmotifs.py#L117) |
| `TwoPaths.trace` | [128](../scripts/vntrmotifs.py#L128) |
| `TwoPaths.audit` | [139](../scripts/vntrmotifs.py#L139) |
| `base_metrics` | [150](../scripts/vntrmotifs.py#L150) |
| `compound_annotations` | [170](../scripts/vntrmotifs.py#L170) |
| `observed_compounds` | [198](../scripts/vntrmotifs.py#L198) |

### vntrreads.py

| 定义 | 当前行号 |
|---|---:|
| `cigar_interval` | [37](../scripts/vntrreads.py#L37) |
| `near_anchor` | [80](../scripts/vntrreads.py#L80) |
| `extract_record` | [95](../scripts/vntrreads.py#L95) |
| `row_molecule_key` | [145](../scripts/vntrreads.py#L145) |
| `read_bam` | [150](../scripts/vntrreads.py#L150) |
| `collect_reads` | [218](../scripts/vntrreads.py#L218) |
| `annotation_results` | [231](../scripts/vntrreads.py#L231) |
| `motif_templates` | [241](../scripts/vntrreads.py#L241) |
| `annotate` | [266](../scripts/vntrreads.py#L266) |
| `annotate.update` | [294](../scripts/vntrreads.py#L294) |
| `write_table` | [343](../scripts/vntrreads.py#L343) |
| `svg_read` | [348](../scripts/vntrreads.py#L348) |
| `svg_read.text` | [352](../scripts/vntrreads.py#L352) |
| `write_outputs` | [427](../scripts/vntrreads.py#L427) |
| `run_locus` | [490](../scripts/vntrreads.py#L490) |
| `main` | [530](../scripts/vntrreads.py#L530) |

### vntrstr.py

| 定义 | 当前行号 |
|---|---:|
| `parse_structure` | [21](../scripts/vntrstr.py#L21) |
| `validate_config` | [43](../scripts/vntrstr.py#L43) |
| `reconstruct` | [53](../scripts/vntrstr.py#L53) |
| `format_blocks` | [59](../scripts/vntrstr.py#L59) |
| `reverse_blocks` | [65](../scripts/vntrstr.py#L65) |
| `reverse_blocks.rc` | [67](../scripts/vntrstr.py#L67) |
| `_blocks` | [80](../scripts/vntrstr.py#L80) |
| `annotate` | [94](../scripts/vntrstr.py#L94) |
| `path_record` | [148](../scripts/vntrstr.py#L148) |
| `annotate_rows` | [153](../scripts/vntrstr.py#L153) |
| `svg_blocks` | [168](../scripts/vntrstr.py#L168) |
| `write_paths_svg` | [194](../scripts/vntrstr.py#L194) |
| `write_outputs` | [226](../scripts/vntrstr.py#L226) |

### vntrstrgraph.py

| 定义 | 当前行号 |
|---|---:|
| `compact_walk` | [21](../scripts/vntrstrgraph.py#L21) |
| `_token` | [28](../scripts/vntrstrgraph.py#L28) |
| `build_graph` | [34](../scripts/vntrstrgraph.py#L34) |
| `write_gfa` | [140](../scripts/vntrstrgraph.py#L140) |
| `_lines` | [156](../scripts/vntrstrgraph.py#L156) |
| `write_topology` | [161](../scripts/vntrstrgraph.py#L161) |
| `write_paths_svg` | [179](../scripts/vntrstrgraph.py#L179) |
| `write_outputs` | [236](../scripts/vntrstrgraph.py#L236) |

### vntrstructure.py

| 定义 | 当前行号 |
|---|---:|
| `support_key` | [30](../scripts/vntrstructure.py#L30) |
| `seed_dictionary` | [36](../scripts/vntrstructure.py#L36) |
| `normalize_dictionary` | [49](../scripts/vntrstructure.py#L49) |
| `add_reference_units` | [70](../scripts/vntrstructure.py#L70) |
| `exact_tandem_units` | [106](../scripts/vntrstructure.py#L106) |
| `segment` | [112](../scripts/vntrstructure.py#L112) |
| `segment.update` | [139](../scripts/vntrstructure.py#L139) |
| `segment.convert` | [161](../scripts/vntrstructure.py#L161) |
| `segment_details` | [188](../scripts/vntrstructure.py#L188) |
| `annotate_rows` | [193](../scripts/vntrstructure.py#L193) |
| `learn_variants` | [202](../scripts/vntrstructure.py#L202) |
| `run_blocks` | [262](../scripts/vntrstructure.py#L262) |
| `compound_blocks` | [275](../scripts/vntrstructure.py#L275) |
| `summarize` | [281](../scripts/vntrstructure.py#L281) |
| `summarize.node` | [284](../scripts/vntrstructure.py#L284) |
| `summarize.add` | [287](../scripts/vntrstructure.py#L287) |
| `summarize.edge` | [290](../scripts/vntrstructure.py#L290) |
| `summarize.finish` | [327](../scripts/vntrstructure.py#L327) |
| `color` | [338](../scripts/vntrstructure.py#L338) |
| `graph_svg` | [349](../scripts/vntrstructure.py#L349) |
| `rows_svg` | [373](../scripts/vntrstructure.py#L373) |
| `write_gfa` | [399](../scripts/vntrstructure.py#L399) |
| `emit` | [414](../scripts/vntrstructure.py#L414) |
| `main` | [506](../scripts/vntrstructure.py#L506) |

### vntrsummary.py

| 定义 | 当前行号 |
|---|---:|
| `identity` | [18](../scripts/vntrsummary.py#L18) |
| `blocks` | [23](../scripts/vntrsummary.py#L23) |
| `structure_text` | [36](../scripts/vntrsummary.py#L36) |
| `reference_report` | [42](../scripts/vntrsummary.py#L42) |
| `is_primary` | [66](../scripts/vntrsummary.py#L66) |
| `strong_boundary` | [71](../scripts/vntrsummary.py#L71) |
| `local_candidates` | [76](../scripts/vntrsummary.py#L76) |
| `structural_summary` | [111](../scripts/vntrsummary.py#L111) |
| `structural_summary.counts` | [167](../scripts/vntrsummary.py#L167) |
| `summary_svg` | [214](../scripts/vntrsummary.py#L214) |

### vntrterminal.py

| 定义 | 当前行号 |
|---|---:|
| `validate_terminals` | [18](../scripts/vntrterminal.py#L18) |
