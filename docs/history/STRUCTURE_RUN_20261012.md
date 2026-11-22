# r4 结构图重跑记录

2026-10-12 完成主程序部署和 Slurm 提交。最后检查时三个作业均为 `PENDING (Priority)`，尚未生成本轮全量结果。

| 点位 | 作业 ID | BED 候选单元长度 |
| --- | --- | --- |
| IL1RN | 159116 | 86 bp |
| PER3 | 159117 | 54 bp、54 bp（不同序列） |
| ABCA7 | 159118 | 25 bp、26 bp |

以上候选种子在 r3 提取的 HiFi 序列中都有直接观测依据：按 BAM/read 行计，包含精确种子的行数分别为 IL1RN 251；PER3 188、183；ABCA7 66、85。这个种子来源统计不是最终图的去重支持数。

本地 61 项单元测试通过；三个点位各取前 12 条 HiFi 行的有限页面试算通过。PER3 该试算未解释比例为 0.913%，不代表本轮全量结果或确定的生物学正确率。浏览器分页/读段清单、搜索、变体切换、放大和节点点击已检查。

服务器主项目：`/home/data/user/stzhu/重复区域分析/vntr-point-graph`

- 本轮不可变执行快照：`runs/r4-structure-20261012/`
- 结果：`runs/r4-structure-20261012/results/`
- 日志：`runs/r4-structure-20261012/logs/`
- 作业表：`runs/r4-structure-20261012/structure_jobs.tsv`
- 更新前备份：`runs/r4-before-update.tgz`
- r3 结果未覆盖。

完成后优先打开每个位点的 `位点.structure.html`，不要把 r3 的 `位点.reads.html` 当作新版。
下载时保留整个位点目录，可同时使用全部 SVG 分页、GFA 和原始记录链接。

每个作业会先执行 61 项回归测试及合成 BAM 集成测试，完成真实数据后再验证 `verification.json`。
验收要求是 read 清单和提取序列相对 r3 完全一致、精确序列 GFA 重构全部通过、共享节点支持数按不同 sample/QNAME 重新核算。

```bash
squeue -j 159116,159117,159118
sacct -j 159116,159117,159118 --format=JobID,State,ExitCode,Elapsed
```

详细原因、算法、参数和限制见 `STRUCTURE_VIEW.md`；种子输入见 `catalogs/VNTR_structure_seeds.bed`。
