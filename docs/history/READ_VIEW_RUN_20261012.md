# 全部 read 版：部署与运行记录

## 当前状态

2026-10-12 已部署到服务器主项目：
`/home/data/user/stzhu/重复区域分析/vntr-point-graph`。

本地 50 项单元测试通过；浏览器合成页面检查通过：23 条 read 的 20+3 分页、搜索、完整节点表及脚本执行均正常。
这些是代码和合成数据验证，不代表真实 BAM 已经跑完。

Slurm 已接受且仅提交以下三个点位；最后检查均为 `PENDING (Priority)`，CPU 队列节点全部占用。

| 点位 | 作业 ID | 状态 |
| --- | --- | --- |
| IL1RN | 159082 | 排队 |
| PER3 | 159083 | 排队 |
| ABCA7 | 159084 | 排队 |

真实 BAM 全量结果和真实保留率尚未生成，不能报告为验收完成。
IL1RN 作业在实际数据分析前会先执行完整测试和真实格式的合成 BAM 集成测试。

## 快照与保护

- 独立执行快照：`runs/r3-all-reads-20261012/`
- 作业清单：`runs/r3-all-reads-20261012/reads_jobs.tsv`
- 输出：`runs/r3-all-reads-20261012/results/`
- 日志：`runs/r3-all-reads-20261012/logs/`
- 更新前源码备份：`runs/r3-before-update.tgz`
- 旧的 `runs/r2/results/` 未覆盖。

## 图生成后如何看

将一个点位的完整结果目录下载到本机，双击对应的 `IL1RN.reads.html`、`PER3.reads.html` 或 `ABCA7.reads.html`。
每页 20 条只是分页显示，不是只选 20 条；可搜索 read ID、按样本和平台筛选。
展开某条 read 的“全部节点与完整碱基”，即可查看有序节点、坐标和真实碱基。
静态图见 `svg/index.html`，所有分页都保留；GFA 与 HTML/SVG 使用相同节点 ID。

## 服务器查看状态

```bash
squeue -j 159082,159083,159084
sacct -j 159082,159083,159084 --format=JobID,State,ExitCode,Elapsed
cd '/home/data/user/stzhu/重复区域分析/vntr-point-graph/runs/r3-all-reads-20261012'
tail -n 30 logs/*.out logs/*.err
```

完成时应核查每个位点的 `status.json`、Slurm 的 `COMPLETED/0:0`，再独立比对每个 BAM 的 read 数与输出清单。
不要重复执行提交脚本；现有三个作业已排队，重复提交会重复计算。

具体提取和节点分解逻辑见同目录 `READ_VIEW.md`。
