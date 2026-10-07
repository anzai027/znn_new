# M8 / M8a 预测项比较

按 2026-10-07 读取的「融合GNN与VFCR-ZNN」对话实现两个新模型，并与已保存的 M3、M5、M6、M7 和 VFCR-S 轨迹比较。

13 次新增积分已完成，独立验收通过。M8a 的尖峰明显改善；更高精度下 VFCR-S 仍更准。

- [报告](报告.md)：主比较、小对照、积分精度复核和改进建议。
- [系统说明](系统说明.md)：公式、参数、命名和实验边界。
- [聊天来源](sources/聊天建议.md)：本轮依据的完整回答快照。
- `m8_system.m`：分离预测和纠错阻尼。
- `m8a_system.m`：直接解预测方程，保留 M6 纠错。
- `data/metrics.csv`：全部 9 次新增积分的指标。
- `data/comparison.csv`：新旧模型严格容差对照及来源路径。
- `precision/metrics.csv`：更严格容差的四组配对复核。
- `verification/`：独立指标、来源与加密采样验收。
- `figures/`：图表及可编辑 FIG。

在仓库根目录运行 MATLAB：

```matlab
addpath('experiments/gnn_predictor_20261006','-begin');
addpath('experiments/gnn_m8_20261007','-begin');
run_m8('data');
precision('precision');
m8_verify('data');
m8_pverify('precision');
probe;
audit;
figures;
figures('precision');
```

已有结果时脚本拒绝覆盖；复跑请改成新的文件夹名称。依赖仓库根目录的 `example2_problem.m` 和旧实验目录的 `pp_parts.m`、`pp_config.m`、`pp_reference.m`，实际使用文件的 SHA256 随数据保存。

主实验只有 Example 2 + zero noise。两模型各做普通和严格容差；额外五组分别检查相同纠错、两个预测阻尼和缩小差分步长。旧 GNN 模型通过来源核验后复用；VFCR-S 另做了一次更高精度积分，避免把积分器的误差误认为模型性能。

`figures` 已保存的图拒绝覆盖。原不平滑 VFCR 的已保存结果只运行到 1.56755 s，完整误差对比采用 VFCR-S；它使用解析系数导数，新 GNN 仅用系数值差分。
