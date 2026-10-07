# M8a-BDF2、delta 与噪声实验

本轮继续 `gnn_m8_20261007`，研究 Example 2。旧目录和原始系统代码没有修改。

- [实验报告](报告.md)：结果、原因和下一步建议。
- [系统说明](系统说明.md)：公式、启动规则、噪声位置和比较限制。
- `results/metrics.csv`：37 组主实验的时间加权统计。
- `results/run_XX.mat`：每组实际分段 ODE 轨迹、配置和误差数组。
- `extra/metrics.csv`：8 组容差、积分步长和正则化基线复核。
- `verification.csv`：重新读取实际轨迹、更密网格和参考根修正后的比较。
- `reference_check.csv`：三个 delta 的参考残差、病态程度和根速度。
- `bias_tight.csv`：用自适应积分收敛复核后的原题 RMSE。
- `other_check.csv`：其余 10 组轨迹和指标的独立复算。
- `joint/metrics.csv`：小 delta 与相同测量噪声的联合检查。
- `prediction.csv`：在相同参考状态上比较差分预测误差。
- `iid.csv`：独立样本噪声的差分放大测试，**不是白噪声 ODE 实验**。
- `audit.csv`：公式检查及原 VFCR 代码等价性检查。
- `figures/`：真实实验数据的 PNG 与可编辑 FIG。

在仓库根目录启动 MATLAB：

```matlab
addpath('experiments/gnn_bdf2_20261007');
bd_setup;
bd_audit;
bd_study('new_results');
```

现有数据不得覆盖。`bd_extra`、`bd_verify`、`bd_figures` 使用本轮 `results` 的固定路径；要对新一轮数据运行这些函数，应先调整目录。`bd_probe` 可以单独重新计算预测诊断。

单独调用模型：

```matlab
bd_setup;
c = bd_config();
p = example2_problem();
dg = m8a_bdf2_system(1,c.g0,p,c);
```

`data/` 是最初调试时的记录：积分到达终点，但读取空采样区间时报错；修复后完整实验在 `results/`。不要将 `data/` 中的保存失败当成模型失败。

代码依赖旧目录的 `pp_parts`、`pp_reference`、`pn_config` 和 `fsmooth`。这些依赖是公式与基线的共同来源。
