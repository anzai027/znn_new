# 预测型 GNN 实验使用指南

实验结果和原因分析见 [实验报告](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/报告.md)。

与旧 VFCR 的统一口径补充比较见 [VFCR对照](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/VFCR对照.md)。新增 pp_vfcr 补跑保存于 vfcr；pp_vcheck 只核对已有数据，pp_vplots 绘制比较图。VFCR-S 两组完成，原式一组达到调用预算；原式没有完整区间指标。pp_vfcr 会拒绝覆盖已有 vfcr/all.mat。

这组实验按 [聊天建议](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/sources/聊天建议.md) 比较 M3、M5、M6、M7，先使用 Example 2 和零噪声。旧模型和旧实验结果不修改。模型公式与诊断字段见 [系统说明](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/系统说明.md)。

## 运行实验

在 MATLAB 中进入本文件夹，执行：

```matlab
cd('C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006')
run_predictor
```

默认保存到 `data`。重新实验时，传入一个尚未使用的目录名，避免覆盖已有结果：

```matlab
run_predictor('data_run2')
```

程序先运行四个模型的主实验，再运行四个模型的严格精度检查，共八组。若目录已有 `all.mat`，程序会拒绝运行；尚未完成的目录仍可能含逐组数据，所以重跑也应使用新名字。

参数集中在 [pp_config.m](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/pp_config.m)。主实验容差为 `RelTol=1e-7、AbsTol=1e-9`，严格检查为 `1e-9、1e-11`。默认时间为 0～10 秒；一般输出间隔 0.0005 秒，两处参考速度峰附近加密到 1e-5 和 1e-6 秒。每组最多 300 秒或 500000 次右端函数调用。

数据目录保存 `reference.mat`、`run_01.mat`～`run_08.mat`、`all.mat`、`metrics.csv`、`precision.csv` 和 `run.log`。MAT 文件包含原始轨迹、运行设置和源码校验值。完整积分只表示程序跑到终点，求解是否准确仍要看误差、残差和约束。

## 补充对照、绘图与验收

主实验保存后，可执行：

```matlab
pp_extra
pp_figures
verify_results('data')
```

[pp_extra.m](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/pp_extra.m) 运行两组严格容差对照：M6 的 kappa 改为 0.01；M7 的 scale 改为 [1,1,1]，即 D=I。其他设置沿用已保存的主实验配置，结果写到 `extra`。这个函数固定读取 `data`；再次执行会重跑并覆盖 `extra`。

[pp_figures.m](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/pp_figures.m) 默认读取 `data/all.mat` 中 main 的四组轨迹，保存 `main、switch、prediction` 三张 PNG 和 FIG 到 `figures`。它不会自动绘制严格容差或补充对照。若图目录已有同名图，程序拒绝覆盖，可显式指定新目录：`pp_figures('data','figures_run2')`。

[verify_results.m](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/verify_results.m) 独立重算参考解和指标，检查源码校验值、轨迹、M7 共同限速及精度一致性；它不重新运行 ODE，验收输出写到 `verification`。

若主实验另存到 `data_run2`，这些工具不会自动切换数据目录。`pp_extra` 仍固定使用 `data`；绘图和验收必须显式传入名称，例如 `pp_figures('data_run2','figures_run2')`、`verify_results('data_run2')`。不要把默认 `pp_extra` 的结果当成 `data_run2` 的补充对照。

## 单独调用系统

```matlab
addpath('C:/Users/anzai/Documents/znn_new')
c = pp_config();
p = example2_problem();
t = c.start;
g = c.g0;
[dg,a] = pp_system(t,g,p,c,"M3");
[dg,a] = m5_system(t,g,p,c);
[dg,a] = m6_system(t,g,p,c);
[dg,a] = m7_system(t,g,p,c);
```

`dg` 是当前状态速度，`a` 是残差、雅可比、预测速度、阻尼和限速等诊断。状态顺序为 `g=[x;mu1;mu2]`。需要只积分一个模型时，可以调用 `r=pp_run("M5",p,c)`，其结果不会自动写入文件。

- M3：预条件平滑反馈，固定阻尼，没有预测项。
- M5：在 M3 上加入题目时间变化的差分预测。
- M6：在 M5 上按照原 J 的最小奇异值调整阻尼。
- M7：在 M6 上加入变量尺度和共同限速，阻尼依据 `B=JD`。

[pp_run.m](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/pp_run.m) 只把 `G、h、P、u、Q、v` 六个题目系数字段传给系统。解析导数 `dG、dh、dP、du、dQ、dv` 只用于参考速度和独立预测诊断，不进入模型。时间差分本身仍然是一种导数信息估计，不能称为完全不使用导数信息。

## 初始化与速度预算

预测差分固定同一个 g。`t>=h` 时读取 t 和 t-h；初始化 `0<=t<h` 时读取 t 和 t+h。这一小段前向差分使用未来题目系数，是非因果处理，系统用 `a.boundary` 标记。默认 h=1e-5 秒。系统不保存历史，不依赖求解器调用时间的先后顺序。

M7 在 z 空间使用 gamma=200，默认 `D=diag(1,1,10,10,10,10,10)`。转换回物理 g 空间后，反馈速度范数的上界可达到 2000；M3、M5、M6 的物理反馈上界为 200。因此四个模型没有相同的物理反馈速度预算，比较结果时必须注明这个差别。预测项还会增加总速度。

M7 先计算包含预测项的整个总速度，再乘同一个限速系数 alpha；三个块不会分别归一化。默认 x、mu1、mu2 的块速度上限为 100、5000、5000。共同缩放在静态题目中保留反馈下降方向，不能据此声称时变总能量始终下降。

## 公式检查与失败记录

执行 `audit` 可重做短公式检查。2026-10-06 的 MATLAB R2025b 检查已通过 320 个测点：稳定 PFB 雅可比最大相对差分误差 2.12672e-10，独立矩阵公式最大相对误差 3.4702e-13，固定状态时间差分与直接重算差异为 0。检查还覆盖移动零点、无等式或无不等式、M7 共同限速，并禁止读取六个解析导数字段。证据在 `audit.csv`、`audit.mat` 和 [audit.log](C:/Users/anzai/Documents/znn_new/experiments/gnn_predictor_20261006/audit.log)。

`launch_failed` 保留初次启动时的八个调用数为 0 的记录。它们来自实验脚本未预先声明 `sol` 的程序错误，积分没有开始，因此排除出模型性能比较。该脚本问题已修复，正式实验已在新的 `data` 中重跑。这些启动错误不代表模型无法收敛。

正式实验中的 `solver_failure`、`budget` 或 `exception` 仍需保存并单独说明；未完成轨迹不计算整个观察区间的性能指标。
