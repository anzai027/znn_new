# VFCR-ZNN 示例1

请将以下文件放在同一个MATLAB目录中：

- `run_example1.m`
- `example1_problem.m`
- `my_system.m`
- `znn_ode.m`
- `vfcr_residual.m`

运行`run_example1.m`。

ODE状态为：

```matlab
Y = [g; z]
```

其中，`g=[x;mu1;mu2]`包含7个元素，历史积分状态`z`也包含7个元素。
因此，公式（20）作为一个14维ODE进行积分。

默认设置`phi_eps=1e-4`，仅在零点附近对NFTAF进行平滑近似。
这样可以避免`ode15s`在不满足利普希茨条件的零点附近采用极小步长。
如果要使用严格的公式（18），请在`run_example1.m`中设置`phi_eps=0`，但运行时间可能明显增加。

论文没有报告随机种子、ODE求解器、误差容限和零点处的数值处理方法。
因此，复现目标是得到相同的定性固定时间收敛行为，并与时变KKT解保持一致，而不是让初始瞬态与论文逐像素相同。
