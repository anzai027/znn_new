# VFCR-ZNN Example 2：Fig. 3复现说明

## 运行方法

将下列文件放在同一MATLAB目录：

- `run_example2_fig3.m`
- `example2_problem.m`
- `my_system.m`
- `vfcr_phi.m`
- `vfcr_residual.m`

运行：

```matlab
run_example2_fig3
```

程序将比较BSAF、SBPAF、NSBPAF和NFTAF在以下两种情形的残差：

1. 无噪声，对应论文Fig. 3(a)；
2. 余弦噪声 `zeta_i(t)=cos(t)`，对应论文Fig. 3(b)。

## Example 1与Example 2的关系

### 通用代码

- `my_system.m`：VFCR-ZNN式(20)的ODE右端；
- `vfcr_phi.m`：论文列出的四种激活函数；
- `vfcr_residual.m`：计算KKT/PFB残差；
- 参数 `r1=r2=a=lambda1=lambda2=1`、`p=q=0.5`、`delta=0.0001`；
- `G,h,P,u,Q`及其时间导数。

### Example 2专用代码

- `example2_problem.m`：把 `v`设置为 `1.2*ones(4,1)`；
- `run_example2_fig3.m`：四种激活函数、无噪声/余弦噪声的循环与Fig. 3绘图。

### 核心差别

Example 1使用 `v=1e8*ones(4,1)`，边界几乎不起限制作用。

Example 2使用 `v=1.2*ones(4,1)`，对应：

```text
-1.2 <= x1(t), x2(t) <= 1.2
```

## 论文参数与本组数值设置的区别

论文明确给出：

- 四种激活函数及其公式；
- NFTAF参数 `a=1,p=0.5,q=0.5`；
- 其他激活函数参数 `alpha=4,beta=0.4`；
- 初始 `g(0)`从 `[0,1]`随机选取；
- 无噪声与 `cos(t)`噪声两种情形。

论文未公开随机种子、ODE求解器、采样点数和误差容限。本代码固定随机种子1，并使用 `ode15s`及显式容差，以保证本次复现实验可重复。因此初始瞬态不保证与论文逐点完全相同。

## 关于Fig. 4

论文Fig. 4比较OZNN、PTCR-ZNN和VFCR-ZNN在六种噪声下的结果。当前项目的 `my_system.m`只实现VFCR-ZNN，因此本代码不冒充完整Fig. 4复现。若要继续复现Fig. 4，还需分别实现OZNN式(15)和PTCR-ZNN，并核对PTCR-ZNN参考文献中的完整参数。
