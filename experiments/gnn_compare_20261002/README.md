# FTCGNN system

本阶段只编写 system，没有运行 ODE 实验，没有生成实验结果。

已使用 MATLAB R2025b 的 `checkcode` 完成静态语法检查，没有诊断信息；记录在 `logs/syntax.log`。

文件 `ftcgnn_system.m` 实现附件《GNN Model With Robust Finite-Time Convergence for Time-Varying Systems of Linear Equations》中的原版 FTCGNN。

## 公式与代码

论文印刷第 4788 页（PDF 第 3 页）式（9）为：

\[
\|A^T(Ax-b)\|_2^p\dot x=-\gamma A^T(Ax-b),\qquad \gamma>0,\quad 0<p<2.
\]

令 `xi=A*x-b`、`s=A.'*xi`、`r=norm(s,2)`，当 `r>0` 时：

```matlab
dx = -gamma*s/r^p;
```

其中 `dx` 是状态此刻的变化速度；它不是整个时间段的求解结果。

论文印刷第 4790 页（PDF 第 5 页）式（18）加入噪声后为：

\[
\dot x=-\gamma\frac{A^T(Ax-b)}{\|A^T(Ax-b)\|_2^p}+\sigma(t).
\]

代码中的 `problem.noise(t)` 对应这里的 \(\sigma(t)\)，它加在整个归一化梯度项之后；不传 `noise` 字段时表示无噪声。

## 输入和调用

函数形式与现有 system 一致：时间、状态、状态长度、题目系数、模型参数依次传入。

```matlab
dx = ftcgnn_system(t,x,l,problem,gamma,p);
```

| 输入 | 含义 |
|---|---|
| `t` | 当前时间 |
| `x` | 长度为 `l` 的状态向量 |
| `l` | 状态变量个数，也就是 `A` 的列数 |
| `problem.A(t)` | 当前时刻的系数矩阵 |
| `problem.b(t)` | 当前时刻的右端向量 |
| `problem.noise(t)` | 可选的加性噪声，标量会扩展为 `l` 个相同分量 |
| `gamma` | 正的增益参数 |
| `p` | 分母中梯度范数的指数 |

`A`、`b` 和 `x` 应为有限的实数数据；`A` 可以是矩形矩阵，但维数必须相容。

论文第 4791 页（PDF 第 6 页）的 Remark 3 指出：静态问题取 `0<p<2`，时变问题的收敛定理取 `p=1`；保证收敛还需要相应的矩阵秩、解的存在性、增益和变化速度条件。

`gamma` 和 `p` 没有写死，因为论文的不同实验使用不同参数；例如时变无噪声对比取 `gamma=100,p=1`，时变有噪声对比取 `gamma=200,p=1`（第 4793 页，PDF 第 8 页）。

## 零点处理

论文式（9）在 `s=0` 时变成 `0=0`，不会唯一指定状态速度；程序在这个点把无噪声的梯度项取为零，避免 `0/0`。

这只是零点处的一种数值取值；对于时变问题，它不能代替对不连续动力学和到达后持续跟踪的严格论证。

程序没有给分母增加小量，没有平滑梯度，也没有加入积分、系数导数或额外激活函数。

## 与之后的 TVQP 对比的关系

这个原版 system 求解的是 `A(t)*x(t)=b(t)`；现有新模型 `gnn_system.m` 求解的是包含不等式互补条件的 KKT/PFB 方程，输入题目类型不同。

VFCR 中的 `theta` 含有依赖状态的平方根项，所以不能直接把 `H` 当成原版 FTCGNN 的 `A`、把 `-theta` 当成外部给定的 `b`。

如果之后把原版中的 `A.'*(A*x-b)` 换成 `J.'*xi`，那就是向 TVQP 的扩展；相同增益、相同 `p` 下，它也正是你现在提出的新模型，因此不能把两份相同的动力学当成两个独立方法来比较。

原 VFCR-ZNN 和新 GNN 文件均未修改；后续实验会直接调用现有代码。

## 当前文件夹

- `ftcgnn_system.m`：本阶段写好的模型。
- `sources/ftcgnn.pdf`：用户提供的 FTCGNN 论文副本。
- `sources/vfcr.pdf`：找到的 VFCR-ZNN 原论文副本。
- `sources/*.txt`：论文文字提取，用于定位公式。
- `sources/ftcgnn_p*.png`：核对公式时使用的论文页面图。
- `data`、`figures`、`logs`：留给后续实验，当前没有实验结果。

后续运行实验时，才会使用 MATLAB `save` 保存数据，并保存可打开的 `.fig` 和 `.png` 图像。
