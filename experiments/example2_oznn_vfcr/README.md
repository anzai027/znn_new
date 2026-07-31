# Example 2：OZNN 与 VFCR-ZNN 对比

这个文件夹是一组独立实验，根目录代码不会被覆盖。

## OZNN 推导

论文公式（14）为

\[
\xi(t)=H(t)g(t)+\vartheta(t).
\]

PFB 项中的 \(\sigma\) 也含有 \(g\)，所以完整求导后是

\[
\dot\xi(t)=J(t)\dot g(t)+M(t)g(t)+\varsigma(t).
\]

将论文公式（15）

\[
\dot\xi(t)=-\lambda\Psi(\xi(t))
\]

代入，得到

\[
\boxed{
J(t)\dot g(t)
=-M(t)g(t)-\varsigma(t)
-\lambda\Psi\!\left(H(t)g(t)+\vartheta(t)\right)
}.
\]

论文 Table I 规定 OZNN 使用线性激活函数 \(\Psi(\xi)=\xi\)，因此本实验实现为

\[
\dot g(t)=J(t)^{-1}
\left[-M(t)g(t)-\varsigma(t)-\lambda\xi(t)\right].
\]

代码使用 `J\rhs`，没有直接计算矩阵逆。

## 公平设置

- 两个模型共用同一个 `g0`。
- VFCR-ZNN 的额外积分初值为 `z0=zeros(7,1)`。
- 两个模型共用同一个 Example 2、`delta`、时间点、求解器容差和无噪声条件。
- VFCR-ZNN 使用论文的 NFTAF，OZNN 使用 Table I 的线性激活函数。
- 论文没有给出 Example 2 的 OZNN 参数 `lambda`，本实验固定为 `10`，这是可复现设置，不冒充论文明示参数。

## 运行与输出

在 MATLAB 中进入本文件夹并运行：

```matlab
run_compare
```

运行结束后会生成：

- `example2_oznn_vfcr.fig`：MATLAB 原始图。
- `example2_oznn_vfcr.png`：图像预览。
- `all.mat`：由 `save all` 保存的全部实验数据。

图中纵轴是完整残差 \(\|\xi(t)\|_2\)。脚本还会检查 \(t\ge5\) 时两种模型的残差和状态差。
