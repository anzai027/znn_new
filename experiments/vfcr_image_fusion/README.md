# VFCR-ZNN 图像融合消噪复现

这个文件夹只复现论文第 VI-A 节中的 VFCR-ZNN，不包含 FFT-ZNN 和 PTCR-ZNN。

## 实验内容

1. BSDS300：20 张带噪图像融合。
   - 高斯噪声：均值 0，方差 0.05、0.08、0.10。
   - 椒盐噪声：密度 0.05、0.10、0.20。
2. 512×512 船图：方差 0.05 的高斯噪声，分别融合 20 张和 50 张。
3. BrainWeb：第 40 层使用高斯噪声，第 100 层使用 Rician 噪声，各融合 50 张。

## 运行方法

先在 PowerShell 中下载数据：

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\download_data.ps1
```

然后在 MATLAB 中运行：

```matlab
run_all
```

也可以在 PowerShell 中直接运行：

```powershell
matlab -batch "cd('此文件夹的绝对路径'); run_all"
```

## 输出文件

- `output/fig`：MATLAB 可编辑的 `.fig` 图。
- `output/image`：便于直接查看的 `.png` 图。
- `output/table`：PSNR、SSIM、残差、权重约束等表格。
- `output/all.mat`：主脚本结束时用 `save` 保存的全部工作区数据。
- `output/run.log`：完整运行日志。

## 与论文一致的部分

- 图像矩阵按论文式 (1)-(5) 构造。
- 二次规划矩阵使用 `G=W'*W/S`。
- 权重满足和为 1、非负的约束。
- VFCR-ZNN 使用论文式 (18) 的 NFTAF 和式 (20)。
- 参数为 `r1=r2=a=lambda1=lambda2=1`、`p=q=0.5`、`delta=1e-4`。
- 代码直接使用 NFTAF，不使用平滑替代函数。

## 论文没有公开的细节

- 论文没有公开随机种子、ODE 求解器和容差，本复现固定随机种子并使用 `ode15s`。
- 论文没有说明 512×512 原图的文件名；根据图 8 外观采用 USC-SIPI `boat.512`。
- 论文没有给出 MRI 外加噪声强度；本复现用论文单张噪声图的 PSNR 反推方差 0.06 和 0.0385。
- 论文 Table III 的 20 张结果与方差 0.05 的常规计算不一致，因此不把逐像素一致当作成功标准。

成功标准是：KKT 残差下降、权重和接近 1、最小权重非负，并且融合图的 PSNR 和 SSIM 明显高于单张带噪图。
