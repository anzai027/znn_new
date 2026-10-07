function bd_report()
% 根据实际结果生成可追溯的中文实验报告。
folder = bd_setup(); s = load(fullfile(folder,'results','all.mat'),'stats'); z = s.stats;
s = load(fullfile(folder,'extra','all.mat'),'stats'); e = s.stats;
s = load(fullfile(folder,'joint','all.mat'),'stats'); j = s.stats;
v = readtable(fullfile(folder,'verification.csv'));
check = readtable(fullfile(folder,'reference_check.csv'));
peak = readtable(fullfile(folder,'peak_check.csv'));
adaptive = readtable(fullfile(folder,'bias_tight.csv'));
prior = readtable(fullfile(folder,'bias_adaptive.csv'));
assert(all(z.complete) && all(e.complete) && all(j.complete) && height(v)==37);
f = fopen(fullfile(folder,'报告.md'),'w','n','UTF-8'); assert(f>=0);
guard = onCleanup(@() fclose(f));
fprintf(f,'# M8a-BDF2、delta study 与 noise study\n\n');
fprintf(f,'本轮完成 **37 组主实验、8 组复核、2 组 delta 与噪声联合检查，共 47 组完整 ODE 实验**。另有公式审计、参考解检查、差分预测诊断和独立样本噪声放大测试。日期 2026-10-07。\n\n');
fprintf(f,'**核心结论：二阶后向差分成功降低无噪声跟踪误差；减小 delta 改善原题精度，却使乘子尖峰更敏感；测量噪声出现后，高精度优势基本被噪声淹没，固定预测阻尼并没有成为可靠的抗噪改进。**\n\n');
fprintf(f,'依据[最新聊天建议](https://chatgpt.com/g/g-p-6ab9211b1de081919276a3f4e721199e-gnngai-jin/shared/c/6a9e7bd7-8730-83ec-9027-9b4f35cd5f38)，沿用 Example 2。公式和噪声位置见[系统说明](系统说明.md)，复现方法见[README](README.md)。\n\n');
fprintf(f,'## 1. 无噪声：BDF2 有效，但不能只比较最后几位数字\n\n');
fprintf(f,'下表统计 5–10 s，容差均为 1e-11/1e-13、delta=1e-4。xd 是相对扰动根的 x RMSE，gpeak 是完整状态最大误差；VFCR-S 使用解析导数，其 h 不参与计算。\n\n');
print(f,z(1:7,{'model','order','h','xd','gpeak','respeak','calls'}));
fprintf(f,'同样 h=1e-5，M8a-BDF2 的 x 跟踪 RMSE 比一阶 M8a 改善约 **%.0f 倍**。BDF2 的 h 从 1e-4 改成 1e-5 时改善约 %.1f 倍；再改成 1e-6 则没有改善，已受到浮点相消、积分误差和参考误差影响。不能继续套用“h 越小越好”。\n\n',z.xd(1)/z.xd(3),z.xd(2)/z.xd(3));
fprintf(f,'预测项的独立审计中，步长减半时误差比的中位数约 4，符合二阶截断误差。第二个切换时刻，BDF2 的预测速度误差从 h=1e-4 的约 2.21e-4 降至 h=1e-5 的约 2.01e-6；h=1e-6 反而约 7.66e-6。完整记录见[prediction.csv](prediction.csv)。\n\n');
fprintf(f,'![无噪声比较](figures/zero.png)\n\n');
fprintf(f,'### 更严格容差复核\n\n');
fprintf(f,'将两种 FD2 模型和解析 VFCR-S 的容差再收紧到 1e-12/1e-14：\n\n');
print(f,e([1,2,8],{'model','rtol','atol','xd','gpeak','respeak','calls'}));
fprintf(f,'本配置下 M8a-BDF2 的计算结果仍小于两个 VFCR 对照，但它的误差已经接近数值底部，尖峰还随容差发生变化。这里不能宣布它超过 VFCR 的数学精度极限。GNN 与 VFCR 的纠错增益、速度预算不同；本轮实现的是“同信息、同初值、同题目”的比较，而非相同物理预算的全面公平比较。\n\n');
fprintf(f,'上一轮 M8a 与 VFCR-S 的主要差距来自一阶差分。本轮不再只有“减小一阶 h”这个办法，BDF2 已消除了大部分差分偏差。这是本轮最确定的改进。\n\n');
fprintf(f,'## 2. Delta study：原题精度提高，完整状态变得更难跟踪\n\n');
fprintf(f,'保持 BDF2 的 h=1e-5，分别改变 delta。x 是相对原始 QP 最优解的 RMSE，bias 是扰动根与原题解之间的 RMSE。\n\n');
print(f,z([3,6,7,8:13],{'model','delta','xd','x','bias','gpeak'}));
fprintf(f,'上表是原始固定网格统计。原题误差在其他窄切换处仍会受网格影响，因此又用自适应积分专门复核，最终应使用下面的原题 RMSE：\n\n');
print(f,adaptive);
fprintf(f,'biasbound/xbound 是平方误差时间积分的估计误差，不是 RMSE 本身的误差界。\n\n');
fprintf(f,'自适应积分容差从 1e-6 收紧到 1e-8，并改变初始分段，x RMSE 最大相对变化 %.3g%%。M8a-BDF2 的原题 RMSE 从 %.6g 降到 %.6g，改善约 **%.1f 倍**；xd 始终很小，说明此时主要限制来自 PFB 扰动，而不是追不上扰动根。三点扫描不足以声称普适的 delta 幂律。\n\n', ...
    100*max(abs(adaptive.x./prior.x-1)),adaptive.x(1),adaptive.x(3),adaptive.x(1)/adaptive.x(3));
fprintf(f,'但完整状态误差尖峰增大。参考根的病态程度与速度也明显恶化：\n\n');
print(f,check(:,{'delta','sigma','speed','residual','polished','gshift','xshift'}));
fprintf(f,'sigma 是参考轨迹上最小奇异值的最小值；speed 是完整状态最大速度。gshift/xshift 表示用稳定残差进行最多三次 Newton 校正时参考解的变化量，用来判断参考解的浮点误差；不是模型误差。小 delta 下参考乘子的误差更值得检查，不能将 gpeak 的全部变化都简单解释成跟踪误差。\n\n');
fprintf(f,'![delta 比较](figures/delta.png)\n\n');
fprintf(f,'图中的原题 RMSE 使用自适应积分，其余曲线来自保存网格。对检测到的尖峰进一步局部加密至 5e-9 s，并使用校正后的参考根，结果如下：\n\n');
print(f,peak);
fprintf(f,'delta=1e-6 的 M8a-BDF2 完整状态峰值应理解为约 6e-8 的量级；最后几位仍受参考根和浮点误差影响。局部加密不改变小 delta 更敏感的结论。\n\n');
fprintf(f,'## 3. Noise study：误差由观测噪声主导\n\n');
fprintf(f,'本轮注入残差测量噪声 xi_obs=xi+n，频率 20–200 Hz，16 个频率叠加，固定三个高斯系数种子，振幅 A=1e-5。各模型看到同一条噪声，噪声进入纠错、时间差分和 VFCR 积分，评价目标仍为干净根。不是在每次 RHS 上重新抽随机数，也不是原 VFCR 的 RHS 外部扰动噪声。\n\n');
noise = table();
for model = ["M8a-FD1","M8a-BDF2","M8-reg","VFCR-FD"]
    part = z(z.group=="noise" & z.model==model,:);
    row = table(model,mean(part.xd),median(part.xd),min(part.xd),max(part.xd), ...
        min(part.gpeak),max(part.gpeak),mean(part.residual), ...
        'VariableNames',{'model','meanxd','medianxd','minxd','maxxd','mingpeak','maxgpeak','meanresidual'});
    noise = [noise;row];
end
print(f,noise); writetable(noise,fullfile(folder,'noise_summary.csv'));
fprintf(f,'一阶、二阶 GNN 和 VFCR-FD2 的平均 x 误差非常接近。无噪声时 BDF2 的极小截断误差，在这个噪声水平下被观测误差淹没。这是三个种子的探索性结果，不是大样本统计显著性结论。\n\n');
fprintf(f,'![带噪比较](figures/noise.png)\n\n');
fprintf(f,'### 固定正则化的代价\n\n');
fprintf(f,'M8-reg 仅将预测阻尼设为 1e-6，纠错与 M8a 相同。它的无噪声基线：\n\n');
print(f,e(7,{'model','xd','gpeak','respeak'}));
fprintf(f,'其带噪 x RMSE 没有明显改善，完整状态尖峰却较大。无噪声时就存在较大的预测阻尼误差，说明它的坏结果不能全部归咎于噪声。固定阻尼会一并削弱目标运动和噪声变化，不能自动得到好的抗噪模型。\n\n');
fprintf(f,'### 振幅、差分步长与常量噪声\n\n');
print(f,z([15:17,26:37],{'model','noise','amp','h','xd','gpeak','residual','obs'}));
fprintf(f,'residual 是干净残差 RMS；obs 是带噪观测残差 RMS。观测残差小，不代表原始问题残差小：可能只是跟着带噪根运动。常量噪声的差分为零，但仍使纠错目标偏移，所以“差分不放大常量噪声”不等于“常量噪声没有影响”。\n\n');
fprintf(f,'本轮有限带宽噪声在 h 很小时接近自身的时间导数，不呈现独立样本那种无界 1/h 放大；改变 h 的实际 ODE 结果只能解释本轮噪声带宽。另对独立高斯样本做代数测试，得到 BDF2 导数噪声约为一阶的 **1.803 倍**，且两者都随 1/h 增长，记录见[iid.csv](iid.csv)。该测试不等于完成了白噪声 ODE 实验。\n\n');
fprintf(f,'![噪声敏感性](figures/sensitivity.png)\n\n');
fprintf(f,'## 4. Delta 与噪声联合检查\n\n');
fprintf(f,'固定种子 1、A=1e-5、h=1e-5，仅将 delta 从 1e-4 改为 1e-6：\n\n');
base = z([15,17],{'model','delta','xd','x','gpeak','residual','obs'});
print(f,[base;j(:,base.Properties.VariableNames)]);
fprintf(f,'减小 delta 仍降低原题平均误差，但完整状态尖峰明显放大；M8a-BDF2 的 gpeak 放大约 %.1f 倍。此时 x 跟踪误差已经与 delta 的原题偏差相当，小 delta 不再能单独决定最终精度。这里只做同一种子的机制检查，不代表所有噪声分布。\n\n',j.gpeak(1)/z.gpeak(15));
fprintf(f,'## 5. 验证和可复现性\n\n');
audit = readtable(fullfile(folder,'audit.csv')); print(f,audit);
fprintf(f,'公式审计 100 个随机状态；解析 VFCR 的稳定实现与原 my_system.m 等价，最大相对差约 %.3g。因果启动和重复调用已检查。所有 value-only 动力学运行时删除六个解析导数字段。\n\n',audit.vfcr);
fprintf(f,'37 组主轨迹均重新读取并独立计算时间加权 RMSE；deval 与保存轨迹差为零。加密输出网格后，扰动根 x 跟踪 RMSE 最大相对变化 %.3g%%，完整状态峰值最大相对变化 %.3g%%。详见[verification.csv](verification.csv)，其中 xdp、gpeakp 使用 Newton 校正后的参考根，专门检查原参考误差对排名的影响。\n\n',100*max(v.rmsechange),100*max(v.peakchange));
fprintf(f,'参考根校正后，主配置 M8a-BDF2 的加密网格 x RMSE 为 %.6g、完整状态峰值 %.6g；因此最后几位数字应按数值量级理解。\n\n',v.xdp(3),v.gpeakp(3));
fprintf(f,'额外 8 组复核和 2 组联合检查也重新读取并复算，轨迹、误差数组和抽查的干净残差一致，见[other_check.csv](other_check.csv)。源码与共用依赖的校验值见[sources.csv](sources.csv)。\n\n');
fprintf(f,'噪声组另收紧容差、减半最大积分步长，仍然完成：\n\n');
print(f,e(3:6,{'model','rtol','atol','h','xd','gpeak','calls'}));
fprintf(f,'该表前两组收紧容差，后两组把最大积分步长减为 0.00025 s。四幅图均重新打开 FIG 核对曲线数据，并逐幅目视检查。\n\n');
fprintf(f,'运行秒数受同时进行的 MATLAB 检查和每次求值成本影响，不能用来宣布普适速度优势。每组原始配置、轨迹、误差和预算结果都保留。调试时的空采样保存错误保存在 data；修复后完整主数据在 results。\n\n');
fprintf(f,'## 6. 下一步建议\n\n');
fprintf(f,'1. **无噪声主模型采用 M8a-BDF2，先用 h=1e-5。** 这一步已有明确实验支持；本题继续缩小 h 没有收益。\n');
fprintf(f,'2. **无噪声高原题精度用 delta=1e-6，带噪时一起检查乘子尖峰。** 本轮噪声下先以 delta=1e-4 作为稳定研究基线，而不是只看 x 的平均误差。\n');
fprintf(f,'3. **抗噪下一步先试时间滤波或固定采样的导数估计。** 从带噪历史残差中估计时间偏导，并单独比较滤波误差、延迟与尖峰。不要同时改变太多参数。该改进本轮尚未实现。\n');
fprintf(f,'4. **保留 VFCR-FD2 作为 value-only 对照。** 如果要比较硬件控制能力，应另外匹配实际物理速度预算；如果要研究采样噪声，应固定采样率和历史插值，再做多种子闭环实验。\n');
fprintf(f,'5. **暂不加入固定时间激活或 M7 缩放。** 先处理观测噪声和小 delta 的病态放大，当前数据没有证明一般鲁棒性或固定时间收敛。\n\n');
fprintf(f,'[完整主指标](results/metrics.csv) · [高精度与噪声复核](extra/metrics.csv) · [联合检查](joint/metrics.csv) · [系统说明](系统说明.md)\n');
end
function print(f,t)
% 将实际数值按六位有效数字写入表格。
names = t.Properties.VariableNames;
fprintf(f,'| %s |\n',strjoin(names,' | ')); fprintf(f,'| %s |\n',strjoin(repmat({'---'},1,numel(names)),' | '));
for k = 1:height(t)
    values = strings(1,width(t));
    for m = 1:width(t)
        a = t{k,m};
        if isnumeric(a), values(m)=string(sprintf('%.6g',a)); else, values(m)=string(a); end
    end
    fprintf(f,'| %s |\n',strjoin(values,' | '));
end
fprintf(f,'\n');
end
