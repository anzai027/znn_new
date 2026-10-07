function report(stats,results,c,speeds,audit,folder)
% 把改善和不足分别写入中文报告。
fid = fopen(fullfile(folder,'结论.md'),'w','n','UTF-8');
guard = onCleanup(@() fclose(fid));
fprintf(fid,'# 第二版 GNN：建议能否改善性能\n\n');
fprintf(fid,'共尝试 %d 组；完整返回 %d 组，中断 %d 组；完整但速度检查异常 %d 组。M1 原无噪声沿用保存的异常轨迹，不能作为性能优势证据。\n\n', ...
    height(stats),sum(stats.complete),sum(~stats.complete),sum(stats.complete & stats.motion == 0));
fprintf(fid,'## 先回答改进是否有效\n\n');
fprintf(fid,'M2 是旧平滑梯度方向，M3 在此基础上改善方向，M4 再把变量分块分配速度；它们的公式见 模型说明.md。下面分别记录平滑版和不平滑版的结果，不用中断时很短的运行时间冒充速度优势。\n\n');
a = find(stats.group == "scan" & stats.model == "M2" & stats.gamma == 200 & stats.eps == 1e-2,1);
b = find(stats.group == "scan" & stats.model == "M3" & stats.gamma == 200 & stats.eps == 1e-2 & stats.lambda == 1e-6,1);
d = find(stats.group == "scan" & stats.model == "M4" & abs(stats.gamma-200) < 1e-8 & stats.eps == 1e-2 & stats.lambda == 1e-6,1);
fprintf(fid,'在 gamma=200、eps=1e-2、同初值同题目的无噪声比较中，M2 的 delta 解 RMSE 为 %.6g，M3 为 %.6g，误差比为 %.4g；对应 RHS 调用数 %d/%d。方向改进的效果可以与单纯提高增益分开观察。\n\n', ...
    stats.rmsed(a),stats.rmsed(b),stats.rmsed(a)/stats.rmsed(b),stats.calls(a),stats.calls(b));
fprintf(fid,'相同总速度上限 200 时，分块 M4 的 delta 解 RMSE 为 %.6g；不能预先认定分块更好。主 M4 的 rates=[100,200,200] 总速度上限为 300，不能把它与 gamma=200 的差异全部归因于分块。\n\n',stats.rmsed(d));
fprintf(fid,'本次证据支持继续使用 M3 的预条件方向；所选 M4 分块设置没有带来额外精度优势。无噪声时 VFCR-S 仍比主 M3 更准确；共同状态噪声下另看逐种结果，不能把某一噪声入口的优势扩大到所有抗噪条件。\n\n');
fprintf(fid,'主设置预先固定：M2/M3 gamma=200，eps=1e-2；M3/M4 lambda=1e-6，M4 rates=[100,200,200]。主实验没有依据扫描结果择优重设参数，全部扫描参数保留。\n\n');
fprintf(fid,'M3/M4 在零误差处仍给出零控制速度；本次只能检查是否改善所测跟踪精度和运行表现，不能证明严格零误差移动根跟踪、有限时间或固定时间性质。\n\n');
fprintf(fid,'## 统一条件与指标\n\n');
fprintf(fid,'MATLAB R2025b、ode15s、0～10 秒、输出间隔 0.005 秒、RelTol=1e-5、AbsTol=1e-7、MaxStep=0.02，每组限 180 秒或 50 万调用。随机噪声每 0.01 秒保持并分段积分，两个随机分布各用种子 1～5。\n\n');
fprintf(fid,'VFCR 直接调用原 my_system，论文参数保持不变；VFCR-S 通过已有 problem.phi 入口平滑激活函数。主比较统一在状态导数加入相同噪声，不给 VFCR 积分状态额外加噪声；不混用论文内部噪声入口。\n\n');
fprintf(fid,'误差和约束统计使用 5～10 秒；Tc 表示残差最后持续低于阈值且至少观察 0.5 秒。完整只表示积分返回有限数，速度检查也只是必要条件。NaN 是缺失，Inf 是未持续达标。\n\n');
fprintf(fid,'真解为独立一维区间二次规划解，delta 解由障碍驻点方程独立求出；不会用某一神经模型当参考解。\n\n');
fprintf(fid,'任务结构中的 lambda 只有 M3/M4 使用，gamma/rates 只有 GNN 系列使用；VFCR 实际传入的参数是 config 中的 r1=r2=a=lambda1=lambda2=1、p=q=0.5、delta=1e-4，统一任务结构的占位字段不改变 VFCR。\n\n');
fprintf(fid,'## 全部机制扫描\n\n');
fprintf(fid,'| 模型 | gamma 或总上限 | rates | eps | lambda | 完整 | delta RMSE | 完整状态 RMSE | 最大残差 | 秒 | 调用数 |\n');
fprintf(fid,'|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for id = find(stats.group == "scan").'
    s = stats(id,:);
    fprintf(fid,'| %s | %.4g | %s | %.0e | %.0e | %d | %.5g | %.5g | %.5g | %.4g | %d |\n', ...
        s.model,s.gamma,s.rates,s.eps,s.lambda,s.complete,s.rmsed,s.full,s.residual,s.seconds,s.calls);
end
fprintf(fid,'\n## 统一噪声下的主实验\n\n');
fprintf(fid,'随机噪声指标为完整且通过检查的五次结果中位数；失败不混入指标。\n\n');
fprintf(fid,'| 噪声 | 模型 | 完整且通过检查 | 真解 RMSE | delta RMSE | delta 最大误差 | 最大残差 | 等式违反 | 不等式违反 | Tc 1e-6 | Tc 1e-4 |\n');
fprintf(fid,'|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for kind = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"]
    for model = ["M2","M3","M4","VFCR","VFCR-S"]
        s = stats(stats.group == "main" & stats.noise == kind & stats.model == model & stats.eps > 0,:);
        good = s.complete & (isnan(s.motion) | s.motion == 1); v = s(good,:);
        fprintf(fid,'| %s | %s | %d/%d | %.5g | %.5g | %.5g | %.5g | %.5g | %.5g | %.5g | %.5g |\n', ...
            kind,model,sum(good),height(s),median(v.rmse),median(v.rmsed),median(v.peakd), ...
            median(v.residual),median(v.equality),median(v.inequality),median(v.tc6),median(v.tc4));
    end
end
fprintf(fid,'\n## 三次计时\n\n');
fprintf(fid,'只使用完整且通过必要检查的本次运行；用时表示计算十秒轨迹的墙钟耗时，不是数学收敛时间。\n\n');
fprintf(fid,'| 模型 | 噪声 | 可用次数 | 秒中位数 | RHS 次数中位数 |\n');
fprintf(fid,'|---|---|---:|---:|---:|\n');
for kind = ["zero","cosine"]
    for model = ["M2","M3","M4","VFCR","VFCR-S"]
        s = stats(ismember(stats.group,["main","timing"]) & stats.model == model ...
            & stats.noise == kind & stats.eps > 0,:);
        good = s.complete & (isnan(s.motion) | s.motion == 1); v = s(good,:);
        fprintf(fid,'| %s | %s | %d/3 | %.5g | %.5g |\n',model,kind,sum(good),median(v.seconds),median(v.calls));
    end
end
fprintf(fid,'\n## 严格容差检查\n\n');
for model = ["M2","M3","M4","VFCR","VFCR-S"]
    i = find(stats.group == "main" & stats.model == model & stats.noise == "zero",1);
    j = find(stats.group == "precision" & stats.model == model,1);
    a = results{stats.id(i)}; b = results{stats.id(j)}; gap = NaN;
    if a.complete && b.complete, gap = max(vecnorm(a.y(:,1:2)-b.y(:,1:2),2,2)); end
    fprintf(fid,'%s：普通/严格容差完整情况 %d/%d，delta RMSE %.6g/%.6g，两个数值轨迹最大距离 %.6g。\n\n', ...
        model,a.complete,b.complete,a.metrics.rmsed,b.metrics.rmsed,gap);
end
for model = ["M3","M4"]
    i = find(stats.group == "main" & stats.model == model & stats.noise == "zero" & stats.eps == 0,1);
    j = find(stats.group == "precision" & stats.model == model & stats.eps == 0,1);
    a = results{stats.id(i)}; b = results{stats.id(j)};
    fprintf(fid,'%s 不平滑版：普通/严格容差完整情况 %d/%d，最后积分时间 %.6g/%.6g 秒；未取得完整轨迹，不能据此计算稳态精度差异。\n\n', ...
        model,a.complete,b.complete,a.last,b.last);
end
fprintf(fid,'## 移动参考解速度复查\n\n');
fprintf(fid,'| 差分网格 | g 最大速度 | x 最大速度 | mu1 最大速度 | mu2 最大速度 | 最大割线速度 |\n');
fprintf(fid,'|---:|---:|---:|---:|---:|---:|\n');
for k = 1:3
    fprintf(fid,'| %.4g | %.6g | %.6g | %.6g | %.6g | %.6g |\n',speeds.steps(k),speeds.peaks(k,:),speeds.secant(k));
end
fprintf(fid,'\n网格细化后，完整参考 g 的速度估计继续上升，不能把粗网格的 141.1 当作真实峰值上界。最细网格的最大割线速度为 %.6g，超过本次 100、200、300 的控制速度上限；在无噪声下，这足以排除这些限速模型在整个区间对完整 g 精确跟踪。它不直接证明单独 x 也受到相同限制；差分峰值仍不是严格全局导数上界。\n\n',speeds.secant(end));
fprintf(fid,'## Example 1 初值和边界诊断\n\n');
fprintf(fid,'reduced 删除了原本不存在的不等式，GNN 状态为三维，VFCR 为六维；large 保留 1e8 写法，状态仍为七维/十四维。两者数学真实最优 x 相同，但扰动参考 g 和状态空间不同。\n\n');
fprintf(fid,'| 设置 | 模型 | eps | 初值范围 | 完整 | 最后时间 | delta RMSE | 最大残差 | Tc 1e-4 | 失败原因 |\n');
fprintf(fid,'|---|---|---:|---:|---:|---:|---:|---:|---:|---|\n');
for id = find(ismember(stats.group,["initial","bounds"])).'
    s = stats(id,:);
    fprintf(fid,'| %s | %s | %.0e | %.4g | %d | %.5g | %.5g | %.5g | %.5g | %s |\n', ...
        s.variant,s.model,s.eps,s.scale,s.complete,s.last,s.rmsed,s.residual,s.tc4,s.reason);
end
fprintf(fid,'\n## 最终判断的边界\n\n');
fprintf(fid,'预条件方向能否改善，优先看相同 gamma、eps 的 M2/M3 对照；分块能否改善，优先看相同总速度上限的 M3/M4 对照。单次扫描较短耗时不足以替代三次主实验计时。\n\n');
fprintf(fid,'如果新模型仍未达到持续 1e-4/1e-6 阈值，就只能说显著降低误差，不能说实现了原预期的高精度或零误差求解。较小 x 误差也不能掩盖完整 g、残差或约束违反量。\n\n');
fprintf(fid,'结构方面保留无系数导数、无积分状态；新增阻尼最小二乘求解，原来“不用线性方程求解”的优点已不再保留。\n\n');
fprintf(fid,'Example 1 的大初值不能用固定的 5～10 秒窗口直接解释成已经进入稳态；M3 初值范围 1000 时该窗口仍含明显过渡过程，必须保留这一未满足表现。固定速度上限还会使大初值所需时间增长，不能因此声称全局固定时间收敛。\n\n');
fprintf(fid,'平滑纯反馈、不同初值曲线和噪声样本都不能证明严格有限时间、固定时间或任意噪声鲁棒性。自适应增益和预测项属于建议的后续阶段，本轮未混入这些模块。\n\n');
fprintf(fid,'## 不平滑版的对应测试\n\n');
fprintf(fid,'两种版本均已测试；M3/M4 的不平滑版在全部主噪声条件、三次计时、严格容差和不同初值实验中均中断。平滑版在 15 组主条件中均能跑完，但仍未满足持续残差 1e-4 的要求。因此，平滑缓解了当前积分方式的失败，并未完成原定的高精度和有限/固定时间目标。\n\n');
i = find(stats.group == "main" & stats.model == "M3" & stats.eps == 0 & stats.noise == "zero",1);
fprintf(fid,'不平滑 M3 的无噪声残差曾在最后接受时间 %.8g 秒降到 %.6g，随后积分中断。它已经非常接近瞬时目标，不能说从未接近目标；但缺少后续轨迹，不能把这一个瞬间替代持续收敛判据。\n\n', ...
    stats.last(i),stats.lastresidual(i));
fprintf(fid,'M3/M4 的 eps=0 为不平滑版，没有另加分母小量；方向恰为零的块取零速度。lambda 仍是方向阻尼，delta 仍是 PFB 参数，二者不是本次比较的梯度归一化平滑尺度。普通 ode15s 的失败不等同于微分包含意义下的数学发散。\n\n');
fprintf(fid,'| 模型 | eps | 主实验次数 | 完整且通过速度检查 | 中断 | 完整但速度异常 |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|\n');
for model = ["M3","M4"]
    for eps = [1e-2,0]
        s = stats(stats.group == "main" & stats.model == model & stats.eps == eps,:);
        fprintf(fid,'| %s | %.0e | %d | %d | %d | %d |\n',model,eps,height(s), ...
            sum(s.complete & s.motion == 1),sum(~s.complete),sum(s.complete & s.motion == 0));
    end
end
fprintf(fid,'\n| 噪声 | 模型 | eps | 可用次数 | delta RMSE | 最大残差 | 等式违反 | 不等式违反 |\n');
fprintf(fid,'|---|---|---:|---:|---:|---:|---:|---:|\n');
for kind = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"]
    for model = ["M3","M4"]
        for eps = [1e-2,0]
            s = stats(stats.group == "main" & stats.noise == kind & stats.model == model & stats.eps == eps,:);
            good = s.complete & s.motion == 1; v = s(good,:);
            fprintf(fid,'| %s | %s | %.0e | %d/%d | %.5g | %.5g | %.5g | %.5g |\n', ...
                kind,model,eps,sum(good),height(s),median(v.rmsed),median(v.residual),median(v.equality),median(v.inequality));
        end
    end
end
fprintf(fid,'\n| 模型 | eps | 噪声 | 可用次数 | 三次用时中位数 | 调用次数中位数 |\n');
fprintf(fid,'|---|---:|---|---:|---:|---:|\n');
for model = ["M3","M4"]
    for eps = [1e-2,0]
        for kind = ["zero","cosine"]
            s = stats(ismember(stats.group,["main","timing"]) & stats.model == model ...
                & stats.eps == eps & stats.noise == kind,:);
            good = s.complete & s.motion == 1; v = s(good,:);
            fprintf(fid,'| %s | %.0e | %s | %d/3 | %.5g | %.5g |\n', ...
                model,eps,kind,sum(good),median(v.seconds),median(v.calls));
        end
    end
end
fprintf(fid,'\n不平滑版若没有完整且通过检查的轨迹，其误差和耗时不参与性能排名；某些通过速度检查的完整轨迹仍可能误差很大，不能等同于准确求解。\n\n');
fprintf(fid,'全部 %d 对 FIG/PNG 已重开读取，噪声与初值一致性、失败后的缺失值、原文件校验均通过。\n',numel(audit.names));
end
