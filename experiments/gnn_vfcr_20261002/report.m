function report(s,refs,c,meta,audit,folder)
% 根据实际结果生成中文结论。
fid = fopen(fullfile(folder,'结论.md'),'w','n','UTF-8');
assert(fid >= 0,'不能保存结论。');
guard = onCleanup(@() fclose(fid));
fprintf(fid,'# 新 GNN 与原 VFCR-ZNN：实验结果\n\n');
fprintf(fid,'实验共尝试 %d 组，完整计算到 10 秒的有 %d 组，其余 %d 组如实记录为中断。\n\n', ...
    height(s),sum(s.complete),sum(~s.complete));
fprintf(fid,'“完整”只表示返回了有限数。完整轨迹中还有 %d 组未通过速度上限检查，其数值只作为异常诊断，不用于宣称性能优势；其余完整轨迹也仍需看精度和约束。\n\n', ...
    sum(s.complete & s.tested & s.motion == 0));
main = s(s.group == "main",:);
gnn = main(main.model == "GNN",:);
vfcr = main(main.model == "VFCR-ZNN",:);
fprintf(fid,'## 先看结论\n\n');
fprintf(fid,'在 Example 2 的统一状态噪声主实验中，新 GNN 返回完整轨迹 %d/%d 组，其中 %d 组速度异常；原 VFCR-ZNN 返回完整轨迹 %d/%d 组。', ...
    sum(gnn.complete),height(gnn),sum(gnn.complete & gnn.motion == 0),sum(vfcr.complete),height(vfcr));
if ~all(gnn.complete)
    fprintf(fid,'因此当前原式、gamma=100、p=1 的 MATLAB 实现尚不能被本次实验确认是可靠完成这些 TVQP 测试的求解器。');
else
    fprintf(fid,'新模型完成了所测轨迹，但是否更准确仍需看下面的误差和约束指标。');
end
if ~all(vfcr.complete)
    fprintf(fid,'原 VFCR 也有未完成的运行，必须把两者的数值实现问题都计入判断。');
end
fprintf(fid,'求解器中断不能直接证明数学模型不存在收敛解；同样，中断前误差很小也不能证明之后能够一直跟踪。\n\n');
fprintf(fid,'**求解是否准确：** 主参数 gamma=100、p=1 尚未通过。部分结果数值异常，随机噪声下没有完整结果；因此不能用这些输出证明其已经保留 VFCR 的 TVQP 求解能力。原 VFCR 无噪声也受调用上限影响，不能把其局部小误差当成完整通过。\n\n');
fprintf(fid,'**哪个更快：** 目前没有可信的主参数 GNN 用时优势证据。无噪声和余弦噪声的三次 GNN 完整输出均存在速度异常，其短耗时不能表示求解更快。\n\n');
fprintf(fid,'**哪个更抗噪：** 在本次共同状态噪声、ode15s 和主参数下，VFCR 的运行可靠性更好；Gamma 与高斯共十组均跑完整，新 GNN 十组全部中断。该判断针对本次数值实现与样本，不是对任意噪声的数学保证。\n\n');
fprintf(fid,'**结构是否更简单：** 已经得到代码和调用检查支持。新 GNN 使用 7 个状态，VFCR 使用 14 个状态；新 GNN 没有积分状态，不读取六个系数导数字段，也不执行 J\\右端项。它仍要计算 J 和 J^T*xi。原 VFCR 文件校验未变。\n\n');
fprintf(fid,'当前单幂次新模型在一般意义上不是固定时间模型；有限时间与抗噪跟踪仍需额外条件及零点处的严谨解释，不能只靠这次曲线作保证。\n\n');

fprintf(fid,'## 主实验误差与抗噪性\n\n');
fprintf(fid,'以下均是两个模型共用状态噪声和初值的结果，新 GNN 固定 gamma=100、p=1。RMSE 与最大误差使用 5～10 秒的数据；随机噪声统计五个种子中完整轨迹的逐组指标中位数，并列出完整次数。速度异常行的数值仅用于诊断，不能解释为数学模型真实误差。NA 表示没有完整运行，Inf 表示完整输出未满足持续达标判据。\n\n');
fprintf(fid,'| 噪声 | 模型 | 完整次数 | 速度检查 | 真解 RMSE | delta 解 RMSE | delta 解最大误差 | 最大残差 | Tc(1e-6) | Tc(1e-4) |\n');
fprintf(fid,'|---|---|---:|---|---:|---:|---:|---:|---:|---:|\n');
for kind = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"]
    for model = ["VFCR-ZNN","GNN"]
        rows = main(main.noise == kind & main.model == model,:);
        good = rows(rows.complete,:);
        fprintf(fid,'| %s | %s | %d/%d | %s | %s | %s | %s | %s | %s | %s |\n', ...
            label0(kind),model,height(good),height(rows),check0(good),number0(median0(good.rmse)), ...
            number0(median0(good.rmsed)),number0(median0(good.peakd)), ...
            number0(median0(good.residual)),number0(median0(good.tc6)),number0(median0(good.tc4)));
    end
end
fprintf(fid,'\n| 噪声 | 模型 | 最大等式违反量 | 最大不等式违反量 |\n|---|---|---:|---:|\n');
for kind = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"]
    for model = ["VFCR-ZNN","GNN"]
        good = main(main.noise == kind & main.model == model & main.complete,:);
        fprintf(fid,'| %s | %s | %s | %s |\n',label0(kind),model, ...
            number0(median0(good.equality)),number0(median0(good.inequality)));
    end
end
fprintf(fid,'\n真正解与 delta 近似解的最大距离：Example 1 为 %.6e，Example 2 为 %.6e。这项偏差来自 delta，不等同于动力学的跟踪误差。\n\n', ...
    max(vecnorm(refs{1}.x-refs{1}.xd,2,2)),max(vecnorm(refs{2}.x-refs{2}.xd,2,2)));
fprintf(fid,'目标值不能单独证明求解正确：不满足约束的点可能拥有更低的目标值。完整 CSV 另外记录等式、不等式违反量与目标差。\n\n');

fprintf(fid,'## 论文原噪声入口的 VFCR 复现\n\n');
fprintf(fid,'VFCR 的论文噪声加在 J*g_dot 的右端；主实验的共同噪声直接加在 g_dot 上。两种入口的结果单列，不混用来宣布谁更抗噪。\n\n');
fprintf(fid,'| 噪声 | 完整次数 | delta 解 RMSE | 最大残差 |\n|---|---:|---:|---:|\n');
for kind = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"]
    rows = s(s.group == "paper" & s.noise == kind,:);
    good = rows(rows.complete,:);
    fprintf(fid,'| %s | %d/%d | %s | %s |\n',label0(kind),height(good),height(rows), ...
        number0(median0(good.rmsed)),number0(median0(good.residual)));
end

fprintf(fid,'\n## 不同初值与参数\n\n');
fprintf(fid,'| 初值范围 | 模型 | 状态 | 实际到达时间 | 最后残差 | Tc(1e-6) |\n|---|---|---|---:|---:|---:|\n');
rows = s(s.group == "initial",:);
for k = 1:height(rows)
    r = rows(k,:);
    fprintf(fid,'| [0,%g] | %s | %s | %.8g | %.3e | %s |\n',r.scale,r.model, ...
        r.status,r.last,r.lastresidual,number0(r.tc6));
end
fprintf(fid,'\n新模型补充参数如下，其他设置保持相同；p 是新模型分母的指数，VFCR 的 p=0.5 是另一种激活函数参数。\n\n');
fprintf(fid,'| gamma | p | 噪声 | 状态 | 速度检查 | 实际到达时间 | delta 解 RMSE | 最大残差 |\n|---:|---:|---|---|---|---:|---:|---:|\n');
rows = s(s.group == "parameter",:);
for k = 1:height(rows)
    r = rows(k,:);
    fprintf(fid,'| %g | %g | %s | %s | %s | %.8g | %s | %s |\n',r.gamma,r.p,label0(r.noise), ...
        r.status,check0(r),r.last,number0(r.rmsed),number0(r.residual));
end
fprintf(fid,'\n有参数版本成功并不等于主参数版本已经通过；较小误差也可能是增益变大造成的，因此结论必须注明参数。\n\n');
candidate = s(s.group == "parameter" & s.gamma == 100 & s.p == 0.5,:);
zero = candidate(candidate.noise == "zero",:);
cosine = candidate(candidate.noise == "cosine",:);
baseline = main(main.model == "VFCR-ZNN" & main.noise == "cosine",:);
fprintf(fid,'gamma=100、p=0.5 在无噪声与共同余弦噪声下均跑完整。相对 delta 解的 RMSE 分别为 %.6g 和 %.6g；余弦下原 VFCR 为 %.6g，误差更小。新模型余弦下最大等式违反量 %.6g，原 VFCR 为 %.6g。p=0.5 因此有近似跟踪的实验支持，但没有达到 1e-6 的残差要求，也没有接受更严格容差或全部随机噪声测试。无噪声 VFCR 未跑完整，不能作完整区间排名。\n\n', ...
    zero.rmsed,cosine.rmsed,baseline.rmsed,cosine.equality,baseline.equality);

fprintf(fid,'## 计算成本\n\n');
fprintf(fid,'下表统计主实验和用时重复实验的三次运行。完整次数仍全部列出，用时中位数排除未完成与速度检查失败的轨迹；未适用速度检查不等于已证明其数值正确。排除的原始耗时保留在 CSV 中。精度不同的完整结果也不能直接解释为求解速度排名。\n\n');
fprintf(fid,'| 噪声 | 模型 | 完整次数 | 可统计次数 | 用时中位数 (s) | 调用次数中位数 |\n|---|---|---:|---:|---:|---:|\n');
for kind = ["zero","cosine"]
    for model = ["VFCR-ZNN","GNN"]
        take = any(s.group == ["main","timing"],2) & s.noise == kind & s.model == model & s.seed == 1;
        rows = s(take,:);
        good = rows(rows.complete & ~(rows.tested & rows.motion == 0),:);
        fprintf(fid,'| %s | %s | %d/%d | %d/%d | %s | %s |\n',label0(kind),model,sum(rows.complete), ...
            height(rows),height(good),height(rows),number0(median0(good.seconds)),number0(median0(good.calls)));
    end
end
fprintf(fid,'\n以上是 ODE 求解的实际墙钟时间，不含参考解、指标计算、保存和绘图；机器负载可能影响数值，结构更简单不自动意味着运行更快。\n\n');

fprintf(fid,'## 精度与数值失败\n\n');
fprintf(fid,'p=1 时，原式在非零梯度处满足 ||g_dot||=gamma，在零梯度处程序取零。因此，在噪声范数不超过 B 时，真实轨迹必须满足 ||g(t)-g(0)|| <= (gamma+B)*t，相邻点的距离也不能超过 (gamma+B) 乘以时间差。\n\n');
fprintf(fid,'这里对接受点和输出网格同时检查，并额外允许一百倍配置容差对应的误差，避免把微小数值误差认成严重异常。此检查只用于排除明显不可能的输出，通过它仍不等于求解准确。\n\n');
fprintf(fid,'| 速度异常的运行编号 | 参数 | 噪声 | 最大超限距离 | 扣除容差后超限 |\n|---:|---|---|---:|---:|\n');
rows = s(s.tested & s.motion == 0,:);
for k = 1:height(rows)
    r = rows(k,:);
    fprintf(fid,'| %d | gamma=%g, p=%g | %s | %.6e | %.6e |\n', ...
        r.id,r.gamma,r.p,label0(r.noise),r.excess,r.margin);
end
fprintf(fid,'\n主参数无噪声和余弦两组的严重超限在求解器实际接受点中也存在，并不只是画图或插值的现象。其他超限值见 MAT 内分项诊断。因此，大误差不能直接归因于数学模型必然发散，短耗时也不能证明模型求解更快。\n\n');
fprintf(fid,'| 模型 | RelTol | AbsTol | 状态 | 实际到达时间 | 最后残差 | delta 解 RMSE |\n|---|---:|---:|---|---:|---:|---:|\n');
rows = s((s.group == "main" & s.noise == "zero") | s.group == "precision",:);
for k = 1:height(rows)
    r = rows(k,:);
    fprintf(fid,'| %s | %.0e | %.0e | %s | %.8g | %.3e | %s |\n',r.model,r.rtol,r.atol, ...
        r.status,r.last,r.lastresidual,number0(r.rmsed));
end
fprintf(fid,'\n每组允许最多 %.0f 秒或 %d 次函数调用。达到限额、求解器步长不足、非有限导数等均记录为未完成；之后的状态保存为 NaN。中断运行的 Tc 和后半段指标也为 NaN。\n\n',c.seconds,c.calls);
fprintf(fid,'| 中断类型 | 次数 |\n|---|---:|\n');
ids = unique(s.identifier(~s.complete));
for id = ids.'
    fprintf(fid,'| %s | %d |\n',replace(id,'|','/'),sum(~s.complete & s.identifier == id));
end
fprintf(fid,'\n小残差附近的步长、梯度范数和雅可比最小奇异值均已保存；其中最小奇异值只是采样诊断，不能被当成所有状态上的统一下界。没有添加 epsilon，也没有提前冻结或停止跟踪。\n\n');
fprintf(fid,'Example 1 使用论文的 1e8 边界，PFB 项的大数相减会受到浮点精度影响；该设置没有改动。\n\n');

fprintf(fid,'## 为什么当前模型不是一般的固定时间模型\n\n');
fprintf(fid,'取最简单的静态问题 min(x^2/2)，没有约束。此时 xi=x、J=1，新模型成为 x_dot=-gamma*sign(x)*abs(x)^(1-p)。它的到达时间为：\n\n');
fprintf(fid,'T = abs(x(0))^p / (gamma*p)。\n\n');
fprintf(fid,'只要 p>0，初值越来越大时，这个时间就没有统一上限。以 gamma=100、p=1 为例，初值 1、10、100、1000 的到达时间分别是 0.01、0.1、1、10 秒。这给出了当前通用单幂次模型不具备全局固定时间性质的明确反例。\n\n');
fprintf(fid,'静态有限时间需要合适的雅可比下界；时变问题还需要控制问题变化速度，并解释 p=1 在零点的不连续动力学。当前程序在零点取零速度，仅是数值取值，不能替代到达后持续跟踪的证明。\n\n');
fprintf(fid,'更直接地说：在残差为零的根上，当前程序返回零速度，而时变参考根通常仍在移动。因此，这个普通 ODE 的零点分支本身不能保证到达后一直精确跟踪。若用 Filippov 广义解把不连续点解释成一组允许的速度，就必须另行说明这种定义、到达条件和到达后的跟踪；这次 ode15s 输出没有完成这个证明。\n\n');

fprintf(fid,'## 最终判断与可复查材料\n\n');
fprintf(fid,'新模型已经具备“沿 KKT/PFB 残差梯度运动、没有显式系数导数、没有积分状态”的明确结构。');
if any(s.complete & s.model == "GNN")
    fprintf(fid,'部分新模型设置完成了测试，其精度应按上表具体参数评价。');
end
fprintf(fid,'但是这不自动等于已经保留 VFCR 的完整时变约束求解能力，也不自动等于更快、精确有限时间或更抗噪。');
if all(gnn.complete)
    fprintf(fid,'主设置的完成情况提供了支持，准确度与鲁棒性结论仍限定在这些题目和噪声样本。');
else
    fprintf(fid,'由于主设置存在未完成轨迹，本次不能得出新模型整体优于原 VFCR-ZNN 的结论。');
end
fprintf(fid,'当前版本不是通用的固定时间求解器，进一步实现固定时间需要改变动力学，而不是仅更换初值或调大 gamma。\n\n');
fprintf(fid,'Gamma 和高斯分布本身没有统一有限幅值上界，本次只检验 10 秒内保存的有限样本，不能推出任意噪声强度下的保证。\n\n');
fprintf(fid,'全部数据见 data/all.mat 与 data/metrics.csv；每组原始结果见 data/run_*.mat；图像见 figures；所有 %d 张 FIG 和 PNG 已重新打开或读取检查。\n\n',numel(audit.figures));
fprintf(fid,'MATLAB：%s；原文件 SHA-256 已记录于 data/meta.mat，运行前后均验证一致。\n',meta.release);
end

function value = median0(x)
if isempty(x), value = NaN; else, value = median(x); end
end

function value = check0(rows)
if isempty(rows)
    value = 'NA';
elseif ~any(rows.tested)
    value = '不适用';
elseif any(rows.motion == 0)
    value = sprintf('%d 组异常',sum(rows.motion == 0));
else
    value = '未发现此类异常';
end
end

function value = number0(x)
if isnan(x), value = 'NA'; elseif isinf(x), value = 'Inf'; else, value = sprintf('%.3e',x); end
end

function name = label0(kind)
switch kind
    case "zero", name = '无噪声';
    case "constant", name = '常数';
    case "linear", name = '线性';
    case "cosine", name = '余弦';
    case "gamma", name = 'Gamma';
    case "gaussian", name = '高斯';
    otherwise, name = '谐波';
end
end
