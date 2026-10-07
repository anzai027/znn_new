function conclude(results,refs,c,meta,folder)
% 保存图像并写出平滑试验结论。
stats = readtable(fullfile(folder,'data','metrics.csv'),'TextType','string');
old = readtable(fullfile(meta.prior,'data','metrics.csv'),'TextType','string');
old = old(ismember(old.id,stats.base),:);
writetable(old,fullfile(folder,'data','baseline.csv'),'Encoding','UTF-8');
for model = ["GNN","VFCR-ZNN"]
    for kind = ["zero","cosine"]
        fig = figure('Visible','off','Color','w','Position',[50,50,1200,470]);
        ax = [subplot(1,2,1),subplot(1,2,2)];
        hold(ax(1),'on'); hold(ax(2),'on');
        ids = find(stats.group == "scan" & stats.model == model & stats.noise == kind);
        for id = ids.'
            r = results{stats.id(id)};
            label = sprintf('eps=%g',r.job.eps);
            if ~r.complete
                label = [label,' (interrupted)'];
            end
            info = struct('last',r.last,'complete',r.complete);
            semilogy(ax(1),r.t,positive(r.metrics.e),'LineWidth',1.1,'DisplayName',label,'UserData',info);
            semilogy(ax(2),r.t,positive(r.metrics.ed),'LineWidth',1.1,'DisplayName',label,'UserData',info);
        end
        yline(ax(1),1e-4,'k:','DisplayName','threshold 1e-4');
        yline(ax(1),1e-6,'k--','DisplayName','threshold 1e-6');
        for a = ax
            set(a,'YScale','log');
            xlabel(a,'t (s)'); xlim(a,[0,10]); grid(a,'on');
            legend(a,'Location','best');
        end
        ylabel(ax(1),'KKT/PFB residual'); ylabel(ax(2),'error to delta reference');
        sgtitle(sprintf('%s | %s | smooth scale scan',model,kind),'Interpreter','none');
        savepic(fig,lower(model)+"_"+kind+"_scan");
    end
end
for eps = [1e-3,1e-2]
for kind = ["zero","cosine"]
    fig = figure('Visible','off','Color','w','Position',[50,50,1200,470]);
    for j = 1:2
        ax = subplot(1,2,j); hold(ax,'on');
        plot(ax,refs{2}.t,refs{2}.xd(:,j),'k--','LineWidth',1.5,'DisplayName','delta reference');
        for model = ["GNN","VFCR-ZNN"]
            id = find(stats.group == "scan" & stats.model == model ...
                & stats.noise == kind & stats.eps == eps,1);
            r = results{stats.id(id)};
            label = model;
            if ~r.complete
                label = label+" (interrupted)";
            end
            info = struct('last',r.last,'complete',r.complete);
            plot(ax,r.t,r.y(:,j),'LineWidth',1.1,'DisplayName',label,'UserData',info);
        end
        xlabel(ax,'t (s)'); ylabel(ax,sprintf('x_%d',j));
        grid(ax,'on'); legend(ax,'Location','best'); xlim(ax,[0,10]);
    end
    sgtitle(sprintf('Example 2 | %s | eps=%g',kind,eps));
    savepic(fig,"states_"+kind+"_"+string(eps));
end
end
for eps = [1e-3,1e-2]
fig = figure('Visible','off','Color','w','Position',[50,50,1200,470]);
for j = 1:2
    model = ["GNN","VFCR-ZNN"]; model = model(j);
    ax = subplot(1,2,j); hold(ax,'on');
    ids = find(stats.group == "random" & stats.model == model & stats.eps == eps);
    for id = ids.'
        r = results{stats.id(id)};
        semilogy(ax,r.t,positive(r.metrics.ed),'LineWidth',0.9, ...
            'DisplayName',sprintf('seed %d',r.job.seed));
    end
    set(ax,'YScale','log'); xlabel(ax,'t (s)'); ylabel(ax,'error to delta reference');
    title(ax,model); grid(ax,'on'); legend(ax,'Location','best'); xlim(ax,[0,10]);
end
sgtitle(sprintf('Gaussian state noise | eps=%g',eps));
savepic(fig,"gaussian_"+string(eps));
end
audit = auditfig(folder,results,c);
save(fullfile(folder,'data','audit.mat'),'audit');
save(fullfile(folder,'data','all.mat'),'audit','-append');
file = fullfile(folder,'结论.md');
fid = fopen(file,'w','n','UTF-8');
guard = onCleanup(@() fclose(fid));
fprintf(fid,'# 平滑能否缓解求解失败\n\n');
fprintf(fid,'共尝试 %d 组，完整到 10 秒 %d 组，中断 %d 组；GNN 速度检查异常 %d 组。\n\n', ...
    height(stats),sum(stats.complete),sum(~stats.complete),sum(stats.tested & stats.motion == 0));
fprintf(fid,'原公式结果来自上一次保存的数据；本次不修改原 system，也不覆盖旧结果。\n\n');
fprintf(fid,'**结论：平滑确实缓解了数值失败，但尚未让新 GNN 达到所要求的高精度 TVQP 跟踪。** 在 eps=1e-2 的对应 13 组中，GNN 从原式的 7 组中断、6 组完整但速度异常，变成 12 组完整且通过速度检查、1 组达到调用上限；VFCR 从 8/13 组完整变成 13/13。\n\n');
fprintf(fid,'GNN eps=1e-2 无噪声的 delta 解 RMSE 约 0.0304、5～10 秒最大残差约 0.185，没有持续达到 1e-4；严格容差所得 RMSE 相同，因此该误差并非单纯放宽容差导致。仍然中断的是 Example 1；该题将无穷边界替换为 1e8，具体失败原因尚需进一步定位。\n\n');
fprintf(fid,'对 VFCR，平滑后无噪声能完整运行并保持较高精度；对新 GNN，平滑更适合作为继续研究的数值诊断，当前不能宣称已满足准确求解、固定时间收敛或全面抗噪的目标。\n\n');
fprintf(fid,'## 平滑方式与设置\n\n');
fprintf(fid,'GNN：令 s=J^T xi，将分母 norm(s)^p 换成 (norm(s)^2+eps^2)^(p/2)，保持 gamma=100、p=1。\n\n');
fprintf(fid,'VFCR：eps 外严格保持原激活函数；eps 内用奇三次函数连接，使函数和一阶导数在边界连续，通过 problem.phi 调用原 my_system。\n\n');
fprintf(fid,'四尺度为 1e-4、1e-3、1e-2、1e-1；主尺度 1e-3 是运行前指定的诊断设置，扫描发现 GNN 无噪声未完成后，追加 1e-2 的同样补充检查，这是依据扫描选择的探索设置。两模型尺度作用于不同量，不能理解为同一种参数。\n\n');
fprintf(fid,'采用 Example 2、共同状态噪声、0～10 秒、ode15s、RelTol=1e-5、AbsTol=1e-7、MaxStep=0.02；严格检查为 1e-7/1e-9。\n\n');
fprintf(fid,'高斯噪声沿用原种子 1～5 的样本，每 0.01 秒保持并分段积分；初值、delta=1e-4 和所有题目系数沿用原实验。\n\n');
fprintf(fid,'每次限 180 秒或 50 万次调用；中断之后标为缺失，完整运行不等同于准确求解。\n\n');
fprintf(fid,'## 四尺度扫描\n\n');
fprintf(fid,'误差、违反量和残差使用 5～10 秒，Tc 要求最后持续达标且至少观察 0.5 秒。\n\n');
fprintf(fid,'| 模型 | 噪声 | eps | 完整 | 秒 | 调用数 | 真解 RMSE | delta 解 RMSE | 最大残差 | 等式违反 | 不等式违反 | Tc 1e-6 | Tc 1e-4 |\n');
fprintf(fid,'|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
ids = find(stats.group == "scan");
for id = ids.'
    s = stats(id,:);
    fprintf(fid,'| %s | %s | %.0e | %d | %.4g | %d | %.4g | %.4g | %.4g | %.4g | %.4g | %.4g | %.4g |\n', ...
        s.model,s.noise,s.eps,s.complete,s.seconds,s.calls,s.rmse,s.rmsed, ...
        s.residual,s.equality,s.inequality,s.tc6,s.tc4);
end
fprintf(fid,'\n## 与原式失败情况比较\n\n');
fprintf(fid,'两个尺度各比较相同题目、初值、噪声、种子、容差和重复次数，每个尺度共 26 组；原式对应同一批基线，不能当成两批独立样本。\n\n');
fprintf(fid,'| 模型 | eps | 组数 | 原式中断 | 原式完整但速度异常 | 平滑中断 | 平滑完整但速度异常 |\n');
fprintf(fid,'|---|---:|---:|---:|---:|---:|---:|\n');
for eps = [1e-3,1e-2]
for model = ["GNN","VFCR-ZNN"]
    s = stats(stats.model == model & stats.eps == eps,:);
    b = old(ismember(old.id,s.base),:);
    fprintf(fid,'| %s | %.0e | %d | %d | %d | %d | %d |\n',model,eps,height(s), ...
        sum(~b.complete),sum(b.complete & b.tested & b.motion == 0), ...
        sum(~s.complete),sum(s.complete & s.tested & s.motion == 0));
end
end
fprintf(fid,'\n## 三次计时与精度\n\n');
fprintf(fid,'耗时是完整积分到 10 秒的耗时，不是数学收敛时间；原式和本次来自不同批次，不是同时运行的硬件基准。\n\n');
fprintf(fid,'| 模型 | eps | 噪声 | 原式完整且通过检查 | 原式耗时中位数 | 平滑完整且通过检查 | 平滑耗时中位数 | 平滑 delta RMSE 中位数 |\n');
fprintf(fid,'|---|---:|---|---:|---:|---:|---:|---:|\n');
for eps = [1e-3,1e-2]
for model = ["GNN","VFCR-ZNN"]
    for kind = ["zero","cosine"]
        s = stats(stats.model == model & stats.noise == kind & stats.eps == eps ...
            & ismember(stats.group,["scan","timing"]),:);
        b = old(ismember(old.id,s.base),:);
        good = s.complete & (~s.tested | s.motion == 1);
        good0 = b.complete & (~b.tested | b.motion == 1);
        fprintf(fid,'| %s | %.0e | %s | %d/3 | %.4g | %d/3 | %.4g | %.4g |\n', ...
            model,eps,kind,sum(good0),median(b.seconds(good0)),sum(good), ...
            median(s.seconds(good)),median(s.rmsed(good)));
    end
end
end
fprintf(fid,'\n## 严格容差检查\n\n');
for eps = [1e-3,1e-2]
for model = ["GNN","VFCR-ZNN"]
    i = find(stats.model == model & stats.group == "scan" & stats.noise == "zero" & stats.eps == eps,1);
    j = find(stats.model == model & stats.group == "precision" & stats.eps == eps,1);
    a = results{stats.id(i)}; b = results{stats.id(j)};
    error = NaN;
    if a.complete && b.complete
        error = max(vecnorm(a.y(:,1:2)-b.y(:,1:2),2,2));
    end
    fprintf(fid,'%s eps=%g：普通/严格完整情况 %d/%d，delta RMSE %.5g/%.5g，两个数值轨迹最大距离 %.5g。\n\n', ...
        model,eps,a.complete,b.complete,a.metrics.rmsed,b.metrics.rmsed,error);
end
end
fprintf(fid,'## 高斯噪声检查\n\n');
for eps = [1e-3,1e-2]
for model = ["GNN","VFCR-ZNN"]
    s = stats(stats.model == model & stats.group == "random" & stats.eps == eps,:);
    good = s.complete & (~s.tested | s.motion == 1);
    fprintf(fid,'%s eps=%g：完整且通过检查 %d/5；这些组的 delta RMSE 中位数 %.5g、最大残差中位数 %.5g。\n\n', ...
        model,eps,sum(good),median(s.rmsed(good)),median(s.residual(good)));
end
end
fprintf(fid,'## 可以得出的结论与边界\n\n');
fprintf(fid,'这次用于判断数值中断能否缓解；是否准确必须同时看残差、误差和约束，不能只看完成次数。速度上限检查只是必要检查，通过也不构成充分证明。\n\n');
fprintf(fid,'GNN 平滑后仍只有七个状态，没有积分状态，不读取系数导数，不执行 J\\右端项；原文件校验全部未变。\n\n');
fprintf(fid,'平滑改变了公式，不能当成原式严格复现，也不能直接沿用原有限时间理论。静态一维 p=1 例子变成 dx/dt=-gamma*x/sqrt(x^2+eps^2)，零点附近速度约为 -gamma*x/eps，因此通常是逐渐趋近零，不再具有原式的严格有限时间到达性质。\n\n');
fprintf(fid,'本次没有检查其他单幂指数、全部噪声类型和大初值，也不能证明对任意时变题目都收敛或固定时间收敛。\n\n');
fprintf(fid,'全部 %d 对 FIG/PNG 已重新打开读取；每次数据用 save 保存，汇总在 data/all.mat，逐组指标在 data/metrics.csv。\n',numel(audit.names));

    function savepic(fig,name)
        name = replace(name,"-","_");
        savefig(fig,fullfile(folder,'figures',name+".fig"));
        exportgraphics(fig,fullfile(folder,'figures',name+".png"),'Resolution',140);
        close(fig);
    end
end

function audit = auditfig(folder,results,c)
files = dir(fullfile(folder,'figures','*.fig'));
audit.names = strings(numel(files),1);
audit.lines = zeros(numel(files),1);
for k = 1:numel(files)
    fig = openfig(fullfile(files(k).folder,files(k).name),'invisible');
    axes = findall(fig,'Type','axes');
    for ax = axes.'
        assert(~isempty(ax.XLabel.String) && ~isempty(ax.YLabel.String));
        assert(~isempty(findall(ax,'Type','line')));
    end
    assert(~isempty(findall(fig,'Type','legend')));
    lines = findall(fig,'Type','line');
    for line = lines.'
        info = line.UserData;
        if isstruct(info) && isfield(info,'complete') && ~info.complete
            assert(all(isnan(line.YData(line.XData > info.last+1e-10))));
        end
    end
    audit.lines(k) = numel(findall(fig,'Type','line'));
    close(fig);
    [~,name] = fileparts(files(k).name);
    pic = imread(fullfile(folder,'figures',[name,'.png']));
    assert(~isempty(pic) && max(pic(:)) > min(pic(:)));
    audit.names(k) = string(name);
end
for k = 1:numel(results)
    r = results{k};
    assert(isequaln(r.metrics.motion,motion(r,c)));
    assert(all(isfinite(r.raw.y),'all') && all(diff(r.raw.t) > 0));
    if ~r.complete
        assert(all(isnan(r.y(c.t > r.last+1e-10,:)),'all'));
    end
end
audit.passed = true;
end

function y = positive(y)
% 保留缺失值并限制对数坐标的下界。
y(isfinite(y) & y < 1e-16) = 1e-16;
end
