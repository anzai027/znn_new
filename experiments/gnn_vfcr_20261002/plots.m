function plots(results,refs,c,folder)
% 绘制成功和中断的轨迹，缺失部分不补造。
set(groot,'DefaultAxesFontName','Microsoft YaHei');
set(groot,'DefaultTextFontName','Microsoft YaHei');
stats = collect(results);
colors = [0.10,0.38,0.72;0.86,0.22,0.16];
models = ["VFCR-ZNN","GNN"];
test = ["zero","cosine"];
f = fresh('不同初值：状态轨迹',[100,50,1100,1250]);
layout = tiledlayout(4,2,'TileSpacing','compact');
for k = 1:4
    for j = 1:2
        ax = nexttile;
        plot(c.t,refs{1}.x(:,j),'k--','LineWidth',1.2,'DisplayName','真正最优解');
        hold on
        for z = 1:2
            id = find(stats.group == "initial" & stats.scale == c.scales(k) & stats.model == models(z),1);
            line0(ax,results{stats.id(id)},'x',j,colors(z,:),false);
        end
        title(sprintf('初值范围 [0,%g]，x_%d',c.scales(k),j));
        xlabel('时间 (s)'); ylabel(sprintf('x_%d',j));
        finish(ax);
    end
end
title(layout,'论文 Example 1：实线只覆盖实际计算到的时间');
save0(f,'initial_states');

f = fresh('不同初值：残差',[120,100,1100,760]);
layout = tiledlayout(2,2,'TileSpacing','compact');
for k = 1:4
    ax = nexttile;
    hold on
    for z = 1:2
        id = find(stats.group == "initial" & stats.scale == c.scales(k) & stats.model == models(z),1);
        line0(ax,results{stats.id(id)},'e',1,colors(z,:),true);
    end
    yline(c.threshold(1),'k:','1e-6','HandleVisibility','off');
    title(sprintf('初值范围 [0,%g]',c.scales(k)));
    xlabel('时间 (s)'); ylabel('残差 ||xi||_2');
    finish(ax);
end
title(layout,'Example 1：横线是共同阈值，叉号是中断位置');
save0(f,'initial_residuals');

f = fresh('主实验状态',[120,100,1100,760]);
layout = tiledlayout(2,2,'TileSpacing','compact');
kinds = ["zero","cosine"];
for k = 1:2
    for j = 1:2
        ax = nexttile;
        plot(c.t,refs{2}.x(:,j),'k--','DisplayName','真正最优解','LineWidth',1.2);
        hold on
        plot(c.t,refs{2}.xd(:,j),':','Color',[0.35,0.35,0.35], ...
            'DisplayName','delta 近似解','LineWidth',1.1);
        for z = 1:2
            r = main0(models(z),kinds(k),1);
            line0(ax,r,'x',j,colors(z,:),false);
        end
        title(sprintf('%s，x_%d',label0(kinds(k)),j));
        xlabel('时间 (s)'); ylabel(sprintf('x_%d',j));
        finish(ax);
    end
end
title(layout,'Example 2：gamma=100、p=1；速度异常的 GNN 输出仅用于诊断');
save0(f,'main_states');

kinds = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"];
fields = ["e","ex","ed","grad"];
names = ["残差 ||xi||_2","与真正解的距离","与 delta 解的距离","梯度 ||J^T xi||_2"];
files = ["main_residuals","main_error","main_delta","main_gradient"];
for j = 1:4
    f = fresh(names(j),[80,40,1250,1050]);
    layout = tiledlayout(3,3,'TileSpacing','compact');
    for k = 1:7
        ax = nexttile;
        hold on
        for z = 1:2
            line0(ax,main0(models(z),kinds(k),1),fields(j),1,colors(z,:),true);
        end
        if j == 1
            yline(c.threshold(1),'k:','HandleVisibility','off');
        end
        title(label0(kinds(k)));
        xlabel('时间 (s)'); ylabel(names(j));
        finish(ax);
    end
    ax = nexttile;
    axis(ax,'off');
    text(ax,0,0.85,{'蓝线：原 VFCR-ZNN';'红线：新 GNN，gamma=100、p=1'; ...
        '两者共用状态噪声和初值';'随机噪声这里显示种子 1'; ...
        '叉号：求解器中断';'未完成的区间不画曲线';'零值只在对数图上显示为 1e-16'}, ...
        'FontSize',12,'VerticalAlignment','top','Interpreter','none');
    title(layout,'Example 2：保存的求解器输出；速度异常的曲线仅用于诊断');
    save0(f,char(files(j)));
end

f = fresh('参数对比',[80,80,1200,850]);
layout = tiledlayout(2,2,'TileSpacing','compact');
colors0 = lines(6);
for k = 1:2
    for j = 1:2
        ax = nexttile;
        hold on
        field = 'e';
        name = '残差';
        if j == 2
            field = 'ed';
            name = '与 delta 解的距离';
        end
        line0(ax,main0("VFCR-ZNN",test(k),1),field,1,colors0(1,:),true);
        line0(ax,main0("GNN",test(k),1),field,1,colors0(2,:),true);
        for z = 1:size(c.pairs,1)
            id = find(stats.group == "parameter" & stats.noise == test(k) ...
                & stats.gamma == c.pairs(z,1) & stats.p == c.pairs(z,2),1);
            r = results{stats.id(id)};
            line0(ax,r,field,1,colors0(z+2,:),true);
        end
        title(sprintf('%s：%s',label0(test(k)),name));
        xlabel('时间 (s)'); ylabel(name);
        finish(ax);
    end
end
title(layout,'补充参数保持原公式，只改变 gamma 或 p');
save0(f,'parameters');

f = fresh('p=0.5 的状态轨迹',[100,80,1100,800]);
layout = tiledlayout(2,2,'TileSpacing','compact');
for k = 1:2
    id = find(stats.group == "parameter" & stats.noise == test(k) ...
        & stats.gamma == 100 & stats.p == 0.5,1);
    for j = 1:2
        ax = nexttile;
        plot(c.t,refs{2}.x(:,j),'k--','DisplayName','真正最优解','LineWidth',1.2);
        hold on
        plot(c.t,refs{2}.xd(:,j),':','Color',[0.4,0.4,0.4], ...
            'DisplayName','delta 近似解','LineWidth',1.1);
        line0(ax,main0("VFCR-ZNN",test(k),1),'x',j,colors(1,:),false);
        line0(ax,results{stats.id(id)},'x',j,colors(2,:),false);
        title(sprintf('%s，x_%d',label0(test(k)),j));
        xlabel('时间 (s)'); ylabel(sprintf('x_%d',j));
        finish(ax);
    end
end
title(layout,'补充参数 gamma=100、p=0.5：有跟踪误差，无噪声 VFCR 中断');
save0(f,'supplement_states');

f = fresh('p=0.5 的约束违反量',[100,80,1100,800]);
layout = tiledlayout(2,2,'TileSpacing','compact');
for k = 1:2
    id = find(stats.group == "parameter" & stats.noise == test(k) ...
        & stats.gamma == 100 & stats.p == 0.5,1);
    for j = 1:2
        ax = nexttile;
        hold on
        field = 'eq';
        name = '等式违反量';
        if j == 2
            field = 'bound';
            name = '不等式违反量';
        end
        line0(ax,main0("VFCR-ZNN",test(k),1),field,1,colors(1,:),true);
        line0(ax,results{stats.id(id)},field,1,colors(2,:),true);
        title(sprintf('%s：%s',label0(test(k)),name));
        xlabel('时间 (s)'); ylabel(name);
        finish(ax);
    end
end
title(layout,'约束违反量：零值仅在对数图中显示为 1e-16');
save0(f,'constraints');

f = fresh('精度检查',[80,80,1100,800]);
layout = tiledlayout(2,2,'TileSpacing','compact');
for k = 1:2
    for j = 1:2
        ax = nexttile;
        hold on
        field = 'e';
        name = '残差';
        if j == 2
            field = 'ed';
            name = '与 delta 解的距离';
        end
        line0(ax,main0(models(k),"zero",1),field,1,colors(k,:),true);
        id = find(stats.group == "precision" & stats.model == models(k),1);
        line0(ax,results{stats.id(id)},field,1,[0.15,0.60,0.25],true);
        title(sprintf('%s：%s',models(k),name));
        xlabel('时间 (s)'); ylabel(name);
        finish(ax);
    end
end
title(layout,'绿色为更严格的容差，模型公式没有变化');
save0(f,'precision');

f = fresh('新 GNN 数值诊断',[100,80,1100,780]);
layout = tiledlayout(2,2,'TileSpacing','compact');
for k = 1:2
    r = main0("GNN",test(k),1);
    ax = nexttile;
    line0(ax,r,'grad',1,colors(2,:),true);
    xlim([0,max(r.last,1e-6)]);
    title(sprintf('%s：梯度与中断',label0(test(k))));
    xlabel('时间 (s)'); ylabel('||J^T xi||_2'); grid on
    ax = nexttile;
    step = diff(r.raw.t);
    take = unique(round(linspace(1,numel(step),min(4000,numel(step)))));
    t = r.raw.t(2:end);
    semilogy(ax,t(take),step(take),'Color',colors(2,:));
    title(sprintf('实际步长，状态：%s',r.status));
    xlabel('时间 (s)'); ylabel('相邻接受点的时间差 (s)'); grid on
end
title(layout,'诊断曲线不改变原模型，也不把中断当成收敛');
save0(f,'diagnostics');

f = fresh('p=1 的速度上限检查',[100,80,1100,800]);
layout = tiledlayout(2,2,'TileSpacing','compact');
for k = 1:2
    r = main0("GNN",test(k),1);
    for j = 1:2
        ax = nexttile((j-1)*2+k);
        t = r.t;
        g = r.y(:,1:c.l);
        name = '输出网格';
        if j == 2
            take = unique(round(linspace(1,numel(r.raw.t),min(4000,numel(r.raw.t)))));
            t = r.raw.t(take);
            g = r.raw.y(take,1:c.l);
            name = '求解器接受点';
        end
        distance = vecnorm(g-r.job.g0.',2,2);
        limit = r.metrics.motion.rate*t;
        semilogy(ax,t(2:end),max(distance(2:end),1e-16),'r-', ...
            'DisplayName','保存的移动距离','LineWidth',1.2);
        hold on
        semilogy(ax,t(2:end),limit(2:end),'k--', ...
            'DisplayName','公式允许的距离上限','LineWidth',1.2);
        title(sprintf('%s：%s',label0(test(k)),name));
        xlabel('时间 (s)'); ylabel('||g(t)-g(0)||_2');
        finish(ax);
    end
end
title(layout,'p=1：严重超限说明数值轨迹未忠实执行原微分方程');
save0(f,'motion');

f = fresh('关键指标表',[100,100,1350,950]);
ax = axes(f);
axis(ax,'off');
text(ax,0,1,'主实验：完整不等于正确，带 * 的行为速度异常','FontSize',16,'FontWeight','bold');
text(ax,0,0.955,'noise       model       done    RMSE(true)   RMSE(delta)  peak(delta)   Tc(1e-6)', ...
    'FontName','Consolas','FontSize',12,'Interpreter','none');
y = 0.91;
for kind = kinds
    for model = models
        rows = stats(stats.group == "main" & stats.noise == kind & stats.model == model,:);
        good = rows(rows.complete,:);
        a = [median0(good.rmse),median0(good.rmsed),median0(good.peakd),median0(good.tc6)];
        count = sprintf('%d/%d',height(good),height(rows));
        if any(good.tested & good.motion == 0)
            count = [count,'*'];
        end
        line = sprintf('%-11s %-11s %-7s %-12s %-12s %-12s %s', ...
            kind,model,count,number0(a(1)),number0(a(2)),number0(a(3)),number0(a(4)));
        text(ax,0,y,line,'FontName','Consolas','FontSize',12,'Interpreter','none');
        y = y-0.042;
    end
end
text(ax,0,y-0.015,'NA：运行未完成；Inf：完整运行但未满足持续达标判据。', ...
    'FontSize',12,'Interpreter','none');
text(ax,0,y-0.067,'*：仅为异常诊断，不能用于性能排名；真正解与 delta 解分别记录。', ...
    'FontSize',12,'Interpreter','none');
save0(f,'scores');

f = fresh('运行时间和成功次数',[100,100,1250,620]);
ax = axes(f);
axis(ax,'off');
text(ax,0,0.98,'三次重复：排除速度异常后再统计用时','FontSize',16,'FontWeight','bold');
text(ax,0,0.88,'noise       model       returned counted   median seconds   median RHS calls', ...
    'FontName','Consolas','FontSize',12,'Interpreter','none');
y = 0.75;
for kind = ["zero","cosine"]
    for model = models
        take = any(stats.group == ["main","timing"],2) & stats.noise == kind ...
            & stats.model == model & stats.seed == 1;
        rows = stats(take,:);
        good = rows(rows.complete & ~(rows.tested & rows.motion == 0),:);
        line = sprintf('%-11s %-11s %d/%d      %d/%d       %-16s %s',kind,model, ...
            sum(rows.complete),height(rows),height(good),height(rows), ...
            number0(median0(good.seconds)),number0(median0(good.calls)));
        text(ax,0,y,line,'FontName','Consolas','FontSize',12,'Interpreter','none');
        y = y-0.13;
    end
end
text(ax,0,0.14,'全部原始耗时见 CSV；数值异常或未完成的短耗时不能说明求解更快。', ...
    'FontSize',12,'Interpreter','none');
save0(f,'runtime');

    function value = main0(model,kind,seed)
        index = find(stats.group == "main" & stats.model == model & stats.noise == kind & stats.seed == seed,1);
        value = results{stats.id(index)};
    end

    function save0(fig,name)
        savefig(fig,fullfile(folder,'figures',[name,'.fig']));
        exportgraphics(fig,fullfile(folder,'figures',[name,'.png']),'Resolution',180);
        close(fig);
    end
end

function f = fresh(name,position)
f = figure('Visible','off','Color','w','Name',name,'Position',position);
end

function line0(ax,r,field,column,color,logarithm)
if r.complete
    m = r.metrics;
else
    m = r.metrics.raw;
end
y = m.(field);
y = y(:,column);
if logarithm
    y(y == 0) = 1e-16;
    set(ax,'YScale','log');
end
name = string(r.job.model);
style = '-';
if r.job.model == "GNN"
    name = sprintf('GNN gamma=%g p=%g',r.job.gamma,r.job.p);
end
if r.job.group == "precision"
    name = string(name)+" 严格容差";
end
if r.metrics.motion.tested && r.metrics.motion.passed == 0
    name = string(name)+"，速度异常";
    style = ':';
end
if ~r.complete
    name = string(name)+sprintf('，中断 %.3gs',r.last);
end
plot(ax,m.t,y,'Color',color,'LineStyle',style,'LineWidth',1.15,'DisplayName',name);
if ~r.complete
    plot(ax,m.t(end),y(end),'x','Color',color,'LineWidth',1.5, ...
        'MarkerSize',8,'HandleVisibility','off');
end
end

function finish(ax)
xlim(ax,[0,10]);
grid(ax,'on');
box(ax,'on');
legend(ax,'Location','best','FontSize',8,'Interpreter','none');
end

function name = label0(kind)
switch kind
    case "zero", name = '无噪声';
    case "constant", name = '常数噪声';
    case "linear", name = '线性噪声';
    case "cosine", name = '余弦噪声';
    case "gamma", name = 'Gamma 噪声';
    case "gaussian", name = '高斯噪声';
    otherwise, name = '谐波噪声';
end
end

function value = median0(x)
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end

function value = number0(x)
if isnan(x)
    value = 'NA';
elseif isinf(x)
    value = 'Inf';
else
    value = sprintf('%.3e',x);
end
end
