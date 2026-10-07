function pp_figures(source,name)
% 从完整新结果生成三张科学比较图。
if nargin < 1, source = 'data'; end
if nargin < 2, name = 'figures'; end
folder = fileparts(mfilename('fullpath'));
out = fullfile(folder,name); file = fullfile(folder,source,'all.mat');
assert(isfile(file),'主实验结果尚未保存。');
saved = load(file,'results','ref','c');
results = saved.results(1:4); ref = saved.ref; c = saved.c;
models = ["M3","M5","M6","M7"];
t = ref.t(:);
for k = 1:4
    r = results{k};
    assert(r.complete && r.status == "complete",'只绘制实际完成的主实验。');
    assert(r.model == models(k) && isequal(r.t,t),'模型顺序或时间网格不一致。');
end
files = {'main','switch','prediction'};
for k = 1:numel(files)
    assert(~isfile(fullfile(out,[files{k},'.png'])) && ...
        ~isfile(fullfile(out,[files{k},'.fig'])),'图已存在，请传入新的图文件夹名称。');
end
if ~isfolder(out), mkdir(out); end
colors = [0.04,0.35,0.72;0.87,0.43,0.08;0.39,0.55,0.17;0.68,0.23,0.47];
styles = {'-','--','-.',':'}; grey = [0.18,0.18,0.18];
fields = {'residual','xd','x','g'};
titles = {'KKT/PFB residual','Decision error to perturbed PFB root', ...
    'Decision error to original QP solution','Full-state error to perturbed PFB root'};
labels = {'||xi(g,t)||','||x - x_delta*||','||x - x*||','||g - g_delta*||'};
fig = canvas(1400,900); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j = 1:4
    ax = nexttile(tiles); values = valuesof(results,fields{j});
    draw(ax,t,values,models,colors,styles);
    style(ax,titles{j},labels{j},true); xlim(ax,[c.start,c.finish]);
    xlabel(ax,'Time t (s)','Interpreter','none');
    if j == 1
        legend(ax,'Location','southwest','NumColumns',2,'Interpreter','none','FontSize',10);
    end
end
title(tiles,{sprintf('Example 2, zero noise | gamma = %g, delta = %.0e',c.gamma,c.delta), ...
    sprintf('Saved main trajectories: RelTol = %.0e, AbsTol = %.0e',c.rtol,c.atol)}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'main',4);

center = c.peak(2); take = abs(t-center) <= 0.015+1e-12;
assert(nnz(take) > 100,'切换附近的采样点不足。');
tau = 1000*(t(take)-center);
fields = {'residual','g','xd','predictor'};
titles = {'KKT/PFB residual','Full-state error to perturbed PFB root', ...
    'Decision error to perturbed PFB root','Predictor error at the evolving model state'};
labels = {'||xi(g,t)||','||g - g_delta*||','||x - x_delta*||','||dp(g,t) - root velocity||'};
fig = canvas(1400,900); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j = 1:4
    ax = nexttile(tiles); values = valuesof(results,fields{j});
    draw(ax,tau,values(take,:),models,colors,styles);
    xline(ax,0,':','Color',[0.55,0.55,0.55],'HandleVisibility','off');
    style(ax,titles{j},labels{j},true); xlim(ax,[-15,15]); xticks(ax,-15:5:15);
    xlabel(ax,'Time from peak (ms)','Interpreter','none');
    if j == 1
        legend(ax,'Location','southwest','NumColumns',2,'Interpreter','none','FontSize',10);
    end
end
title(tiles,{sprintf('Dense switching window | peak time = %.12f s',center), ...
    'Independent PFB root and barrier-KKT velocity; M3 uses dp = 0'}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'switch',4);

fig = canvas(1500,1200); tiles = tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
for j = 1:2
    center = c.peak(j); take = abs(t-center) <= 0.015+1e-12;
    tau = 1000*(t(take)-center);
    ax = nexttile(tiles);
    values = ref.speed(take);
    for k = 2:4, values(:,k) = vecnorm(results{k}.metrics.dp(take,:),2,2); end
    names = ["Root velocity","M5 dp(g,t)","M6 dp(g,t)","M7 dp(g,t)"];
    draw(ax,tau,values,names,[grey;colors(2:4,:)],{'-',styles{2:4}});
    style(ax,sprintf('Peak %d: root speed and actual predictor speed',j),'Velocity norm (state units/s)',false);
    xlim(ax,[-15,15]); xticks(ax,-15:5:15); ylim(ax,[0,5500]);
    xlabel(ax,sprintf('Time from %.9f s (ms)',center),'Interpreter','none');
    xline(ax,0,':','Color',[0.55,0.55,0.55],'HandleVisibility','off');
    legend(ax,'Location','northeast','NumColumns',2,'Interpreter','none','FontSize',9);
end
center = c.peak(2); take = abs(t-center) <= 0.015+1e-12;
tau = 1000*(t(take)-center);
fields = {'oracle','predictor','sigma','lambda'};
titles = {'Predictor error evaluated at the reference root', ...
    'Predictor error evaluated at the evolving model state', ...
    'Smallest singular value of physical-state Jacobian J', ...
    'Damping used along the evolving trajectories'};
labels = {'||dp(g_delta*,t) - root velocity||','||dp(g,t) - root velocity||', ...
    'sigma_min(J)','lambda'};
for j = 1:4
    ax = nexttile(tiles); values = valuesof(results,fields{j});
    if j == 3
        draw(ax,tau,[ref.sigma(take),values(take,:)], ...
            ["Root J",models],[grey;colors],{'-',styles{:}});
    else
        draw(ax,tau,values(take,:),models,colors,styles);
    end
    style(ax,titles{j},labels{j},true); xlim(ax,[-15,15]); xticks(ax,-15:5:15);
    xlabel(ax,'Time from second peak (ms)','Interpreter','none');
    xline(ax,0,':','Color',[0.55,0.55,0.55],'HandleVisibility','off');
    if j == 3
        legend(ax,'Location','southwest','NumColumns',3,'Interpreter','none','FontSize',9);
    else
        legend(ax,'Location','southwest','NumColumns',2,'Interpreter','none','FontSize',9);
    end
end
title(tiles,{'Both peak speeds; remaining diagnostics near the second peak', ...
    'Predictors remain damped; M3 and M5 fixed-lambda curves coincide'}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'prediction',6);
fprintf('三张 PNG 与 FIG 已保存至 %s。\n',out);
end

function fig = canvas(width,height)
% 使用统一字体和隐藏图窗。
fig = figure('Visible','off','Color','white','Position',[50,50,width,height], ...
    'DefaultAxesFontName','Arial','DefaultAxesFontSize',11, ...
    'DefaultTextInterpreter','none','DefaultAxesTickLabelInterpreter','none', ...
    'DefaultLegendInterpreter','none');
end

function values = valuesof(results,name)
% 读取同一指标的四个真实模型序列。
values = zeros(numel(results{1}.t),4);
for k = 1:4, values(:,k) = results{k}.metrics.(name); end
assert(all(isfinite(values),'all') && all(values >= 0,'all'),'图中指标必须为有限非负数。');
end

function draw(ax,t,values,names,colors,styles)
% 固定每个模型的颜色和线型。
hold(ax,'on');
for k = 1:size(values,2)
    plot(ax,t,values(:,k),'Color',colors(k,:),'LineStyle',styles{k}, ...
        'LineWidth',1.4,'DisplayName',names(k));
end
end

function style(ax,name,label,log)
% 保留清楚的轴标签和浅网格。
title(ax,name,'Interpreter','none','FontSize',11,'FontWeight','bold');
ylabel(ax,label,'Interpreter','none');
set(ax,'Box','on','GridColor',[0.65,0.65,0.65],'GridAlpha',0.18, ...
    'MinorGridAlpha',0.08,'LineWidth',0.8,'TickDir','out');
grid(ax,'on');
if log, set(ax,'YScale','log'); end
end

function saveplot(fig,folder,name,count)
% 保存后重开图文件并核对输出尺寸。
png = fullfile(folder,[name,'.png']); file = fullfile(folder,[name,'.fig']);
exportgraphics(fig,png,'Resolution',150);
savefig(fig,file); close(fig);
info = imfinfo(png); assert(info.Width >= 1000 && info.Height >= 700);
check = openfig(file,'invisible');
assert(numel(findall(check,'Type','axes')) >= count);
close(check);
end
