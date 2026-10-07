function pp_vplots(name)
% 用统一加密网格比较新模型和新跑的 VFCR。
if nargin < 1, name = 'figures'; end
folder = fileparts(mfilename('fullpath'));
vfile = fullfile(folder,'vfcr','all.mat');
if isfile(vfile)
    v = load(vfile,'results');
else
    v.results = cell(2,1);
    for k = 1:2
        file = fullfile(folder,'vfcr',sprintf('run_%02d.mat',k));
        assert(isfile(file),'VFCR-S 新精度组尚未保存。');
        saved = load(file,'r'); v.results{k} = saved.r;
    end
end
main = load(fullfile(folder,'data','all.mat'),'results','ref','c');
extra = load(fullfile(folder,'extra','all.mat'),'results');
c = main.c; ref = main.ref; t = ref.t(:);
results = cell(4,1);
for k = 1:numel(v.results)
    r = v.results{k};
    if isempty(r) || ~r.complete, continue, end
    model = modelof(r);
    if model ~= "VFCR-S", continue, end
    [rtol,atol] = tolerance(r);
    if rtol == 1e-7 && atol == 1e-9
        assert(isempty(results{1}),'VFCR-S 主精度组重复。'); results{1} = r;
    elseif rtol == 1e-9 && atol == 1e-11
        assert(isempty(results{2}),'VFCR-S 严格精度组重复。'); results{2} = r;
    end
end
assert(~isempty(results{1}) && ~isempty(results{2}),'两个 VFCR-S 精度组必须完整返回。');
results{3} = main.results{6}; results{4} = extra.results{1};
assert(modelof(results{3}) == "M5" && modelof(results{4}) == "M6",'新模型结果顺序不一致。');
assert(results{4}.c.kappa == 0.01,'M6 对照必须使用 kappa=0.01。');
fields = {'residual','xd','x','g'};
for k = 1:4
    r = results{k};
    assert(r.complete && string(r.status) == "complete",'只绘制完整新结果。');
    assert(isequal(r.t(:),t),'共同时间网格不一致。');
    assert(size(r.g,2) == 7 && size(r.g,1) == numel(t),'g 必须是共同七维状态。');
    assert(norm(r.g(1,:).'-c.g0) <= 1e-13,'共同初值不一致。');
    assert(r.c.delta == c.delta,'共同 PFB 扰动参数不一致。');
    [rtol,atol] = tolerance(r);
    if k > 1, assert(rtol == 1e-9 && atol == 1e-11,'公平主比较必须使用严格精度。'); end
    for field = fields
        value = r.metrics.(field{1});
        assert(numel(value) == numel(t) && all(isfinite(value)) && all(value >= 0), ...
            '误差序列必须为共同网格上的有限非负数。');
    end
end
out = fullfile(folder,'vfcr',name);
files = {'comparison','switch'};
for k = 1:numel(files)
    assert(~isfile(fullfile(out,[files{k},'.png'])) && ...
        ~isfile(fullfile(out,[files{k},'.fig'])),'图已存在，请传入新的图文件夹名称。');
end
if ~isfolder(out), mkdir(out); end
labels = ["VFCR-S main","VFCR-S strict","M5 strict","M6 kappa=0.01 strict"];
colors = [0.64,0.51,0.83;0.43,0.27,0.69;0.87,0.43,0.08;0.39,0.55,0.17];
styles = {':','-','--','-.'};
titles = {'KKT/PFB residual','Decision error to perturbed PFB root', ...
    'Decision error to original QP solution','Primal/dual state error to perturbed PFB root'};
ylabels = {'||xi(g,t)||','||x - x_delta*||','||x - x*||','||g - g_delta*||'};
fig = canvas(); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j = 1:4
    ax = nexttile(tiles);
    draw(ax,t,valuesof(results,fields{j}),labels,colors,styles);
    decorate(ax,titles{j},ylabels{j});
    xlim(ax,[c.start,c.finish]); xlabel(ax,'Time t (s)','Interpreter','none');
    if j == 1
        key = legend(ax,'Orientation','horizontal','NumColumns',4,'FontSize',10,'Interpreter','none');
        key.Layout.Tile = 'north';
    end
end
title(tiles,{'Example 2, zero noise | same initial primal/dual state, roots and dense output grid', ...
    'main: RelTol/AbsTol = 1e-7/1e-9 | strict: 1e-9/1e-11'}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'comparison');

center = c.peak(2); take = abs(t-center) <= 0.015+1e-12;
assert(nnz(take) > 100,'切换附近必须保留加密输出。');
tau = 1000*(t(take)-center);
fig = canvas(); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j = 1:4
    ax = nexttile(tiles); values = valuesof(results,fields{j});
    draw(ax,tau,values(take,:),labels,colors,styles);
    decorate(ax,titles{j},ylabels{j}); xlim(ax,[-15,15]); xticks(ax,-15:5:15);
    xline(ax,0,':','Color',[0.55,0.55,0.55],'HandleVisibility','off');
    xlabel(ax,'Time from peak (ms)','Interpreter','none');
    if j == 1
        key = legend(ax,'Orientation','horizontal','NumColumns',4,'FontSize',10,'Interpreter','none');
        key.Layout.Tile = 'north';
    end
end
title(tiles,{sprintf('Dense switching comparison | reference velocity peak at %.12f s',center), ...
    'main: RelTol/AbsTol = 1e-7/1e-9 | strict: 1e-9/1e-11'}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'switch');
fprintf('VFCR 比较图已保存至 %s。\n',out);
end

function name = modelof(r)
% 兼容两种实验记录的模型名称。
if isfield(r,'model'), name = string(r.model);
else, name = string(r.job.model); end
end

function [rtol,atol] = tolerance(r)
% 从保存的实际实验设置读取精度。
if isfield(r,'c') && isfield(r.c,'rtol') && isfield(r.c,'atol')
    rtol = r.c.rtol; atol = r.c.atol;
else
    rtol = r.job.rtol; atol = r.job.atol;
end
end

function values = valuesof(results,name)
% 读取四个新结果的共同误差序列。
values = zeros(numel(results{1}.t),4);
for k = 1:4, values(:,k) = results{k}.metrics.(name); end
end

function fig = canvas()
% 使用隐藏图窗和统一字体。
fig = figure('Visible','off','Color','white','Position',[50,50,1500,950], ...
    'DefaultAxesFontName','Arial','DefaultAxesFontSize',11, ...
    'DefaultTextInterpreter','none','DefaultAxesTickLabelInterpreter','none', ...
    'DefaultLegendInterpreter','none');
end

function draw(ax,t,values,labels,colors,styles)
% 为同一模型保持相同颜色和线型。
hold(ax,'on');
for k = 1:4
    plot(ax,t,values(:,k),'Color',colors(k,:),'LineStyle',styles{k}, ...
        'LineWidth',1.4,'DisplayName',labels(k));
end
end

function decorate(ax,name,label)
% 使用对数误差轴和浅色网格。
title(ax,name,'Interpreter','none','FontSize',11,'FontWeight','bold');
ylabel(ax,label,'Interpreter','none');
set(ax,'YScale','log','Box','on','GridColor',[0.65,0.65,0.65], ...
    'GridAlpha',0.18,'MinorGridAlpha',0.08,'LineWidth',0.8,'TickDir','out');
grid(ax,'on');
end

function saveplot(fig,folder,name)
% 保存后检查图片尺寸并重开可编辑图。
png = fullfile(folder,[name,'.png']); file = fullfile(folder,[name,'.fig']);
exportgraphics(fig,png,'Resolution',150); savefig(fig,file); close(fig);
info = imfinfo(png); assert(info.Width >= 1000 && info.Height >= 700);
check = openfig(file,'invisible'); assert(numel(findall(check,'Type','axes')) >= 4);
close(check);
end
