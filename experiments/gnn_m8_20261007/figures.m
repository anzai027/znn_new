function figures(source,name)
% 从已完成结果绘制严格比较和小对照。
if nargin < 1, source = 'data'; end
if nargin < 2, name = 'figures'; end
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
if strcmp(source,'precision'), higher(folder,name); return, end
old = fullfile(root,'experiments','gnn_predictor_20261006');
file = fullfile(folder,source,'all.mat'); assert(isfile(file),'新实验尚未完整保存。');
saved = load(file,'results','ref','c'); all = saved.results; ref = saved.ref; c = saved.c;
assert(numel(all) == 9,'正式比较需要九组完整结果。');
for k = 1:numel(all), verify(all{k},ref,c,k > 2); end
assert(all{3}.label == "M8" && all{4}.label == "M8a" && ...
    all{5}.label == "M8a fixed" && all{6}.label == "M8 lp1e-9" && ...
    all{7}.label == "M8 lp1e-8" && all{8}.label == "M8 h1e-6" && ...
    all{9}.label == "M8a h1e-6",'正式任务的模型顺序不一致。');
a = load(fullfile(old,'vfcr','run_02.mat'),'r');
b = load(fullfile(old,'data','run_06.mat'),'r');
d = load(fullfile(old,'extra','run_01.mat'),'r');
results = {a.r,b.r,d.r,all{3},all{4}};
for k = 1:numel(results), verify(results{k},ref,c,true); end
out = fullfile(folder,name); files = {'main','switch','ablation'};
for k = 1:numel(files)
    assert(~isfile(fullfile(out,[files{k},'.png'])) && ...
        ~isfile(fullfile(out,[files{k},'.fig'])),'已有图像，请指定新的输出文件夹。');
end
if ~isfolder(out), mkdir(out); end
labels = ["VFCR-S","M5","M6 kappa=0.01","M8","M8a"];
colors = [0.43,0.27,0.69;0.87,0.43,0.08;0.39,0.55,0.17;0.04,0.35,0.72;0.68,0.23,0.47];
styles = {'-','--','-.',':','-'};
t = ref.t(:); fields = {'residual','xd','x','g'};
titles = {'KKT/PFB residual','Decision error to perturbed PFB root', ...
    'Decision error to original QP solution','Primal/dual state error to perturbed PFB root'};
axes = {'||xi(g,t)||','||x - x_delta*||','||x - x*||','||g - g_delta*||'};
fig = canvas(); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j = 1:4
    ax = nexttile(tiles); draw(ax,t,valuesof(results,fields{j}),labels,colors,styles,true);
    decorate(ax,titles{j},axes{j}); xlim(ax,[c.start,c.finish]);
    xlabel(ax,'Time t (s)','Interpreter','none');
    if j == 1, key(ax,5); end
end
title(tiles,{'Example 2, zero noise | same initial primal/dual state, roots and dense grid', ...
    'All curves strict: RelTol/AbsTol = 1e-9/1e-11 | GNN h = 1e-5 s', ...
    'Log-axis limits are for display; original-QP errors include the same delta bias'}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'main');

center = c.peak(2); take = abs(t-center) <= 0.015+1e-12;
assert(nnz(take) > 100,'第二峰附近必须使用加密输出。');
tau = 1000*(t(take)-center);
fig = canvas(); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j = 1:4
    ax = nexttile(tiles); values = valuesof(results,fields{j});
    draw(ax,tau,values(take,:),labels,colors,styles,true);
    decorate(ax,titles{j},axes{j}); local(ax,15);
    if j == 1, key(ax,5); end
end
title(tiles,{sprintf('Strict switching comparison | reference speed peak at %.12f s',center), ...
    'RelTol/AbsTol = 1e-9/1e-11 | decision and multiplier errors use the same PFB root', ...
    'Log-axis limits are for display; original-QP errors include the same delta bias'}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'switch');

take = abs(t-center) <= 0.001+1e-12; tau = 1000*(t(take)-center);
assert(nnz(take) > 100,'小对照必须使用峰中心加密输出。');
fig = canvas(); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
list = {all{3},all{6},all{7}};
names = ["M8 lp=1e-10","M8 lp=1e-9","M8 lp=1e-8"];
blue = [colors(4,:);0.28,0.52,0.87;0.02,0.19,0.46];
ax = nexttile(tiles); values = valuesof(list,'g');
draw(ax,tau,values(take,:),names,blue,{':','--','-.'},false);
decorate(ax,'Prediction damping: same fixed correction and h','||g - g_delta*||');
local(ax,1); legend(ax,'Location','southwest','Interpreter','none','FontSize',10);

list = {all{4},all{5}};
names = ["M8a adaptive correction","M8a fixed correction"];
pink = [colors(5,:);0.44,0.16,0.34];
ax = nexttile(tiles); values = valuesof(list,'g');
draw(ax,tau,values(take,:),names,pink,{'-','--'},false);
decorate(ax,'Direct predictor: correction damping control','||g - g_delta*||');
local(ax,1); legend(ax,'Location','southwest','Interpreter','none','FontSize',10);

list = {all{3},all{8},all{4},all{9}};
names = ["M8 h=1e-5","M8 h=1e-6","M8a h=1e-5","M8a h=1e-6"];
shades = [colors(4,:);0.28,0.52,0.87;colors(5,:);0.44,0.16,0.34];
for j = 1:2
    ax = nexttile(tiles); field = 'xd';
    if j == 2, field = 'oracle'; end
    values = valuesof(list,field); draw(ax,tau,values(take,:),names,shades,{':','--','-','-.'},false);
    if j == 1
        decorate(ax,'Time-difference step: actual decision error','||x - x_delta*||');
    else
        decorate(ax,'Time-difference step: predictor tested at the reference root','||dp(g_delta*,t) - root velocity||');
    end
    local(ax,1); legend(ax,'Location','southwest','NumColumns',2,'Interpreter','none','FontSize',9);
end
title(tiles,{sprintf('Small controls within 1 ms of %.12f s',center), ...
    'All controls strict: RelTol/AbsTol = 1e-9/1e-11 | one parameter changed per comparison', ...
    'Log-axis limits are for display; predictor errors use the independent reference root'}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'ablation');
fprintf('三张严格比较图已保存至 %s。\n',out);
end

function higher(folder,name)
% 只为更高精度结果生成一张独立新图。
file = fullfile(folder,'precision','all.mat'); assert(isfile(file),'更高精度结果尚未保存。');
saved = load(file,'results','ref','c'); list = saved.results; ref = saved.ref; c = saved.c;
assert(numel(list) == 4 && c.tail == 5,'复核需要四组共同尾段结果。');
assert(list{1}.label == "M8 h1e-5" && list{2}.label == "M8a h1e-5" && ...
    list{3}.label == "M8a h1e-6" && list{4}.label == "VFCR-S",'更高精度模型顺序不一致。');
for k = 1:4
    verify(list{k},ref,c,false);
    assert(list{k}.c.rtol == 1e-11 && list{k}.c.atol == 1e-13,'四组必须使用共同的更高精度。');
end
results = {list{4},list{1},list{2},list{3}};
labels = ["VFCR-S","M8 h=1e-5","M8a h=1e-5","M8a h=1e-6"];
colors = [0.43,0.27,0.69;0.04,0.35,0.72;0.68,0.23,0.47;0.44,0.16,0.34];
styles = {'-',':','-','-.'};
out = fullfile(folder,name);
assert(~isfile(fullfile(out,'precision.png')) && ~isfile(fullfile(out,'precision.fig')), ...
    '已有精度图，请指定新的输出文件夹。');
if ~isfolder(out), mkdir(out); end
t = ref.t(:); center = c.peak(2); tail = t >= c.tail;
take = abs(t-center) <= 0.001+1e-12; assert(nnz(take) > 100,'峰中心必须保留密集输出。');
fields = {'xd','g'};
titles = {'Decision error to perturbed PFB root','Primal/dual state error to perturbed PFB root'};
axes = {'||x - x_delta*||','||g - g_delta*||'};
fig = canvas(); tiles = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for j = 1:4
    field = mod(j-1,2)+1; ax = nexttile(tiles); values = valuesof(results,fields{field});
    if j <= 2
        at = t(tail); values = values(tail,:);
    else
        at = 1000*(t(take)-center); values = values(take,:);
    end
    draw(ax,at,values,labels,colors,styles,false);
    line = findall(ax,'Type','line','DisplayName','M8a h=1e-5');
    times = linspace(at(1),at(end),8); ids = zeros(size(times));
    for k = 1:numel(times), [~,ids(k)] = min(abs(at-times(k))); end
    set(line,'Marker','o','MarkerIndices',unique(ids),'MarkerSize',3.5,'MarkerFaceColor','white');
    decorate(ax,titles{field},axes{field});
    if j <= 2
        xlim(ax,[c.tail,c.finish]); xticks(ax,c.tail:c.finish);
        xlabel(ax,'Time t (s)','Interpreter','none');
    else
        local(ax,1);
    end
    if j == 1, key(ax,4); end
end
title(tiles,{'Higher-precision check: same problem, initial state and dense reference grid', ...
    'All curves: RelTol/AbsTol = 1e-11/1e-13 | h labels apply to GNN only', ...
    sprintf('Log-axis limits are for display | tail 5-10 s and peak at %.12f s',center)}, ...
    'FontSize',14,'Interpreter','none');
saveplot(fig,out,'precision');
fprintf('更高精度比较图已保存至 %s。\n',fullfile(out,'precision.png'));
end

function verify(r,ref,c,strict)
% 只比较完整结果和共同的题目设置。
assert(r.complete && string(r.status) == "complete",'结果未完整返回。');
assert(isequal(r.t(:),ref.t(:)) && r.c.delta == c.delta && r.c.tail == c.tail, ...
    '参考时间或统计口径不一致。');
assert(size(r.g,2) == 7 && norm(r.g(1,:).'-c.g0) <= 1e-13,'共同初值或状态维数不一致。');
if strict, assert(r.c.rtol == 1e-9 && r.c.atol == 1e-11,'比较必须使用严格精度。'); end
end

function values = valuesof(results,name)
% 读取实际保存的误差序列。
values = zeros(numel(results{1}.t),numel(results));
for k = 1:numel(results), values(:,k) = results{k}.metrics.(name); end
assert(all(isfinite(values),'all') && all(values >= 0,'all'),'误差必须为有限非负数。');
end

function fig = canvas()
% 使用隐藏图窗和统一字体。
fig = figure('Visible','off','Color','white','Position',[50,50,1500,950], ...
    'DefaultAxesFontName','Arial','DefaultAxesFontSize',11, ...
    'DefaultTextInterpreter','none','DefaultAxesTickLabelInterpreter','none', ...
    'DefaultLegendInterpreter','none');
end

function draw(ax,t,values,names,colors,styles,mark)
% 保持模型颜色并用线型区别参数。
hold(ax,'on');
for k = 1:size(values,2)
    line = plot(ax,t,values(:,k),'Color',colors(k,:),'LineStyle',styles{k}, ...
        'LineWidth',1.4,'DisplayName',names(k));
    if mark && k == 5
        times = linspace(t(1),t(end),8); ids = zeros(size(times));
        for j = 1:numel(times), [~,ids(j)] = min(abs(t-times(j))); end
        set(line,'Marker','o','MarkerIndices',unique(ids),'MarkerSize',3.5,'MarkerFaceColor','white');
    end
end
end

function decorate(ax,name,label)
% 使用对数误差轴和浅网格。
title(ax,name,'Interpreter','none','FontSize',11,'FontWeight','bold');
ylabel(ax,label,'Interpreter','none');
set(ax,'YScale','log','Box','on','GridColor',[0.65,0.65,0.65], ...
    'GridAlpha',0.18,'MinorGridAlpha',0.08,'LineWidth',0.8,'TickDir','out');
grid(ax,'on');
end

function local(ax,limit)
% 将切换窗口标为距峰时间。
xlim(ax,[-limit,limit]);
if limit == 1, xticks(ax,-1:0.5:1); else, xticks(ax,linspace(-limit,limit,7)); end
xline(ax,0,':','Color',[0.55,0.55,0.55],'HandleVisibility','off');
xlabel(ax,'Time from peak (ms)','Interpreter','none');
end

function key(ax,count)
% 将共同图例放在所有面板上方。
label = legend(ax,'Orientation','horizontal','NumColumns',count, ...
    'FontSize',10,'Interpreter','none');
label.Layout.Tile = 'north';
end

function saveplot(fig,folder,name)
% 保存后重开可编辑图并检查图片尺寸。
png = fullfile(folder,[name,'.png']); file = fullfile(folder,[name,'.fig']);
drawnow; state = capture(fig);
exportgraphics(fig,png,'Resolution',150); savefig(fig,file); close(fig);
info = imfinfo(png); assert(info.Width >= 1000 && info.Height >= 700);
check = openfig(file,'invisible'); drawnow;
assert(numel(findall(check,'Type','axes')) >= 4);
assert(isequaln(capture(check),state),'重开 FIG 后的数据或标签改变。'); close(check);
end

function state = capture(fig)
% 核对重开前后的曲线数据和显示属性。
axes = findall(fig,'Type','axes'); names = strings(numel(axes),1);
for k = 1:numel(axes), names(k) = string(axes(k).Title.String); end
[~,ids] = sort(names); axes = axes(ids); state = cell(numel(axes),1);
for k = 1:numel(axes)
    ax = axes(k); lines = findall(ax,'Type','line'); names = strings(numel(lines),1);
    for j = 1:numel(lines), names(j) = string(lines(j).DisplayName); end
    [~,ids] = sort(names); lines = lines(ids); curves = cell(numel(lines),1);
    for j = 1:numel(lines)
        line = lines(j);
        curves{j} = {line.DisplayName,line.XData,line.YData,line.Color, ...
            line.LineStyle,line.LineWidth,line.Marker,line.MarkerIndices};
    end
    state{k} = {ax.Title.String,ax.XLabel.String,ax.YLabel.String, ...
        ax.XLim,ax.YLim,ax.YScale,curves};
end
end
