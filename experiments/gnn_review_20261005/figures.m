function figures()
% 从新结果重建主要比较图。
folder = fileparts(mfilename('fullpath')); out = fullfile(folder,'figures');
if ~isfolder(out), mkdir(out); end
data = load(fullfile(folder,'rerun','data','all.mat'),'results','ref');
models = ["M2","M3","M4","VFCR-S"]; results = cell(4,1);
for k = 1:4
    id = find(cellfun(@(r) r.job.group == "main" && r.job.model == models(k) && r.job.noise == "zero",data.results),1);
    results{k} = data.results{id};
end
fig = figure('Visible','off','Position',[40,40,1350,1050]); tiledlayout(3,2);
for k = 1:6
    nexttile; hold on;
    for j = 1:4
        r = results{j};
        if k <= 2, plot(r.t,r.y(:,k),'LineWidth',1); else
            names = {'e','ed','eg','eq'}; semilogy(r.t,max(r.metrics.(names{k-2}),realmin),'LineWidth',1);
        end
    end
    if k <= 2
        plot(data.ref.t,data.ref.xd(:,k),'k--','LineWidth',1);
        ylabel(sprintf('x%d',k)); legend([models,"delta target"],'Location','best');
    else
        labels = {'Residual','Decision error','Full state error','Equality violation'};
        ylabel(labels{k-2}); set(gca,'YScale','log'); legend(models,'Location','best');
    end
    xlabel('Time (s)'); grid on;
end
exportgraphics(fig,fullfile(out,'main.png'),'Resolution',160);
savefig(fig,fullfile(out,'main.fig')); close(fig);
fig = figure('Visible','off','Position',[40,40,1300,1000]); tiledlayout(3,1);
names = {'grad','sigma','speed'}; labels = {'Gradient norm','Smallest singular value','Control speed'};
for k = 1:3
    nexttile; hold on;
    for j = 1:4, r = results{j}; plot(r.t,r.metrics.(names{k})); end
    if k < 3, set(gca,'YScale','log'); end
    ylabel(labels{k}); xlabel('Time (s)'); legend(models,'Location','best'); grid on;
end
exportgraphics(fig,fullfile(out,'mechanics.png'),'Resolution',160);
savefig(fig,fullfile(out,'mechanics.fig')); close(fig);
data = load(fullfile(folder,'mechanism','dense.mat'),'results');
fig = figure('Visible','off','Position',[50,50,1100,760]); tiledlayout(2,2);
peak = 3.4611091472136+2*pi;
names = {'e','e','eg','ed'}; labels = {'Residual','Residual','Full state error','Decision error'};
for k = 1:4
    nexttile; hold on;
    for j = 1:2, r = data.results{j}; plot(r.t,r.metrics.(names{k})); end
    xlim([peak-0.005,peak+0.01]); ylabel(labels{k}); xlabel('Time (s)'); grid on;
    if k == 2, set(gca,'YScale','log'); end
    legend('gamma=200','gamma=5000','Location','best');
end
exportgraphics(fig,fullfile(folder,'mechanism','dense.png'),'Resolution',160);
savefig(fig,fullfile(folder,'mechanism','dense.fig')); close(fig);
files = dir(fullfile(folder,'**','*.fig'));
for k = 1:numel(files)
    fig = openfig(fullfile(files(k).folder,files(k).name),'invisible');
    assert(~isempty(findall(fig,'Type','axes'))); close(fig);
end
files = dir(fullfile(folder,'**','*.png'));
for k = 1:numel(files), picture = imread(fullfile(files(k).folder,files(k).name)); assert(~isempty(picture)); end
fprintf('已检查 %d 张 PNG，MATLAB FIG 均能重新打开。\n',numel(files));
end
