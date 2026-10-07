function dense()
% 加密活动约束切换附近的输出。
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
v2 = fullfile(root,'experiments','gnn_v2_20261002');
prior = fullfile(root,'experiments','gnn_vfcr_20261002');
cd(root); addpath(v2,'-begin'); addpath(root,'-begin');
out = fullfile(folder,'mechanism'); c = config();
peak = 3.4611091472136;
t = c.t;
for time = [peak,peak+2*pi]
    t = [t;(time-0.005:0.000001:time+0.005).';time];
end
c.t = unique(t); p = example2_problem();
ref = reference(c.t,p,c.delta); list = jobs(c,prior);
base = list(find([list.group] == "main" & [list.model] == "M3" & [list.noise] == "zero" & [list.eps] > 0,1));
rows = table(); results = cell(2,1);
for k = 1:2
    gains = [200,5000]; j = base; j.gamma = gains(k); j.rates = gains(k); j.id = k;
    bank = noises("zero",1,c); a = simulate(j,p,bank,c);
    a.metrics = metrics(a,ref,p,c);
    tail = a.t >= c.tail; t = a.t(tail); m = a.metrics;
    xd = sqrt(trapz(t,m.ed(tail).^2)/(t(end)-t(1)));
    full = sqrt(trapz(t,m.eg(tail).^2)/(t(end)-t(1)));
    row = table(j.gamma,a.complete,a.calls,xd,full,max(m.ed(tail)), ...
        max(m.eg(tail)),max(m.e(tail)),m.tc(2), ...
        'VariableNames',{'gamma','complete','calls','weighted','full','xpeak','gpeak','residual','tc4'});
    rows = [rows;row]; results{k} = a;
    save(fullfile(out,sprintf('dense_%d.mat',k)),'a','ref','-v7.3');
end
writetable(rows,fullfile(out,'dense.csv')); disp(rows);
save(fullfile(out,'dense.mat'),'rows','results','ref','-v7.3');
fig = figure('Visible','off','Position',[50,50,1100,760]); tiledlayout(2,2);
for k = 1:2
    a = results{k}; nexttile(k); hold on;
    for j = 1:2, b = results{j}; plot(b.t,b.metrics.e); end
    xlim([peak+2*pi-0.005,peak+2*pi+0.01]); ylabel('Residual'); grid on;
    if k == 2, set(gca,'YScale','log'); end
    legend('gamma=200','gamma=5000');
end
nexttile; hold on;
for j = 1:2, b = results{j}; plot(b.t,b.metrics.eg); end
xlim([peak+2*pi-0.005,peak+2*pi+0.01]); xlabel('Time (s)'); ylabel('Full state error'); grid on;
nexttile; hold on;
for j = 1:2, b = results{j}; plot(b.t,b.metrics.ed); end
xlim([peak+2*pi-0.005,peak+2*pi+0.01]); xlabel('Time (s)'); ylabel('Decision error'); grid on;
exportgraphics(fig,fullfile(out,'dense.png'),'Resolution',160);
savefig(fig,fullfile(out,'dense.fig')); close(fig);
end
