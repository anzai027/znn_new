function bd_joint()
% 检查减小 delta 后相同测量噪声的影响。
folder = bd_setup(); out = fullfile(folder,'joint'); if ~isfolder(out), mkdir(out); end
assert(~isfile(fullfile(out,'all.mat')),'已有联合检查结果。');
p = example2_problem(); s = load(fullfile(folder,'results','plan.mat'),'plans');
ref = load(fullfile(folder,'results','ref_1e-06.mat'),'ref'); stats = table();
diary(fullfile(out,'run.log')); guard = onCleanup(@() diary('off'));
ids = [15,17];
for k = 1:2
    plan = s.plans{ids(k)}; c = plan.c; c.delta = 1e-6;
    fprintf('开始联合检查 %s\n',plan.model);
    r = bd_run(plan.model,p,c); r.label = plan.model+" delta1e-6 noise"; r.group = "joint";
    [r.metrics,row] = bd_metrics(r,ref.ref,p); stats = [stats;row];
    save(fullfile(out,sprintf('run_%02d.mat',k)),'r','row','-v7.3');
    writetable(stats,fullfile(out,'metrics.csv')); disp(row);
end
save(fullfile(out,'all.mat'),'stats');
end
