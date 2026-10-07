function bd_extra()
% 复核积分精度并补充正则化的无噪声基线。
folder = bd_setup(); out = fullfile(folder,'extra'); if ~isfolder(out), mkdir(out); end
assert(~isfile(fullfile(out,'all.mat')),'已有复核结果。');
diary(fullfile(out,'run.log')); guard = onCleanup(@() diary('off'));
p = example2_problem(); s = load(fullfile(folder,'results','plan.mat'),'plans');
ids = [3,6,15,17,15,17,3,7]; models = ["M8a-BDF2","VFCR-FD", ...
    "M8a-BDF2","VFCR-FD","M8a-BDF2","VFCR-FD","M8-reg","VFCR-S"];
base = load(fullfile(folder,'results','ref_0.0001.mat'),'ref'); stats = table();
for k = 1:numel(ids)
    ci = s.plans{ids(k)}.c;
    if k <= 2 || k == 8, ci.rtol = 1e-12; ci.atol = 1e-14; end
    if k == 3 || k == 4, ci.rtol = 1e-9; ci.atol = 1e-11; end
    if k == 5 || k == 6, ci.step = 0.00025; end
    fprintf('开始复核 %d/%d %s\n',k,numel(ids),models(k));
    r = bd_run(models(k),p,ci); r.group = "check"; r.label = models(k)+" check";
    [r.metrics,row] = bd_metrics(r,base.ref,p); row.source = ids(k);
    stats = [stats;row]; save(fullfile(out,sprintf('run_%02d.mat',k)),'r','row','-v7.3');
    writetable(stats,fullfile(out,'metrics.csv'));
    fprintf('完成复核 %d complete=%d xd=%.9g gpeak=%.9g\n',k,r.complete,row.xd,row.gpeak);
end
save(fullfile(out,'all.mat'),'stats');
end
