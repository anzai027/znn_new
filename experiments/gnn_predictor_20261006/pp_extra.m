function pp_extra()
% 用两个对照区分阻尼和尺度的影响。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(root,'-begin'); addpath(folder,'-begin');
source = load(fullfile(folder,'data','all.mat'),'c','ref','hash');
c = source.c; ref = source.ref; hash = source.hash;
c.rtol = 1e-9; c.atol = 1e-11;
p = example2_problem();
out = fullfile(folder,'extra');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'run.log')); clean = onCleanup(@() diary('off'));
results = cell(2,1); stats = table();
for k = 1:2
    ci = c;
    if k == 1
        group = "kappa01"; model = "M6"; ci.kappa = 0.01;
    else
        group = "unit"; model = "M7"; ci.scale = [1,1,1];
    end
    r = pp_run(model,p,ci);
    assert(r.complete,'额外对照未完成：%s',r.reason);
    [r.metrics,row] = pp_metrics(r,ref,p);
    row = addvars(row,group,'Before',1);
    results{k} = r;
    stats = [stats;row];
    save(fullfile(out,sprintf('run_%02d.mat',k)),'r','hash','-v7.3');
    writetable(stats,fullfile(out,'metrics.csv'));
    fprintf('%s xd=%.9g gpeak=%.9g residual=%.9g calls=%d\n', ...
        group,row.xd,row.gpeak,row.residual,r.calls);
end
save(fullfile(out,'all.mat'),'results','stats','ref','c','hash','-v7.3');
disp(stats);
end
