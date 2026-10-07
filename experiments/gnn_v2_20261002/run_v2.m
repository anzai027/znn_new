function run_v2(limit)
% 分组运行新模型并保存全部数据。
if nargin < 1, limit = Inf; end
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
prior = fullfile(fileparts(folder),'gnn_vfcr_20261002');
before = pwd; guard = onCleanup(@() cd(before)); cd(folder); addpath(root,'-end');
clear my_system gnn_system example1_problem example2_problem
assert(strcmpi(which('my_system'),fullfile(root,'my_system.m')));
assert(strcmpi(which('gnn_system'),fullfile(root,'gnn_system.m')));
c = config(); list = jobs(c,prior);
code = {'my_system.m','gnn_system.m','example1_problem.m','example2_problem.m'};
meta.root = root; meta.matlab = string(version); meta.date = string(datetime('now','TimeZone','Asia/Shanghai'));
meta.hash = cellfun(@(x) digest(fullfile(root,x)),code,'UniformOutput',false);
file = fullfile(folder,'data','settings.mat');
if isfile(file)
    old = load(file,'meta'); assert(isequal(old.meta.hash,meta.hash),'原文件在运行之间改变。');
end
meta.chat = "https://chatgpt.com/g/g-p-6ab9211b1de081919276a3f4e721199e-gnngai-jin/c/6a9e7bd7-8730-83ec-9027-9b4f35cd5f38";
meta.primary = "主参数预先设为 lambda=1e-6、eps=1e-2，M2/M3 gamma=200，M4 rates=[100,200,200]。";
meta.solve = "增广最小二乘实现阻尼 Gauss-Newton，等价于正规方程但不显式求逆。";
cases = ["finite","reduced","large"]; refs = cell(3,1);
for k = 1:3
    file = fullfile(folder,'data',"reference_"+cases(k)+".mat");
    if isfile(file)
        data = load(file,'ref'); ref = data.ref;
    else
        j = list(find([list.variant] == cases(k),1)); ref = reference(c.t,problemof(j),c.delta);
        save(file,'ref');
    end
    refs{k} = ref;
end
checks = check(c,refs);
save(fullfile(folder,'data','settings.mat'),'c','list','meta','refs','checks');
results = cell(numel(list),1); done = 0;
for k = 1:numel(list)
    j = list(k); file = fullfile(folder,'data',sprintf('run_%03d.mat',k));
    if isfile(file)
        data = load(file,'r'); r = data.r; assert(isequaln(r.job,j));
    elseif done < limit
        fprintf('[%d/%d] %s %s ex%d %s %s gamma=%g eps=%g lambda=%g\n', ...
            k,numel(list),j.group,j.model,j.example,j.variant,j.noise,j.gamma,j.eps,j.lambda);
        p = problemof(j); ci = c; ci.l = j.l; bank = noises(j.noise,j.seed,ci);
        if j.reuse ~= ""
            data = load(j.reuse,'r'); r = data.r;
            assert(isequal(r.job.g0,j.g0) && isequal(r.bank,bank)); r.job = j;
        else
            r = simulate(j,p,bank,ci);
        end
        ref = refs{find(cases == j.variant,1)}; r.metrics = metrics(r,ref,p,ci);
        save(file,'r','-v7.3'); done = done+1;
        fprintf('  %s last=%.8g %.3fs calls=%d deltaRMSE=%.4g residual=%.4g motion=%g\n', ...
            r.status,r.last,r.seconds,r.calls,r.metrics.rmsed,r.metrics.residual,r.metrics.motion.passed);
    else, continue
    end
    results{k} = r;
    stats = collect(results); writetable(stats,fullfile(folder,'data','metrics.csv'),'Encoding','UTF-8');
end
stats = collect(results); finished = all(~cellfun(@isempty,results));
for k = 1:numel(code), assert(strcmp(digest(fullfile(root,code{k})),meta.hash{k})); end
meta.unchanged = true;
save(fullfile(folder,'data','all.mat'),'c','list','meta','refs','checks','stats','results','finished','-v7.3');
if finished, conclude(results,refs,c,meta,folder); end
fprintf('已保存 %d/%d 组。\n',height(stats),numel(list));
end
