function run_smooth(limit)
% 续跑平滑试验并保存每次结果。
if nargin < 1
    limit = Inf;
end
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
prior = fullfile(fileparts(folder),'gnn_vfcr_20261002');
before = pwd;
guard = onCleanup(@() cd(before));
cd(folder);
addpath(root,'-end');
clear my_system example1_problem example2_problem
assert(strcmpi(which('my_system'),fullfile(root,'my_system.m')));
c = config();
code = {'my_system.m','gnn_system.m','example1_problem.m','example2_problem.m'};
old = load(fullfile(prior,'data','meta.mat'),'meta');
for k = 1:numel(code)
    assert(strcmp(digest(fullfile(root,code{k})),old.meta.hash{k}));
end
data = load(fullfile(prior,'data','jobs.mat'),'list');
list = series(data.list);
problems = {example1_problem(),example2_problem()};
refs = cell(2,1);
for k = 1:2
    data = load(fullfile(prior,'data',sprintf('reference%d.mat',k)),'ref');
    refs{k} = data.ref;
end
meta.date = string(datetime('now','TimeZone','Asia/Shanghai'));
meta.matlab = string(version);
meta.root = root;
meta.prior = prior;
meta.hash = old.meta.hash;
meta.gnn = "-gamma*s/(sqrt(norm(s)^2+eps^2))^p，gamma=100，p=1。";
meta.vfcr = "仅 abs(x)<eps 内用奇三次函数连接原激活函数，边界值和一阶导数连续。";
meta.selection = "1e-3 为预先指定尺度；扫描后发现其无噪声 GNN 未完成，因此追加验证 1e-2，属于依据扫描选择的探索结果。";
meta.units = "GNN 尺度作用于梯度范数，VFCR 尺度作用于残差分量，数值相同不表示效果相同。";
checks = check(c,problems);
save(fullfile(folder,'data','settings.mat'),'c','meta','list','checks','refs');
results = cell(numel(list),1);
done = 0;
for k = 1:numel(list)
    job = list(k);
    file = fullfile(folder,'data',sprintf('run_%03d.mat',k));
    if isfile(file)
        data = load(file,'r');
        r = data.r;
        assert(isequaln(r.job,job));
    elseif done < limit
        fprintf('[%d/%d] %s ex%d %s %s eps=%g repeat=%d seed=%d\n', ...
            k,numel(list),job.group,job.example,job.model,job.noise, ...
            job.eps,job.repeat,job.seed);
        data = load(fullfile(prior,'data',sprintf('run_%03d.mat',job.base)),'r');
        bank = data.r.bank;
        assert(isequal(bank,noises(job.noise,job.seed,c)));
        assert(isequal(data.r.job.g0,job.g0));
        r = simulate(job,problems{job.example},bank,c);
        r.metrics = metrics(r,refs{job.example},problems{job.example},c);
        save(file,'r','-v7.3');
        done = done+1;
        fprintf('  %s t=%.9g %.3fs calls=%d deltaRMSE=%.4g motion=%g\n', ...
            r.status,r.last,r.seconds,r.calls,r.metrics.rmsed,r.metrics.motion.passed);
    else
        continue
    end
    results{k} = r;
    stats = tableof(results);
    writetable(stats,fullfile(folder,'data','metrics.csv'),'Encoding','UTF-8');
end
stats = tableof(results);
finished = all(~cellfun(@isempty,results));
for k = 1:numel(code)
    assert(strcmp(digest(fullfile(root,code{k})),meta.hash{k}));
end
meta.unchanged = true;
save(fullfile(folder,'data','all.mat'),'c','meta','list','checks','refs', ...
    'results','stats','finished','-v7.3');
if finished
    conclude(results,refs,c,meta,folder);
end
fprintf('本次已保存 %d/%d 组。\n',height(stats),numel(list));
end

function stats = tableof(results)
stats = collect(results);
ids = find(~cellfun(@isempty,results));
stats.eps = cellfun(@(r) r.job.eps,results(ids));
stats.base = cellfun(@(r) r.job.base,results(ids));
end

function checks = check(c,problems)
v = c.vfcr;
checks.connection = zeros(4,2);
for k = 1:4
    eps = 10^(k-5);
    h = eps*1e-5;
    exact = v.a*exp(eps^v.q)*eps^v.p;
    slope = exact*(v.p+v.q*eps^v.q)/eps;
    checks.connection(k,1) = abs(fsmooth(eps,v,eps)-exact);
    left = (fsmooth(eps,v,eps)-fsmooth(eps-h,v,eps))/h;
    right = (fsmooth(eps+h,v,eps)-fsmooth(eps,v,eps))/h;
    checks.connection(k,2) = max(abs([left,right]-slope))/slope;
    assert(checks.connection(k,1) < 1e-12 && checks.connection(k,2) < 1e-4);
    assert(fsmooth(0,v,eps) == 0);
end
problem = rmfield(problems{2},{'dG','dh','dP','du','dQ','dv'});
g = ones(c.l,1)/2;
dy = gsmooth(0,g,c.l,problem,100,1,c.delta,1e-3);
assert(all(isfinite(dy)) && norm(dy) <= 100*(1+1e-12));
checks.derivatives = true;
checks.speed = norm(dy);
end
