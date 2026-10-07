function run_review(limit)
% 独立重跑指定实验并逐组保存结果。
if nargin < 1, limit = Inf; end
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(folder)));
v2 = fullfile(root,'experiments','gnn_v2_20261002');
prior = fullfile(root,'experiments','gnn_vfcr_20261002');
data = fullfile(folder,'data'); logs = fullfile(folder,'logs');
if ~isfolder(data), mkdir(data); end
if ~isfolder(logs), mkdir(logs); end
before = pwd; guard = onCleanup(@() cd(before));
addpath(v2,'-begin'); addpath(root,'-begin'); cd(root);
clear config jobs simulate metrics reference problemof noises collect parts motion
clear my_system gnn_system pg_system gsmooth fsmooth digest
names = {'config','jobs','simulate','metrics','reference','problemof','noises', ...
    'collect','parts','motion','pg_system','gsmooth','fsmooth','digest'};
paths = strings(numel(names)+4,2);
for k = 1:numel(names)
    paths(k,:) = [string(names{k}),string(which(names{k}))];
    assert(strcmpi(which(names{k}),fullfile(v2,[names{k},'.m'])));
end
names = {'my_system','gnn_system','example1_problem','example2_problem'};
for k = 1:numel(names)
    paths(k+14,:) = [string(names{k}),string(which(names{k}))];
    assert(strcmpi(which(names{k}),fullfile(root,[names{k},'.m'])));
end
diary(fullfile(logs,'progress.log')); cleanup = onCleanup(@() diary('off'));
c = config(); full = jobs(c,prior);
models = ["M2","M3","M4","VFCR-S"];
ids = [1:35,find([full.group] == "main" & ismember([full.model],models) & [full.eps] == c.primary.eps), ...
    find([full.group] == "precision" & ismember([full.model],models) & [full.eps] == c.primary.eps)];
assert(numel(ids) == 99 && numel(unique(ids)) == 99);
list = full(ids); source = [list.id].';
for k = 1:numel(list), list(k).id = k; list(k).reuse = ""; end
files = [fullfile(root,{'my_system.m','gnn_system.m','example1_problem.m','example2_problem.m'}), ...
    fullfile(v2,{'config.m','jobs.m','simulate.m','metrics.m','reference.m','problemof.m', ...
    'noises.m','collect.m','parts.m','motion.m','pg_system.m','gsmooth.m','fsmooth.m'}), ...
    {mfilename('fullpath')+".m"}];
meta.root = root; meta.paths = paths; meta.files = files;
meta.hash = cellfun(@digest,files,'UniformOutput',false);
meta.matlab = string(version); meta.date = string(datetime('now','TimeZone','Asia/Shanghai'));
meta.scope = "35 scan + 60 main + 4 precision；全部重新求解，reuse 为空。";
file = fullfile(data,'settings.mat');
if isfile(file)
    old = load(file,'meta','ref');
    assert(isequal(old.meta.hash,meta.hash),'代码校验值改变，不能混入同一批重跑结果。');
    ref = old.ref;
else
    ref = reference(c.t,problemof(list(1)),c.delta);
end
assert(max(ref.check,[],'all') < 1e-10);
save(file,'c','list','source','meta','ref');
writetable(array2table(paths,'VariableNames',{'function','path'}), ...
    fullfile(data,'paths.csv'),'Encoding','UTF-8');
old = readtable(fullfile(v2,'data','metrics.csv'),'TextType','string');
results = cell(numel(list),1); done = 0; timer = tic;
for k = 1:numel(list)
    j = list(k); file = fullfile(data,sprintf('run_%03d.mat',k));
    if isfile(file)
        saved = load(file,'r'); r = saved.r; assert(isequaln(r.job,j));
    elseif done < limit
        fprintf('[%s] [%d/%d] source=%d %s %s noise=%s seed=%d gamma=%g eps=%g lambda=%g\n', ...
            string(datetime('now','TimeZone','Asia/Shanghai')),k,numel(list),source(k), ...
            j.group,j.model,j.noise,j.seed,j.gamma,j.eps,j.lambda);
        p = problemof(j); ci = c; ci.l = j.l; bank = noises(j.noise,j.seed,ci);
        r = simulate(j,p,bank,ci); r.metrics = metrics(r,ref,p,ci);
        r.source = source(k); r.fresh = true;
        prior = load(fullfile(v2,'data',sprintf('run_%03d.mat',source(k))),'r');
        assert(isequaln(bank,prior.r.bank));
        both = isfinite(r.y) & isfinite(prior.r.y);
        delta = abs(r.y-prior.r.y); r.compare.trajectory = max(delta(both),[],'all');
        if isempty(r.compare.trajectory), r.compare.trajectory = NaN; end
        r.compare.samples = nnz(both); r.compare.bank = true;
        r.compare.status = r.status == prior.r.status;
        save(file,'r','-v7.3'); done = done+1;
        fprintf('  %s last=%.8g %.3fs calls=%d deltaRMSE=%.6g residual=%.6g motion=%g olddiff=%.6g\n', ...
            r.status,r.last,r.seconds,r.calls,r.metrics.rmsed,r.metrics.residual, ...
            r.metrics.motion.passed,r.compare.trajectory);
    else
        continue
    end
    results{k} = r;
    stats = collect(results);
    stats.source = source(~cellfun(@isempty,results));
    writetable(stats,fullfile(data,'metrics.csv'),'Encoding','UTF-8');
    compare = comparison(results,old,source);
    writetable(compare,fullfile(data,'comparison.csv'),'Encoding','UTF-8');
end
stats = collect(results); stats.source = source(~cellfun(@isempty,results));
compare = comparison(results,old,source); finished = all(~cellfun(@isempty,results));
meta.seconds = toc(timer);
for k = 1:numel(files), assert(strcmp(digest(files{k}),meta.hash{k})); end
meta.unchanged = true;
save(fullfile(data,'all.mat'),'c','list','source','meta','ref','stats','compare','results','finished','-v7.3');
summary = summarize(stats,compare);
writetable(summary,fullfile(data,'summary.csv'),'Encoding','UTF-8');
fprintf('独立重跑已保存 %d/%d 组，完成到 t=10 的有 %d 组，总耗时 %.1f 秒。\n', ...
    height(stats),numel(list),sum(stats.complete),meta.seconds);
disp(summary);
end

function out = comparison(results,old,source)
% 比较同一配置的新旧指标和轨迹。
rows = repmat(struct(),0,1);
fields = {'rmse','rmsed','residual','equality','inequality','full'};
for k = 1:numel(results)
    if isempty(results{k}), continue, end
    r = results{k}; j = r.job; index = find(old.id == source(k),1);
    assert(~isempty(index));
    row = struct('id',k,'source',source(k),'group',j.group,'model',j.model, ...
        'noise',j.noise,'seed',j.seed,'status',r.status,'oldstatus',old.status(index), ...
        'same',r.compare.status,'complete',r.complete,'oldcomplete',old.complete(index), ...
        'seconds',r.seconds,'oldseconds',old.seconds(index),'calls',r.calls, ...
        'oldcalls',old.calls(index),'callsdiff',r.calls-old.calls(index), ...
        'trajectory',r.compare.trajectory,'samples',r.compare.samples);
    for f = fields
        name = f{1}; row.(name) = r.metrics.(name); row.(['old',name]) = old.(name)(index);
        row.([name,'diff']) = abs(row.(name)-row.(['old',name]));
    end
    if isempty(rows), rows = row; else, rows(end+1,1) = row; end
end
if isempty(rows), out = table(); else, out = struct2table(rows); end
end

function out = summarize(stats,compare)
% 汇总各模型的完成率和速度检查。
rows = repmat(struct(),0,1);
for model = unique(stats.model,'stable').'
    take = stats.model == model; tested = take & isfinite(stats.motion);
    row = struct('model',model,'saved',sum(take),'complete',sum(stats.complete(take)), ...
        'failed',sum(~stats.complete(take)),'motiontested',sum(tested), ...
        'motionpassed',sum(stats.motion(tested) == 1),'samestatus',sum(compare.same(take)), ...
        'callsmaxdiff',max(abs(compare.callsdiff(take))), ...
        'trajectorymaxdiff',max(compare.trajectory(take),[],'omitmissing'), ...
        'rmsedmaxdiff',max(compare.rmseddiff(take),[],'omitmissing'), ...
        'seconds',sum(stats.seconds(take)));
    if isempty(rows), rows = row; else, rows(end+1,1) = row; end
end
out = struct2table(rows);
end
