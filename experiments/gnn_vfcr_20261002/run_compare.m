function run_compare(limit)
% 运行或续跑实验，已有结果不会重复计算。
if nargin < 1
    limit = Inf;
end
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
before = pwd;
guard = onCleanup(@() cd(before));
cd(folder);
addpath(root,'-end');
clear my_system gnn_system example1_problem example2_problem
assert(strcmpi(which('my_system'),fullfile(root,'my_system.m')), ...
    'VFCR-ZNN 必须调用原文件。');
assert(strcmpi(which('gnn_system'),fullfile(folder,'gnn_system.m')), ...
    'GNN 必须调用本文件夹内的原样副本。');
c = config();
code = {'my_system.m','gnn_system.m','example1_problem.m','example2_problem.m'};
file = fullfile(folder,'data','meta.mat');
if isfile(file)
    old = load(file,'meta');
    meta = old.meta;
else
    meta.date = string(datetime('now','TimeZone','Asia/Shanghai'));
    meta.matlab = string(version);
    meta.release = string(version('-release'));
    meta.computer = string(computer);
    meta.root = root;
    meta.source = fullfile(root,code);
    meta.hash = cellfun(@digest,meta.source,'UniformOutput',false);
    meta.copy = digest(fullfile(folder,'gnn_system.m'));
    meta.formulas = "原公式，没有平滑，没有修改分母。";
    meta.sign = "以式（58）（59）的 -cos(4t) 为准，正文 P 的展示符号不同。";
    meta.random = "噪声每 0.01 秒更新且分段积分，随机种子与离散化为本实验补充。";
    save(file,'meta','c');
end
for k = 1:numel(code)
    assert(strcmp(digest(fullfile(root,code{k})),meta.hash{k}), ...
        '原文件已改变，请先核对：%s',code{k});
end
assert(strcmp(digest(fullfile(folder,'gnn_system.m')),meta.hash{2}), ...
    'GNN 副本与原文件不同。');
problems = {example1_problem(),example2_problem()};
refs = cell(2,1);
for k = 1:2
    file = fullfile(folder,'data',sprintf('reference%d.mat',k));
    if isfile(file)
        data = load(file,'ref');
        refs{k} = data.ref;
    else
        ref = reference(c.t,problems{k},c.delta);
        refs{k} = ref;
        save(file,'ref');
    end
end
checks = verify(c,problems,refs);
save(fullfile(folder,'data','checks.mat'),'checks');
fprintf('参考解和无导数调用检查通过，雅可比相对误差 %.3e。\n',checks.jacobian);
list = jobs(c);
save(fullfile(folder,'data','jobs.mat'),'list','c');
results = cell(numel(list),1);
done = 0;
for k = 1:numel(list)
    job = list(k);
    file = fullfile(folder,'data',sprintf('run_%03d.mat',k));
    if isfile(file)
        data = load(file,'r');
        r = data.r;
        assert(isequaln(r.job,job),'已有运行的设置与当前设置不同。');
        diagnostic = motion(r,c);
        if ~isfield(r.metrics,'motion') || ~isequaln(r.metrics.motion,diagnostic)
            r.metrics.motion = diagnostic;
            save(file,'r','-v7.3');
        end
    elseif done < limit
        fprintf('[%d/%d] %s ex%d %s %s seed%d gamma=%g p=%g\n', ...
            k,numel(list),job.group,job.example,job.model,job.noise, ...
            job.seed,job.gamma,job.p);
        bank = noises(job.noise,job.seed,c);
        r = simulate(job,problems{job.example},bank,c);
        r.metrics = metrics(r,refs{job.example},problems{job.example},c);
        save(file,'r','-v7.3');
        done = done+1;
        fprintf('  %s: t=%.8g, %.3fs, calls=%d, residual=%.3e\n', ...
            r.status,r.last,r.seconds,r.calls,r.metrics.last);
    else
        continue
    end
    results{k} = r;
    stats = collect(results);
    writetable(stats,fullfile(folder,'data','metrics.csv'),'Encoding','UTF-8');
    fid = fopen(fullfile(folder,'logs','progress.txt'),'w','n','UTF-8');
    fprintf(fid,'%d/%d 已保存；最后 id=%d，状态=%s，最后时间=%.9g。\n', ...
        height(stats),numel(list),k,r.status,r.last);
    fclose(fid);
end
stats = collect(results);
finished = all(~cellfun(@isempty,results));
for k = 1:numel(code)
    assert(strcmp(digest(fullfile(root,code{k})),meta.hash{k}), ...
        '运行后原文件校验失败。');
end
meta.unchanged = true;
save(fullfile(folder,'data','all.mat'),'c','meta','list','refs', ...
    'checks','results','stats','finished','-v7.3');
if finished
    plots(results,refs,c,folder);
    audit = inspect(folder,results,c);
    save(fullfile(folder,'data','audit.mat'),'audit');
    report(stats,refs,c,meta,audit,folder);
    save(fullfile(folder,'data','all.mat'),'c','meta','list','refs', ...
        'checks','results','stats','finished','audit','-v7.3');
    fprintf('全部 %d 组已尝试并保存，完整运行 %d 组。\n',height(stats),sum(stats.complete));
else
    fprintf('本批次保存 %d/%d 组，调用 run_compare 可继续。\n',height(stats),numel(list));
end
end
