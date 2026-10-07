function bd_study(name)
% 按精度、扰动量和噪声顺序保存每一组实验。
folder = bd_setup(); if nargin == 0, name = 'data'; end
out = fullfile(folder,name); if ~isfolder(out), mkdir(out); end
assert(~isfile(fullfile(out,'all.mat')),'已有结果，请换目录名。');
diary(fullfile(out,'run.log')); guard = onCleanup(@() diary('off'));
p = example2_problem(); c = bd_config(); plans = {}; stats = table(); refs = {};
models = ["M8a-FD1","M8a-BDF2","M8a-BDF2","M8a-BDF2","VFCR-FD","VFCR-FD","VFCR-S"];
hs = [1e-5,1e-4,1e-5,1e-6,1e-5,1e-5,1e-5]; orders = [1,2,2,2,1,2,2];
for k = 1:7
    ci = c; ci.h = hs(k); ci.order = orders(k);
    plans{end+1} = struct('model',models(k),'group',"zero",'c',ci);
end
for delta = [1e-5,1e-6]
    for model = ["M8a-BDF2","VFCR-FD","VFCR-S"]
        ci = c; ci.delta = delta;
        plans{end+1} = struct('model',model,'group',"delta",'c',ci);
    end
end
for seed = 1:3
    for model = ["M8a-FD1","M8a-BDF2","M8-reg","VFCR-FD"]
        ci = noisy(c,seed,1e-5); ci.order = 2;
        if model == "M8a-FD1", ci.order = 1; end
        plans{end+1} = struct('model',model,'group',"noise",'c',ci);
    end
end
for model = ["M8a-BDF2","M8-reg","VFCR-FD"]
    ci = noisy(c,1,1e-7);
    plans{end+1} = struct('model',model,'group',"amplitude",'c',ci);
end
for h = [1e-4,1e-3]
    for model = ["M8a-BDF2","M8-reg","VFCR-FD"]
        ci = noisy(c,1,1e-5); ci.h = h;
        plans{end+1} = struct('model',model,'group',"step",'c',ci);
    end
end
for model = ["M8a-BDF2","M8-reg","VFCR-FD"]
    ci = noisy(c,1,1e-5); ci.kind = "constant"; ci.step = c.step;
    plans{end+1} = struct('model',model,'group',"constant",'c',ci);
end
paths = ["bd_noise.m","bd_config.m","bd_system.m","bd_run.m", ...
    "bd_metrics.m","bd_study.m","bd_setup.m"];
hash = strings(size(paths));
for k = 1:numel(paths), hash(k) = bd_hash(fullfile(folder,paths(k))); end
save(fullfile(out,'plan.mat'),'plans','paths','hash');
for delta = [1e-4,1e-5,1e-6]
    ref = pp_reference(c.t,p,delta); refs{end+1} = ref;
    save(fullfile(out,sprintf('ref_%g.mat',delta)),'ref','-v7.3');
end
for k = 1:numel(plans)
    plan = plans{k}; ci = plan.c;
    ref = refs{find([1e-4,1e-5,1e-6]==ci.delta)};
    label = sprintf('%s d%.0e h%.0e o%d %s A%.0e seed%d', ...
        plan.model,ci.delta,ci.h,ci.order,ci.kind,ci.amp,ci.seed);
    file = fullfile(out,sprintf('run_%02d.mat',k));
    fprintf('开始 %d/%d %s\n',k,numel(plans),label);
    if isfile(file)
        saved = load(file,'r','row'); r = saved.r; row = saved.row;
        assert(isequaln(r.c,ci) && r.model==plan.model);
    else
        r = bd_run(plan.model,p,ci); r.label = string(label); r.group = plan.group;
        [r.metrics,row] = bd_metrics(r,ref,p);
        save(file,'r','row','hash','-v7.3');
    end
    stats = [stats;row]; writetable(stats,fullfile(out,'metrics.csv'));
    fprintf('完成 %d complete=%d calls=%d xd=%.8g x=%.8g gpeak=%.8g %s\n', ...
        k,r.complete,r.calls,row.xd,row.x,row.gpeak,r.reason);
end
for k = 1:numel(paths), assert(hash(k)==bd_hash(fullfile(folder,paths(k)))); end
save(fullfile(out,'all.mat'),'stats','plans','paths','hash','-v7.3'); disp(stats);
end
function ci = noisy(c,seed,amp)
% 固定每个种子的频率和高斯系数。
ci = c; ci.rtol = 1e-8; ci.atol = 1e-10; ci.step = 0.0005;
ci.seed = seed; ci.amp = amp; ci.kind = "colored";
s = RandStream('mt19937ar','Seed',100+seed);
ci.wave.a = randn(s,7,16)/4; ci.wave.b = randn(s,7,16)/4;
end
