function list = jobs(c,prior)
% 先验证机制再运行统一噪声和初值实验。
data = load(fullfile(prior,'data','jobs.mat'),'list');
old = data.list;
base = old(find([old.group] == "main" & [old.model] == "GNN" & [old.noise] == "zero",1));
base.variant = "finite"; base.l = 7; base.eps = c.primary.eps;
base.lambda = c.primary.lambda; base.rates = 200; base.reuse = "";
list = repmat(base,0,1);
push("scan","M1",100,0,0,"zero",1,1,2,"finite",base.g0);
id = find([old.group] == "main" & [old.model] == "GNN" & [old.noise] == "zero",1);
list(end).reuse = string(fullfile(prior,'data',sprintf('run_%03d.mat',old(id).id)));
push("scan","VFCR",100,0,0,"zero",1,1,2,"finite",base.g0);
push("scan","VFCR-S",100,1e-2,0,"zero",1,1,2,"finite",base.g0);
for eps = [1e-3,1e-2]
    for gamma = [100,200]
        push("scan","M2",gamma,eps,0,"zero",1,1,2,"finite",base.g0);
        for lambda = [1e-6,1e-4,1e-2]
            push("scan","M3",gamma,eps,lambda,"zero",1,1,2,"finite",base.g0);
        end
    end
    for lambda = [1e-6,1e-4,1e-2]
        for budget = [200,300]
            rates = [100,200,200]*budget/300;
            push("scan","M4",rates,eps,lambda,"zero",1,1,2,"finite",base.g0);
        end
    end
end
for gamma = [150,300]
    for model = ["M2","M3"]
        push("scan",model,gamma,1e-2,1e-6,"zero",1,1,2,"finite",base.g0);
    end
end
models = ["M2","M3","M4","VFCR","VFCR-S"];
for kind = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"]
    seeds = 1;
    if any(kind == ["gamma","gaussian"]), seeds = 1:5; end
    for seed = seeds
        for model = models
            rates = c.primary.gamma;
            if model == "M4", rates = c.primary.rates; end
            push("main",model,rates,c.primary.eps,c.primary.lambda,kind,seed,1,2,"finite",base.g0);
        end
    end
end
for rep = 2:3
    for kind = ["zero","cosine"]
        for model = models
            rates = c.primary.gamma;
            if model == "M4", rates = c.primary.rates; end
            push("timing",model,rates,c.primary.eps,c.primary.lambda,kind,1,rep,2,"finite",base.g0);
        end
    end
end
for model = models
    rates = c.primary.gamma;
    if model == "M4", rates = c.primary.rates; end
    push("precision",model,rates,c.primary.eps,c.primary.lambda,"zero",1,1,2,"finite",base.g0);
    list(end).rtol = 1e-7; list(end).atol = 1e-9;
end
for scale = c.scales
    j = old(find([old.group] == "initial" & [old.model] == "GNN" & [old.scale] == scale,1));
    for model = ["M1",models]
        rates = c.primary.gamma;
        if model == "M1", rates = 100; end
        if model == "M4", rates = c.primary.rates; end
        push("initial",model,rates,c.primary.eps,c.primary.lambda,"zero",1,1,1,"reduced",j.g0(1:3));
        list(end).scale = scale;
    end
end
j = old(find([old.group] == "initial" & [old.model] == "GNN" & [old.scale] == 1,1));
for model = models
    rates = c.primary.gamma;
    if model == "M4", rates = c.primary.rates; end
    push("bounds",model,rates,c.primary.eps,c.primary.lambda,"zero",1,1,1,"large",j.g0);
end
ids = find(ismember([list.model],["M3","M4"]));
for id = ids
    j = list(id); j.eps = 0;
    same = find([list.model] == j.model & [list.group] == j.group ...
        & [list.eps] == 0 & [list.lambda] == j.lambda & abs([list.gamma]-j.gamma) < 1e-8);
    if j.group == "scan" && ~isempty(same), continue, end
    j.id = numel(list)+1; list(end+1,1) = j;
end

    function push(group,model,rates,eps,lambda,noise,seed,rep,ex,variant,g)
        j = base; j.id = numel(list)+1; j.group = group; j.model = model;
        j.rates = rates; j.gamma = norm(rates); j.eps = eps; j.lambda = lambda;
        j.noise = noise; j.seed = seed; j.repeat = rep; j.example = ex; j.variant = variant;
        j.g0 = g; j.l = numel(g); j.p = 1; j.reuse = ""; j.port = "state";
        if startsWith(model,"VFCR"), j.p = 0.5; end
        list(end+1,1) = j;
    end
end
