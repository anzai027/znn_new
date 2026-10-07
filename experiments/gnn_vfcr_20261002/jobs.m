function list = jobs(c)
% 排列全部实验，并固定共享初值。
base = struct('id',0,'group',"",'example',0,'model',"",'port',"", ...
    'noise',"",'seed',0,'scale',0,'g0',zeros(c.l,1),'gamma',0, ...
    'p',0,'repeat',0,'rtol',0,'atol',0,'step',0);
list = repmat(base,0,1);
stream = RandStream('mt19937ar','Seed',c.seed);
for scale = c.scales
    g = scale*rand(stream,c.l,1);
    for model = ["VFCR-ZNN","GNN"]
        push("initial",1,model,"none","zero",1,scale,g,100,1,1);
    end
end
stream = RandStream('mt19937ar','Seed',c.seed);
g = rand(stream,c.l,1);
kinds = ["zero","constant","linear","cosine","gamma","gaussian","harmonic"];
for kind = kinds
    seeds = 1;
    if any(kind == ["gamma","gaussian"])
        seeds = c.seeds;
    end
    for seed = seeds
        push("paper",2,"VFCR-ZNN","paper",kind,seed,1,g,100,1,1);
        for model = ["VFCR-ZNN","GNN"]
            push("main",2,model,"state",kind,seed,1,g,100,1,1);
        end
    end
end
for k = 1:size(c.pairs,1)
    for kind = ["zero","cosine"]
        push("parameter",2,"GNN","state",kind,1,1,g,c.pairs(k,1),c.pairs(k,2),1);
    end
end
for model = ["VFCR-ZNN","GNN"]
    push("precision",2,model,"state","zero",1,1,g,100,1,1);
    list(end).rtol = 1e-7;
    list(end).atol = 1e-9;
end
for rep = 2:3
    for kind = ["zero","cosine"]
        for model = ["VFCR-ZNN","GNN"]
            push("timing",2,model,"state",kind,1,1,g,100,1,rep);
        end
    end
end

    function push(group,ex,model,port,kind,seed,scale,g,gamma,p,rep)
        if model == "VFCR-ZNN"
            gamma = NaN;
            p = c.vfcr.p;
        end
        job = struct('id',numel(list)+1,'group',group,'example',ex, ...
            'model',model,'port',port,'noise',kind,'seed',seed, ...
            'scale',scale,'g0',g,'gamma',gamma,'p',p,'repeat',rep, ...
            'rtol',c.rtol,'atol',c.atol,'step',c.step);
        list(end+1,1) = job;
    end
end
