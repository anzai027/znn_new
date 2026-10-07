function conclude(results,refs,c,meta,folder)
% 生成可复查的曲线和中文结论。
stats = collect(results);
ids = find(stats.group == "main" & stats.noise == "zero");
models = ["M2","M3","M4","VFCR","VFCR-S"];
labels = ["old smooth GNN","preconditioned GNN","block GNN","original VFCR","smooth VFCR"];
for kind = ["zero","cosine"]
    ids = find(stats.group == "main" & stats.noise == kind & stats.eps > 0);
    fig = makefig(1400,780); ax = gobjects(4,1);
    for j = 1:4, ax(j) = subplot(2,2,j); hold(ax(j),'on'); end
    for model = models
        id = ids(find(stats.model(ids) == model,1)); r = results{stats.id(id)};
        label = model; if ~r.complete, label = label+" (interrupted)"; end
        info = struct('last',r.last,'complete',r.complete);
        values = {r.metrics.e,r.metrics.ed,r.metrics.eg,r.metrics.grad};
        for j = 1:4
            semilogy(ax(j),r.t,positive(values{j}),'LineWidth',1.05,'DisplayName',label,'UserData',info);
        end
    end
    yline(ax(1),1e-4,'k:','DisplayName','threshold 1e-4');
    names = ["KKT/PFB residual","error to delta x","error to delta g","gradient norm"];
    for j = 1:4, decorate(ax(j),names(j)); set(ax(j),'YScale','log'); end
    sgtitle("Example 2 | "+kind+" | primary parameters");
    savepic(fig,"main_"+kind);
    fig = makefig(1250,470);
    for j = 1:2
        ax = subplot(1,2,j); hold(ax,'on');
        plot(ax,refs{1}.t,refs{1}.xd(:,j),'k--','LineWidth',1.4,'DisplayName','delta reference');
        for model = ["M2","M3","M4"]
            id = ids(find(stats.model(ids) == model,1)); r = results{stats.id(id)};
            plot(ax,r.t,r.y(:,j),'LineWidth',1,'DisplayName',model);
        end
        decorate(ax,"x_"+j);
    end
    sgtitle("Example 2 states | "+kind); savepic(fig,"states_"+kind);
end
fig = makefig(1400,780);
for k = 1:4
    ax = subplot(2,2,k); hold(ax,'on');
    ids = find(stats.group == "scan" & stats.eps == 1e-2 & (stats.lambda == 1e-6 | stats.model == "M2") ...
        & abs(stats.gamma-200) < 1e-8 & ismember(stats.model,["M2","M3","M4"]));
    for id = ids.'
        r = results{stats.id(id)};
        values = {r.metrics.sigma,r.metrics.condition,r.metrics.angle,r.metrics.direction};
        y = values{k}; if k ~= 3, y = positive(y); end
        plot(ax,r.t,y,'LineWidth',1,'DisplayName',r.job.model);
    end
    names = ["minimum singular value","Jacobian condition number","direction alignment cosine","direction norm"];
    decorate(ax,names(k)); if k ~= 3, set(ax,'YScale','log'); end
end
sgtitle('Same speed cap 200 | same eps 0.01 | conditioning diagnostics');
savepic(fig,"conditioning");
fig = makefig(1250,470);
for k = 1:2
    ax = subplot(1,2,k); hold(ax,'on');
    for model = ["M2","M3"]
        ids = find(stats.group == "scan" & stats.eps == 1e-2 ...
            & (stats.lambda == 1e-6 | stats.model == "M2") & stats.model == model);
        [gamma,order] = sort(stats.gamma(ids)); ids = ids(order);
        if k == 1, value = stats.rmsed(ids); else, value = stats.residual(ids); end
        semilogy(ax,gamma,value,'o-','LineWidth',1.2,'DisplayName',model);
    end
    xlabel(ax,'gamma'); names = ["delta x RMSE (5-10 s)","maximum residual (5-10 s)"];
    ylabel(ax,names(k)); grid(ax,'on'); legend(ax,'Location','best');
    set(ax,'YScale','log');
end
sgtitle('Gain diagnostic | eps=0.01 | lambda=1e-6 for M3'); savepic(fig,"gains");
fig = makefig(1250,470);
for k = 1:2
    ax = subplot(1,2,k); hold(ax,'on'); eps = [1e-3,1e-2]; eps = eps(k);
    for model = ["M3","M4"]
        for gamma = [200,300]
            ids = find(stats.group == "scan" & stats.eps == eps & abs(stats.gamma-gamma) < 1e-8 ...
                & stats.model == model);
            if isempty(ids), continue; end
            [lambda,order] = sort(stats.lambda(ids)); ids = ids(order);
            loglog(ax,lambda,stats.rmsed(ids),'o-','LineWidth',1.1, ...
                'DisplayName',sprintf('%s cap=%g',model,gamma));
        end
    end
    set(ax,'XScale','log','YScale','log'); xlabel(ax,'lambda'); ylabel(ax,'delta x RMSE (5-10 s)');
    grid(ax,'on'); legend(ax,'Location','best'); title(ax,sprintf('eps=%g',eps));
end
sgtitle('Damping and block gains | all tested settings retained'); savepic(fig,"damping");
speeds = speedof(refs{1},folder,c);
fig = makefig(1400,780);
for k = 1:4
    ax = subplot(2,2,k); hold(ax,'on');
    if k == 1
        for j = 1:3
            plot(ax,speeds.t{j},speeds.values{j}(:,1),'DisplayName',sprintf('reference dt=%g',speeds.steps(j)));
        end
        yline(ax,100,'r--','DisplayName','cap 100'); yline(ax,200,'b--','DisplayName','cap 200');
        yline(ax,300,'k--','DisplayName','cap 300');
        ylabel0 = 'reference full state speed';
    elseif k == 2
        plot(ax,speeds.t{3},speeds.values{3}(:,2),'DisplayName','x');
        plot(ax,speeds.t{3},speeds.values{3}(:,3),'DisplayName','mu1');
        plot(ax,speeds.t{3},speeds.values{3}(:,4),'DisplayName','mu2'); ylabel0 = 'reference block speed (dt=0.001)';
    else
        for model = ["M2","M3","M4"]
            id = find(stats.group == "main" & stats.noise == "zero" & stats.model == model,1);
            r = results{stats.id(id)};
            y = r.metrics.speed; if k == 4, y = r.metrics.sigma; end
            plot(ax,r.t,y,'DisplayName',model);
        end
        ylabel0 = 'control speed'; if k == 4, ylabel0 = 'minimum singular value'; end
    end
    decorate(ax,ylabel0);
end
sgtitle('Moving-root speed check | delta reference | control excludes state noise');
savepic(fig,"speeds");
for kind = ["gamma","gaussian"]
    fig = makefig(1400,780);
    for j = 1:4
        ax = subplot(2,2,j); hold(ax,'on'); model = models(j);
        ids = find(stats.group == "main" & stats.noise == kind & stats.model == model & stats.eps > 0);
        for id = ids.'
            r = results{stats.id(id)};
            semilogy(ax,r.t,positive(r.metrics.ed),'DisplayName',"seed "+r.job.seed);
        end
        set(ax,'YScale','log'); decorate(ax,'error to delta x'); title(ax,model);
    end
    sgtitle(kind+" state noise | shared seeds 1-5"); savepic(fig,kind);
end
fig = makefig(1400,780);
for j = 1:4
    ax = subplot(2,2,j); hold(ax,'on'); model = ["M2","M3","M4","VFCR-S"]; model = model(j);
    ids = find(stats.group == "initial" & stats.model == model & stats.eps > 0);
    for id = ids.'
        r = results{stats.id(id)}; label = sprintf('range [0,%g]',r.job.scale);
        if ~r.complete, label = [label,' (interrupted)']; end
        info = struct('last',r.last,'complete',r.complete);
        semilogy(ax,r.t,positive(r.metrics.e),'DisplayName',label,'UserData',info);
    end
    set(ax,'YScale','log'); decorate(ax,'equality KKT residual'); title(ax,model);
end
sgtitle('Example 1 | remove absent inequalities | shared initial x and mu1'); savepic(fig,"initial");
fig = makefig(1250,470);
for j = 1:2
    ax = subplot(1,2,j); hold(ax,'on');
    for model = ["M2","M3","M4","VFCR-S"]
        ids = find((stats.group == "initial" & stats.scale == 1 | stats.group == "bounds") ...
            & stats.model == model & stats.eps > 0);
        for id = ids.'
            r = results{stats.id(id)}; label = r.job.model+" "+r.job.variant;
            if ~r.complete, label = label+" (interrupted)"; end
            info = struct('last',r.last,'complete',r.complete);
            value = r.metrics.ex; if j == 2, value = r.metrics.e; end
            semilogy(ax,r.t,positive(value),'DisplayName',label,'UserData',info);
        end
    end
    set(ax,'YScale','log'); names = ["true x error","KKT/PFB residual"];
    decorate(ax,names(j));
end
sgtitle('Example 1 | 1e8 bounds versus removed constraints | separate state spaces'); savepic(fig,"bounds");
for kind = ["zero","cosine"]
    fig = makefig(1400,780); ax = gobjects(4,1);
    for j = 1:4, ax(j) = subplot(2,2,j); hold(ax(j),'on'); end
    for model = ["M3","M4"]
        for eps = [0,1e-2]
            id = find(stats.group == "main" & stats.model == model & stats.noise == kind & stats.eps == eps,1);
            r = results{stats.id(id)}; label = sprintf('%s eps=%g',model,eps);
            if ~r.complete, label = [label,' (interrupted)']; end
            if r.metrics.motion.passed == 0, label = [label,' (speed failed)']; end
            info = struct('last',r.last,'complete',r.complete);
            [t,values] = early(r,c);
            for j = 1:4
                style = '-'; if eps == 0, style = 'o-'; end
                semilogy(ax(j),t,positive(values{j}),style,'MarkerSize',3, ...
                    'MarkerIndices',unique(round(linspace(1,numel(t),10))),'DisplayName',label,'UserData',info);
            end
        end
    end
    names = ["KKT/PFB residual","error to delta x","error to delta g","control speed"];
    for j = 1:4, decorate(ax(j),names(j)); set(ax(j),'YScale','log'); xlim(ax(j),[0,0.02]); end
    sgtitle("Smooth versus unsmoothed | "+kind+" | accepted points, first 0.02 s");
    savepic(fig,"plain_"+kind);
end
fig = makefig(1400,780);
for j = 1:4
    ax = subplot(2,2,j); hold(ax,'on');
    model = ["M3","M4","M3","M4"]; model = model(j);
    eps = [1e-2,1e-2,0,0]; eps = eps(j);
    ids = find(stats.group == "initial" & stats.model == model & stats.eps == eps);
    for id = ids.'
        r = results{stats.id(id)}; label = sprintf('range [0,%g]',r.job.scale);
        if ~r.complete, label = [label,' (interrupted)']; end
        if r.metrics.motion.passed == 0, label = [label,' (speed failed)']; end
        info = struct('last',r.last,'complete',r.complete);
        semilogy(ax,r.t,positive(r.metrics.e),'DisplayName',label,'UserData',info);
    end
    decorate(ax,'equality KKT residual'); set(ax,'YScale','log'); title(ax,sprintf('%s eps=%g',model,eps));
end
sgtitle('Example 1 | smooth and unsmoothed | all four initial ranges');
savepic(fig,"plain_initial");
audit = auditof(results,c,meta,folder);
save(fullfile(folder,'data','audit.mat'),'audit');
save(fullfile(folder,'data','all.mat'),'speeds','audit','-append');
report(stats,results,c,speeds,audit,folder);

    function fig = makefig(width,height)
        fig = figure('Visible','off','Color','w','Position',[40,40,width,height]);
    end
    function decorate(ax,label)
        xlabel(ax,'t (s)'); ylabel(ax,label); xlim(ax,[0,10]); grid(ax,'on');
        legend(ax,'Location','best','Interpreter','none');
    end
    function savepic(fig,name)
        savefig(fig,fullfile(folder,'figures',name+".fig"));
        exportgraphics(fig,fullfile(folder,'figures',name+".png"),'Resolution',140);
        close(fig);
    end
end

function y = positive(y)
% 缺失数据保留为缺失。
y(isfinite(y) & y < 1e-16) = 1e-16;
end

function [t,values] = early(r,c)
% 用已接受点显示中断前的短轨迹。
ids = find(r.raw.t <= 0.02);
ids = ids(unique([1:min(50,numel(ids)),round(linspace(1,numel(ids),250))]));
t = r.raw.t(ids); y = r.raw.y(ids,1:r.job.l);
p = problemof(r.job); ref = reference(t,p,c.delta);
e = zeros(numel(t),1); speed = e;
for k = 1:numel(t)
    g = y(k,:).'; e(k) = norm(parts(t(k),g,p,c.delta));
    dy = pg_system(t(k),g,r.job.l,p,r.job.rates,r.job.lambda,r.job.eps,c.delta,r.job.model == "M4");
    speed(k) = norm(dy);
end
values = {e,vecnorm(y(:,1:2)-ref.xd,2,2),vecnorm(y-ref.g,2,2),speed};
end

function speeds = speedof(ref,folder,c)
steps = [0.005,0.0025,0.001]; speeds.steps = steps; speeds.peaks = zeros(3,4);
speeds.secant = zeros(3,1);
for k = 1:3
    dt = steps(k);
    r = ref;
    if k ~= 1, r = reference((0:dt:10).',example2_problem(),c.delta); end
    dy = zeros(size(r.g));
    for j = 1:size(dy,2), dy(:,j) = gradient(r.g(:,j),dt); end
    values = [vecnorm(dy,2,2),vecnorm(dy(:,1:2),2,2),abs(dy(:,3)),vecnorm(dy(:,4:7),2,2)];
    speeds.peaks(k,:) = max(values); speeds.secant(k) = max(vecnorm(diff(r.g),2,2)/dt);
    speeds.t{k} = r.t; speeds.values{k} = values;
    if k == 1, speeds.grid = values; end
end
save(fullfile(folder,'data','speeds.mat'),'speeds');
end

function audit = auditof(results,c,meta,folder)
files = dir(fullfile(folder,'figures','*.fig')); audit.names = strings(numel(files),1);
for k = 1:numel(files)
    fig = openfig(fullfile(files(k).folder,files(k).name),'invisible');
    axes = findall(fig,'Type','axes');
    for ax = axes.'
        assert(~isempty(ax.XLabel.String) && ~isempty(ax.YLabel.String));
        assert(~isempty(findall(ax,'Type','line')));
    end
    assert(~isempty(findall(fig,'Type','legend')));
    lines = findall(fig,'Type','line');
    for line = lines.'
        info = line.UserData;
        if isstruct(info) && isfield(info,'complete') && ~info.complete
            assert(all(isnan(line.YData(line.XData > info.last+1e-10))));
        end
    end
    close(fig); [~,name] = fileparts(files(k).name);
    pic = imread(fullfile(folder,'figures',[name,'.png']));
    assert(~isempty(pic) && max(pic(:)) > min(pic(:))); audit.names(k) = string(name);
end
for k = 1:numel(results)
    r = results{k}; stored = load(fullfile(folder,'data',sprintf('run_%03d.mat',k)),'r');
    assert(isequaln(stored.r.job,r.job)); assert(all(isfinite(r.raw.y),'all'));
    assert(all(diff(r.raw.t) > 0)); assert(isequaln(r.metrics.motion,motion(r,c)));
    if ~r.complete, assert(all(isnan(r.y(c.t > r.last+1e-10,:)),'all')); end
end
stats = collect(results);
for kind = unique(stats.noise(stats.group == "main")).'
    for seed = unique(stats.seed(stats.group == "main" & stats.noise == kind)).'
        ids = find(stats.group == "main" & stats.noise == kind & stats.seed == seed);
        for id = ids.'
            assert(isequal(results{stats.id(id)}.bank,results{stats.id(ids(1))}.bank));
            assert(isequal(results{stats.id(id)}.job.g0,results{stats.id(ids(1))}.job.g0));
        end
    end
end
code = {'my_system.m','gnn_system.m','example1_problem.m','example2_problem.m'};
for k = 1:4, assert(strcmp(digest(fullfile(meta.root,code{k})),meta.hash{k})); end
audit.noise = true; audit.originals = true; audit.passed = true;
end
