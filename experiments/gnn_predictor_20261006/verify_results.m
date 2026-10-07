function checks = verify_results(name,stage)
% 独立验收已有结果，不再运行积分器。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(root,'-begin'); addpath(folder,'-begin');
if nargin == 0, name = 'data'; end
source = fullfile(folder,name);
out = fullfile(folder,'verification');
if ~isfolder(out), mkdir(out); end
if nargin > 1 && string(stage) == "refine"
    checks = refinement(folder,source,out);
    return
end
diary(fullfile(out,'run.log')); clean = onCleanup(@() diary('off'));
data = load(fullfile(source,'all.mat'));
results = data.results; stats = data.stats; ref = data.ref; c = data.c;
assert(numel(results) == 8 && height(stats) == 8,'主验收需要八组结果。');
assert(all(cellfun(@(r) r.complete,results)),'存在未完成结果，不能验收完整性能。');
checks.source = string(source);
checks.hash = true;
files = fieldnames(data.hash);
for k = 1:numel(files)
    file = files{k};
    if strcmp(file,'problem'), path = fullfile(root,'example2_problem.m');
    else, path = fullfile(folder,[file,'.m']); end
    assert(strcmp(digest(path),data.hash.(file)),'源码与运行哈希不同：%s',file);
end
p = example2_problem();
fresh = pp_reference(c.t,p,c.delta);
checks.reference = max(vecnorm(fresh.g-ref.g,2,2));
checks.residual = max(fresh.error);
checks.velocity = max(vecnorm(fresh.velocity-ref.velocity,2,2));
checks.derivative = max(fresh.derivative);
assert(checks.reference < 1e-10 && checks.residual < 1e-9 && checks.velocity < 1e-8);
assert(checks.derivative < 1e-9);
coarse = (c.start:0.005:c.finish).';
ids = round(interp1(c.t,(1:numel(c.t)).',coarse,'nearest'));
assert(max(abs(c.t(ids)-coarse)) < 1e-11,'输出网格未包含稀疏对照点。');
assert(isequal(ref.t,c.t));
models = ["M3","M5","M6","M7"];
rows = repmat(struct(),0,1);
sampling = repmat(struct(),0,1);
limits = repmat(struct(),0,1);
base = rmfield(results{1}.c,{'rtol','atol'});
q = struct();
fields = {'G','h','P','u','Q','v'};
for k = 1:numel(fields), q.(fields{k}) = p.(fields{k}); end
for field = ["dG","dh","dP","du","dQ","dv"], q.(field) = @(t) forbidden(); end
for k = 1:8
    r = results{k}; m = r.metrics;
    fprintf('独立验收 %d/8 %s\n',k,r.model);
    stored = load(fullfile(source,sprintf('run_%02d.mat',k)));
    assert(isequal(stored.hash,data.hash));
    assert(isequaln(rmfield(stored.r,'sol'),rmfield(r,'sol')));
    assert(isequaln(stored.r.sol.x,r.sol.x) && isequaln(stored.r.sol.y,r.sol.y));
    assert(r.model == models(mod(k-1,4)+1) && stats.model(k) == r.model);
    assert(r.status == "complete" && r.last >= c.finish-1e-10);
    assert(isequal(r.t,c.t) && isequal(size(r.g),[numel(c.t),7]));
    assert(all(isfinite(r.g),'all') && all(isfinite(r.raw.g),'all'));
    assert(all(diff(r.raw.t) > 0) && r.raw.t(1) == c.start);
    assert(norm(r.raw.g(1,:).'-c.g0) == 0 && norm(r.g(1,:).'-c.g0) == 0);
    assert(isequaln(r.c.g0,c.g0) && isequaln(rmfield(r.c,{'rtol','atol'}),base));
    if k <= 4
        assert(r.c.rtol == c.rtol && r.c.atol == c.atol && stats.group(k) == "main");
    else
        assert(r.c.rtol == 1e-9 && r.c.atol == 1e-11 && stats.group(k) == "precision");
    end
    reconstructed = deval(r.sol,r.t).';
    interpolation = max(vecnorm(reconstructed-r.g,2,2));
    assert(interpolation < 1e-12);
    xd = vecnorm(r.g(:,1:2)-ref.xd,2,2);
    x = vecnorm(r.g(:,1:2)-ref.x,2,2);
    g = vecnorm(r.g-ref.g,2,2);
    prediction = vecnorm(m.dp-ref.velocity,2,2);
    count = numel(r.t);
    residual = zeros(count,1); equality = residual; inequality = residual;
    for j = 1:count
        time = r.t(j); state = r.g(j,:).';
        G = p.G(time); G = (G+G.')/2; P = p.P(time); Q = p.Q(time);
        slack = p.v(time)-Q*state(1:2); mu = state(4:7);
        sigma = hypot(hypot(slack,mu),sqrt(c.delta));
        joint = slack+mu; phi = joint-sigma; take = joint >= 0;
        phi(take) = (2*slack(take).*mu(take)-c.delta)./(joint(take)+sigma(take));
        xi = [G*state(1:2)+p.h(time)+P.'*state(3)+Q.'*mu;P*state(1:2)-p.u(time);phi];
        residual(j) = norm(xi); equality(j) = abs(xi(3));
        inequality(j) = max([0;-slack]);
    end
    metric = max([max(abs(xd-m.xd)),max(abs(x-m.x)),max(abs(g-m.g)), ...
        max(abs(prediction-m.predictor)),max(abs(residual-m.residual)), ...
        max(abs(equality-m.equality)),max(abs(inequality-m.inequality))]);
    assert(metric < 1e-10);
    take = r.t >= c.tail;
    values = [weight(r.t,xd,c.tail),weight(r.t,x,c.tail), ...
        weight(r.t,g,c.tail),weight(r.t,prediction,c.tail)];
    expected = [stats.xd(k),stats.x(k),stats.g(k),stats.prediction(k)];
    weighted = max(abs(values-expected));
    assert(weighted <= 1e-12*max(1,max(expected)));
    assert(abs(max(residual(take))-stats.residual(k)) < 1e-10);
    assert(abs(max(equality(take))-stats.equality(k)) < 1e-10);
    assert(abs(max(inequality(take))-stats.inequality(k)) < 1e-10);
    item = struct('group',stats.group(k),'model',r.model,'complete',r.complete, ...
        'hash',true,'same',true,'interpolation',interpolation,'metric',metric, ...
        'weighted',weighted,'xd',values(1),'x',values(2),'g',values(3), ...
        'prediction',values(4),'residual',max(residual(take)));
    if isempty(rows), rows = item; else, rows(end+1) = item; end
    sparse = coarse >= c.tail;
    item = struct('group',stats.group(k),'model',r.model, ...
        'dense_xd',values(1),'sparse_xd',weight(coarse,xd(ids),c.tail), ...
        'dense_x',values(2),'sparse_x',weight(coarse,x(ids),c.tail), ...
        'dense_g',values(3),'sparse_g',weight(coarse,g(ids),c.tail), ...
        'dense_gpeak',max(g(take)),'sparse_gpeak',max(g(ids(sparse))), ...
        'dense_residual',max(residual(take)),'sparse_residual',max(residual(ids(sparse))), ...
        'dense_xdpeak',max(xd(take)),'sparse_xdpeak',max(xd(ids(sparse))));
    if isempty(sampling), sampling = item; else, sampling(end+1) = item; end
    if r.model == "M7"
        physical = [vecnorm(m.dg(:,1:2),2,2),abs(m.dg(:,3)),vecnorm(m.dg(:,4:7),2,2)];
        maximum = max(physical);
        assert(all(maximum <= c.limit.*(1+1e-10)));
        alpha = 0; common = 0; metadata = 0;
        for j = 1:count
            [dg,a] = pp_system(r.t(j),r.g(j,:).',q,r.c,r.model);
            raw = a.dp-r.c.gamma*a.dc/hypot(norm(a.dc./a.scale),r.c.eps);
            block = [norm(raw(1:2)),abs(raw(3)),norm(raw(4:7))];
            expected = min([1,r.c.limit./block]);
            alpha = max(alpha,abs(a.alpha-expected));
            common = max(common,norm(dg-expected*raw)/max(1,norm(dg)));
            metadata = max(metadata,max([abs(a.alpha-m.alpha(j)),norm(dg-m.dg(j,:).')]));
        end
        assert(alpha < 1e-12 && common < 1e-10 && metadata < 1e-10);
        flag = double(m.alpha < 1-1e-12);
        duration = sum(diff(r.t).*(flag(1:end-1)+flag(2:end))/2);
        at = r.t(take); tailflag = flag(take);
        tailduration = sum(diff(at).*(tailflag(1:end-1)+tailflag(2:end))/2);
        limited = sum(flag);
        assert(limited == stats.limited(k));
        margin = -Inf;
        for block = 1:3
            blocks = {1:2,3,4:7};
            columns = blocks{block};
            for grid = 1:2
                if grid == 1, times = r.raw.t; state = r.raw.g(:,columns);
                else, times = r.t; state = r.g(:,columns); end
                normg = vecnorm(state,2,2);
                tolerance = 100*(r.c.atol+r.c.rtol*max([ones(numel(times),1),normg],[],2));
                excess = vecnorm(diff(state),2,2)-r.c.limit(block)*diff(times) ...
                    -max(tolerance(1:end-1),tolerance(2:end));
                margin = max(margin,max(excess));
            end
        end
        assert(margin <= 0);
        item = struct('group',stats.group(k),'model',r.model,'x',maximum(1), ...
            'mu1',maximum(2),'mu2',maximum(3),'alpha',alpha,'common',common, ...
            'metadata',metadata,'margin',margin,'limited',limited, ...
            'duration',duration,'fraction',duration/(c.finish-c.start), ...
            'tailduration',tailduration,'tailfraction',tailduration/(c.finish-c.tail));
        if isempty(limits), limits = item; else, limits(end+1) = item; end
    end
end

function checks = refinement(folder,source,out)
% 加密评估网格，不重新积分模型。
root = fileparts(fileparts(folder));
addpath(root,'-begin'); addpath(folder,'-begin');
diary(fullfile(out,'refinement.log')); clean = onCleanup(@() diary('off'));
data = load(fullfile(source,'all.mat'));
c = data.c; p = example2_problem();
t = unique([c.t;(c.start:0.00025:c.finish).']);
ref = pp_reference(t,p,c.delta);
coarse = (c.start:0.005:c.finish).';
sparse = pp_reference(coarse,p,c.delta);
take = coarse >= c.tail;
checks.points = numel(t);
checks.bias = weight(t,ref.bias,c.tail);
checks.oldbias = sqrt(mean(sparse.bias(take).^2));
checks.coarsebias = weight(coarse,sparse.bias,c.tail);
checks.mainbias = weight(c.t,data.ref.bias,c.tail);
rows = repmat(struct(),0,1);
for k = 1:numel(data.results)
    r = data.results{k}; g = deval(r.sol,t).';
    xd = vecnorm(g(:,1:2)-ref.xd,2,2);
    x = vecnorm(g(:,1:2)-ref.x,2,2);
    eg = vecnorm(g-ref.g,2,2);
    values = [weight(t,xd,c.tail),weight(t,x,c.tail),weight(t,eg,c.tail)];
    old = [data.stats.xd(k),data.stats.x(k),data.stats.g(k)];
    difference = abs(values-old)./max(realmin,old);
    gs = deval(r.sol,coarse).';
    xs = vecnorm(gs(:,1:2)-sparse.x,2,2);
    item = struct('group',data.stats.group(k),'model',r.model, ...
        'xd',values(1),'x',values(2),'g',values(3), ...
        'xdrelative',difference(1),'xrelative',difference(2),'grelative',difference(3), ...
        'oldmeanx',sqrt(mean(xs(take).^2)),'coarsex',weight(coarse,xs,c.tail));
    if isempty(rows), rows = item; else, rows(end+1) = item; end
end
refinement = struct2table(rows);
writetable(refinement,fullfile(out,'refinement.csv'));
bias = table(checks.bias,checks.mainbias,checks.coarsebias,checks.oldbias, ...
    'VariableNames',{'refined','main','coarse','oldmean'});
writetable(bias,fullfile(out,'bias.csv'));
extra = table();
file = fullfile(folder,'extra','all.mat');
if isfile(file)
    ext = load(file);
    assert(isequal(ext.hash,data.hash));
    rows = repmat(struct(),0,1);
    for k = 1:numel(ext.results)
        r = ext.results{k}; assert(r.complete);
        assert(isequal(r.t,c.t) && isequal(r.c.g0,c.g0));
        expected = c; expected.rtol = 1e-9; expected.atol = 1e-11;
        if ext.stats.group(k) == "kappa01"
            assert(r.model == "M6" && r.c.kappa == 0.01);
            expected.kappa = 0.01;
        else
            assert(r.model == "M7" && isequal(r.c.scale,[1,1,1]));
            expected.scale = [1,1,1];
        end
        assert(isequaln(r.c,expected),'额外组包含预期之外的配置差异。');
        assert(r.c.rtol == 1e-9 && r.c.atol == 1e-11);
        g = deval(r.sol,t).';
        xd = vecnorm(g(:,1:2)-ref.xd,2,2);
        x = vecnorm(g(:,1:2)-ref.x,2,2);
        eg = vecnorm(g-ref.g,2,2);
        values = [weight(t,xd,c.tail),weight(t,x,c.tail),weight(t,eg,c.tail)];
        original = [weight(r.t,vecnorm(r.g(:,1:2)-data.ref.xd,2,2),c.tail), ...
            weight(r.t,vecnorm(r.g(:,1:2)-data.ref.x,2,2),c.tail), ...
            weight(r.t,vecnorm(r.g-data.ref.g,2,2),c.tail)];
        assert(max(abs(original-[ext.stats.xd(k),ext.stats.x(k),ext.stats.g(k)])) < 1e-12);
        take = t >= c.tail;
        flag = double(r.metrics.alpha < 1-1e-12);
        duration = sum(diff(r.t).*(flag(1:end-1)+flag(2:end))/2);
        physical = [vecnorm(r.metrics.dg(:,1:2),2,2), ...
            abs(r.metrics.dg(:,3)),vecnorm(r.metrics.dg(:,4:7),2,2)];
        if r.model == "M7", assert(all(max(physical) <= r.c.limit.*(1+1e-10))); end
        item = struct('group',ext.stats.group(k),'model',r.model, ...
            'xd',values(1),'x',values(2),'g',values(3),'gpeak',max(eg(take)), ...
            'xdrelative',abs(values(1)-original(1))/original(1), ...
            'xrelative',abs(values(2)-original(2))/original(2), ...
            'grelative',abs(values(3)-original(3))/original(3), ...
            'limited',sum(flag),'duration',duration,'fraction',duration/(c.finish-c.start));
        if isempty(rows), rows = item; else, rows(end+1) = item; end
    end
    extra = struct2table(rows); writetable(extra,fullfile(out,'extra.csv'));
    checks.extra = digest(fullfile(folder,'pp_extra.m'));
end
checks.passed = true;
save(fullfile(out,'refinement.mat'),'checks','refinement','bias','extra');
disp(checks); disp(refinement); disp(bias); disp(extra);
end
integrity = struct2table(rows);
sampling = struct2table(sampling);
limits = struct2table(limits);
rows = repmat(struct(),0,1);
for k = 1:4
    a = results{k}; b = results{k+4};
    x = vecnorm(a.g(:,1:2)-b.g(:,1:2),2,2);
    g = vecnorm(a.g-b.g,2,2);
    item = struct('model',a.model,'xpeak',max(x),'gpeak',max(g), ...
        'xtail',weight(c.t,x,c.tail),'gtail',weight(c.t,g,c.tail), ...
        'main_xd',integrity.xd(k),'strict_xd',integrity.xd(k+4), ...
        'xdabsolute',abs(integrity.xd(k)-integrity.xd(k+4)), ...
        'xdrelative',abs(integrity.xd(k)-integrity.xd(k+4))/max(realmin,integrity.xd(k+4)), ...
        'main_residual',integrity.residual(k),'strict_residual',integrity.residual(k+4));
    if isempty(rows), rows = item; else, rows(end+1) = item; end
end
precision = struct2table(rows);
checks.complete = true;
checks.same = true;
checks.passed = true;
writetable(integrity,fullfile(out,'integrity.csv'));
writetable(sampling,fullfile(out,'sampling.csv'));
writetable(limits,fullfile(out,'limits.csv'));
writetable(precision,fullfile(out,'precision.csv'));
save(fullfile(out,'results.mat'),'checks','integrity','sampling','limits','precision');
disp(checks); disp(integrity); disp(sampling); disp(limits); disp(precision);
fprintf('验收通过；限速 duration 是按时间积分的采样估计，limited 仍是点数。\n');
end

function value = weight(t,v,start)
% 按相邻时间段独立计算均方误差。
take = t >= start; t = t(take); v = v(take);
value = sqrt(sum(diff(t).*(v(1:end-1).^2+v(2:end).^2)/2)/(t(end)-t(1)));
end

function code = digest(path)
% 读取源码的 SHA256 校验值。
file = fopen(path,'rb'); assert(file >= 0);
clean = onCleanup(@() fclose(file)); bytes = fread(file,Inf,'*uint8');
md = java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
code = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function value = forbidden()
% 禁止模型读取解析系数导数。
error('Verify:Derivative','模型读取了禁止使用的解析系数导数。');
value = [];
end
