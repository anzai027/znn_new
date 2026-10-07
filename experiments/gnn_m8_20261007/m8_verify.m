function checks = m8_verify(name)
% 独立复核保存结果，不重新运行积分器。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
old = fullfile(root,'experiments','gnn_predictor_20261006');
addpath(old,'-begin'); addpath(root,'-begin'); addpath(folder,'-begin');
if nargin == 0, name = 'data'; end
out = fullfile(folder,'verification');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'run.log')); clean = onCleanup(@() diary('off'));
source = fullfile(folder,name);
data = load(fullfile(source,'all.mat'));
c = data.c; p = example2_problem();
checkhash(data.hash,folder,root,old);
assert(numel(data.results) == 9 && height(data.stats) == 9);
labels = ["M8","M8a","M8","M8a","M8a fixed", ...
    "M8 lp1e-9","M8 lp1e-8","M8 h1e-6","M8a h1e-6"];
models = ["M8","M8A","M8","M8A","M8A","M8","M8","M8","M8A"];
base = pp_config();
for key = ["g0","t","start","finish","tail","delta","step","gamma","eps","h"]
    assert(isequaln(c.(key),base.(key)),'公共配置发生变化：%s',key);
end
assert(c.rtol == 1e-7 && c.atol == 1e-9);
ref = pp_reference(c.t,p,c.delta);
checks.reference = max(vecnorm(ref.g-data.ref.g,2,2));
checks.residual = max(ref.error);
checks.velocity = max(vecnorm(ref.velocity-data.ref.velocity,2,2));
assert(checks.reference < 1e-10 && checks.residual < 1e-9 && checks.velocity < 1e-8);
assert(isequal(ref.t,c.t));
csv = readtable(fullfile(source,'metrics.csv'),'TextType','string');
comparecsv(csv,data.stats);
rows = table(); sampling = table(); precision = table();
coarse = (c.start:0.005:c.finish).';
ids = round(interp1(c.t,(1:numel(c.t)).',coarse,'nearest'));
assert(max(abs(c.t(ids)-coarse)) < 1e-11);
for k = 1:numel(data.results)
    r = data.results{k};
    fprintf('验收新模型 %d/%d %s\n',k,numel(data.results),r.model);
    stored = load(fullfile(source,sprintf('run_%02d.mat',k)));
    assert(isequal(stored.hash,data.hash));
    assert(isequaln(rmfield(stored.r,'sol'),rmfield(r,'sol')),'逐组 MAT 与汇总 MAT 不一致。');
    assert(isequaln(stored.r.sol.x,r.sol.x) && isequaln(stored.r.sol.y,r.sol.y));
    expected = pn_config(models(k)); group = "main";
    if k > 2, expected.rtol = 1e-9; expected.atol = 1e-11; group = "precision"; end
    if k > 4, group = "ablation"; end
    if k == 5, expected.correct = "fixed"; end
    if k == 6, expected.predict = 1e-9; end
    if k == 7, expected.predict = 1e-8; end
    if k >= 8, expected.h = 1e-6; end
    assert(isequaln(r.c,expected),'该组含预期之外的配置变化。');
    assert(r.label == labels(k) && r.model == models(k) && r.group == group);
    assert(data.stats.label(k) == labels(k) && data.stats.group(k) == group);
    assert(r.complete && r.status == "complete" && r.last >= c.finish-1e-10);
    assert(isequal(r.t,c.t) && isequal(size(r.g),[numel(c.t),7]));
    assert(all(isfinite(r.g),'all') && all(isfinite(r.raw.g),'all'));
    assert(all(diff(r.raw.t) > 0) && r.raw.t(1) == c.start);
    assert(isequal(r.g(1,:).',c.g0) && isequal(r.raw.g(1,:).',c.g0));
    for key = ["g0","t","start","finish","tail","delta","step","gamma","eps"]
        assert(isequaln(r.c.(key),c.(key)),'单组配置发生变化：%s',key);
    end
    assert((r.c.rtol == 1e-7 && r.c.atol == 1e-9) || ...
        (r.c.rtol == 1e-9 && r.c.atol == 1e-11));
    dg = deval(r.sol,c.t).';
    interpolation = max(vecnorm(dg-r.g,2,2));
    assert(interpolation < 1e-12);
    a = measure(r.t,r.g,ref,p,r.c);
    difference = max([max(abs(a.xd-r.metrics.xd)),max(abs(a.x-r.metrics.x)), ...
        max(abs(a.g-r.metrics.g)),max(abs(a.residual-r.metrics.residual)), ...
        max(abs(a.equality-r.metrics.equality)),max(abs(a.inequality-r.metrics.inequality))]);
    assert(difference < 1e-10);
    values = numbers(r.t,a,r.c.tail);
    expected = [data.stats.xd(k),data.stats.x(k),data.stats.g(k), ...
        data.stats.gpeak(k),data.stats.residual(k)];
    error = max(abs(values-expected));
    assert(error < 1e-12*max(1,max(abs(expected))));
    predictor = vecnorm(r.metrics.dp-ref.velocity,2,2);
    assert(max(abs(predictor-r.metrics.predictor)) < 1e-10);
    if ismember('prediction',data.stats.Properties.VariableNames)
        assert(abs(weight(r.t,predictor,r.c.tail)-data.stats.prediction(k)) < 1e-10);
    end
    row = table(string(r.label),r.c.rtol,r.c.atol,r.calls,interpolation,difference,error, ...
        values(1),values(2),values(3),values(4),values(5), ...
        'VariableNames',{'model','rtol','atol','calls','interpolation','metric','weighted', ...
        'xd','x','g','gpeak','residual'});
    rows = [rows;row];
    b = struct();
    for key = ["xd","x","g","residual","equality","inequality"]
        b.(key) = a.(key)(ids);
    end
    sparse = numbers(coarse,b,r.c.tail);
    row = table(string(r.label),r.c.rtol,values(1),sparse(1),values(3),sparse(3), ...
        values(4),sparse(4),values(5),sparse(5), ...
        'VariableNames',{'model','rtol','densexd','coarsexd','denseg','coarseg', ...
        'densegpeak','coarsegpeak','denseresidual','coarseresidual'});
    sampling = [sampling;row];
end
writetable(rows,fullfile(out,'integrity.csv'));
writetable(sampling,fullfile(out,'sampling.csv'));
for k = 1:numel(data.results)
    a = data.results{k};
    if a.c.rtol ~= 1e-7, continue, end
    config = rmfield(a.c,{'rtol','atol'});
    for j = 1:numel(data.results)
        b = data.results{j};
        if b.c.rtol ~= 1e-9 || a.model ~= b.model || ...
                ~isequaln(config,rmfield(b.c,{'rtol','atol'})), continue, end
        x = vecnorm(a.g(:,1:2)-b.g(:,1:2),2,2);
        g = vecnorm(a.g-b.g,2,2);
        row = table(string(a.label),weight(c.t,x,c.tail),weight(c.t,g,c.tail), ...
            abs(rows.xd(k)-rows.xd(j)),abs(rows.xd(k)-rows.xd(j))/max(realmin,rows.xd(j)), ...
            abs(rows.gpeak(k)-rows.gpeak(j))/max(realmin,rows.gpeak(j)), ...
            'VariableNames',{'model','xdifference','gdifference','xdabsolute','xdrelative','peakrelative'});
        precision = [precision;row];
    end
end
writetable(precision,fullfile(out,'precision.csv'));
t = unique([c.t;(c.start:0.00025:c.finish).']);
for peak = c.peak
    t = unique([t;(peak-0.001:5e-7:peak+0.001).']);
end
fine = pp_reference(t,p,c.delta);
assert(max(fine.error) < 1e-9 && max(fine.derivative) < 1e-9);
refinement = table();
for k = 1:numel(data.results)
    r = data.results{k};
    fprintf('加密评估 %d/%d %s\n',k,numel(data.results),r.model);
    g = deval(r.sol,t).';
    a = measure(t,g,fine,p,r.c);
    values = numbers(t,a,c.tail);
    oldvalues = [rows.xd(k),rows.x(k),rows.g(k),rows.gpeak(k),rows.residual(k)];
    relative = abs(values-oldvalues)./max(oldvalues,realmin);
    row = table(string(r.label),r.c.rtol,values(1),values(2),values(3),values(4),values(5), ...
        relative(1),relative(2),relative(3),relative(4),relative(5), ...
        'VariableNames',{'model','rtol','xd','x','g','gpeak','residual', ...
        'xdrelative','xrelative','grelative','peakrelative','residualrelative'});
    refinement = [refinement;row];
end
writetable(refinement,fullfile(out,'refinement.csv'));
baseline = table();
paths = {fullfile(old,'data','run_05.mat'),fullfile(old,'data','run_06.mat'), ...
    fullfile(old,'data','run_07.mat'),fullfile(old,'extra','run_01.mat'), ...
    fullfile(old,'data','run_08.mat'),fullfile(old,'vfcr','run_02.mat')};
labels = ["M3","M5","M6 k0.05","M6 k0.01","M7","VFCR-S"];
for k = 1:numel(paths)
    item = load(paths{k}); r = item.r;
    fprintf('验收旧基准 %d/%d %s\n',k,numel(paths),labels(k));
    checkhash(item.hash,old,root,old);
    assert(r.complete && r.status == "complete");
    for key = ["g0","t","start","finish","tail","delta","step","gamma","eps"]
        assert(isequaln(r.c.(key),c.(key)),'旧基准配置不同：%s',key);
    end
    assert(r.c.rtol == 1e-9 && r.c.atol == 1e-11);
    assert(isequal(r.t,c.t) && isequal(r.g(1,:).',c.g0));
    if k == 6
        assert(size(r.y,2) == 14 && isequal(r.g,r.y(:,1:7)));
        assert(all(r.y(1,8:14) == 0));
    end
    a = measure(r.t,r.g,ref,p,r.c);
    values = numbers(r.t,a,c.tail);
    assert(max([max(abs(a.xd-r.metrics.xd)),max(abs(a.x-r.metrics.x)), ...
        max(abs(a.g-r.metrics.g)),max(abs(a.residual-r.metrics.residual))]) < 1e-10);
    g = deval(r.sol,t).'; g = g(:,1:7);
    a = measure(t,g,fine,p,r.c);
    refined = numbers(t,a,c.tail);
    change = abs(values-refined)./max(values,realmin);
    row = table(labels(k),r.c.rtol,r.c.atol,r.calls,r.seconds,values(1),values(2), ...
        values(3),values(4),values(5),refined(1),refined(2),refined(3),refined(4),refined(5), ...
        max(change(1:3)),max(change(4:5)), ...
        'VariableNames',{'model','rtol','atol','calls','seconds','xd','x','g','gpeak','residual', ...
        'refinedxd','refinedx','refinedg','refinedgpeak','refinedresidual','rmsrelative','peakrelative'});
    baseline = [baseline;row];
end
raw = load(fullfile(old,'vfcr','run_03.mat'));
assert(~raw.r.complete && raw.r.status == "budget" && raw.r.last < c.tail);
checkhash(raw.hash,old,root,old);
writetable(baseline,fullfile(out,'baseline.csv'));
comparison = readtable(fullfile(source,'comparison.csv'),'TextType','string');
stored = load(fullfile(source,'comparison.mat'));
comparecsv(comparison,stored.comparison);
assert(height(comparison) == 13);
for k = 1:height(comparison)
    if k <= 6
        assert(comparison.origin(k) == "old" && comparison.group(k) == "precision");
        assert(comparison.label(k) == baseline.model(k));
        expected = [baseline.xd(k),baseline.x(k),baseline.g(k),baseline.gpeak(k),baseline.residual(k)];
    else
        j = k-4;
        assert(comparison.origin(k) == "new" && comparison.label(k) == data.results{j}.label);
        expected = [rows.xd(j),rows.x(j),rows.g(j),rows.gpeak(j),rows.residual(j)];
    end
    assert(isfile(fullfile(root,comparison.file(k))));
    assert(comparison.rtol(k) == 1e-9 && comparison.atol(k) == 1e-11 && comparison.delta(k) == c.delta);
    actual = [comparison.xd(k),comparison.x(k),comparison.g(k),comparison.gpeak(k),comparison.residual(k)];
    assert(max(abs(actual-expected)) < 1e-12);
end
checks.points = numel(c.t); checks.refined = numel(t);
checks.bias = weight(c.t,ref.bias,c.tail);
checks.refinedbias = weight(t,fine.bias,c.tail);
checks.rawlast = raw.r.last; checks.rawcalls = raw.r.calls;
checks.hash = true; checks.same = true; checks.passed = true;
save(fullfile(out,'results.mat'),'checks','rows','precision','sampling','refinement','baseline','comparison','-v7.3');
disp(checks); disp(rows); disp(precision); disp(refinement); disp(baseline);
end

function a = measure(t,g,ref,p,c)
% 使用独立残差表达式复算误差。
a.xd = vecnorm(g(:,1:2)-ref.xd,2,2);
a.x = vecnorm(g(:,1:2)-ref.x,2,2);
a.g = vecnorm(g-ref.g,2,2);
a.residual = zeros(numel(t),1); a.equality = a.residual; a.inequality = a.residual;
for j = 1:numel(t)
    time = t(j); state = g(j,:).';
    G = p.G(time); G = (G+G.')/2; P = p.P(time); Q = p.Q(time);
    slack = p.v(time)-Q*state(1:2); mu = state(4:7);
    sigma = hypot(hypot(slack,mu),sqrt(c.delta));
    joint = slack+mu; phi = joint-sigma; take = joint >= 0;
    phi(take) = (2*slack(take).*mu(take)-c.delta)./(joint(take)+sigma(take));
    xi = [G*state(1:2)+p.h(time)+P.'*state(3)+Q.'*mu;P*state(1:2)-p.u(time);phi];
    a.residual(j) = norm(xi); a.equality(j) = abs(xi(3));
    a.inequality(j) = max([0;-slack]);
end
end

function values = numbers(t,a,start)
% 分别返回时间加权误差和密集峰值。
take = t >= start;
values = [weight(t,a.xd,start),weight(t,a.x,start),weight(t,a.g,start), ...
    max(a.g(take)),max(a.residual(take))];
end

function value = weight(t,e,start)
% 按相邻时间段独立计算均方误差。
take = t >= start; t = t(take); e = e(take);
value = sqrt(sum(diff(t).*(e(1:end-1).^2+e(2:end).^2)/2)/(t(end)-t(1)));
end

function comparecsv(csv,stats)
% 核对 CSV 与 MAT 中的每一个表格值。
assert(height(csv) == height(stats));
assert(isequal(csv.Properties.VariableNames,stats.Properties.VariableNames));
for key = string(stats.Properties.VariableNames)
    a = csv.(key); b = stats.(key);
    if isnumeric(b) || islogical(b)
        assert(isequal(isnan(a),isnan(b)) && isequal(isinf(a),isinf(b)));
        take = isfinite(a) & isfinite(b);
        assert(all(abs(double(a(take))-double(b(take))) <= 1e-12*max(1,abs(double(b(take))))));
    else
        assert(isequal(string(a),string(b)));
    end
end
end

function checkhash(hash,folder,root,old)
% 按运行记录核对新旧源码的 SHA256。
if istable(hash)
    assert(isequal(hash.Properties.VariableNames,{'file','code'}));
    for k = 1:height(hash)
        path = fullfile(root,hash.file(k));
        assert(isfile(path) && strcmp(digest(path),hash.code(k)), ...
            '源码与运行哈希不同：%s',hash.file(k));
    end
    return
end
keys = fieldnames(hash);
for k = 1:numel(keys)
    key = keys{k};
    if isstruct(hash.(key)), checkhash(hash.(key),folder,root,old); continue, end
    if strcmp(key,'problem'), file = fullfile(root,'example2_problem.m');
    elseif strcmp(key,'my_system'), file = fullfile(root,'my_system.m');
    elseif strcmp(key,'fsmooth'), file = fullfile(root,'experiments','gnn_v2_20261002','fsmooth.m');
    elseif strcmp(key,'driver'), file = fullfile(old,'pp_vfcr.m');
    elseif startsWith(key,'pp_') || any(strcmp(key,{'run_predictor','verify_results','audit','reference_check'}))
        file = fullfile(old,[key,'.m']);
    else, file = fullfile(folder,[key,'.m']); end
    assert(isfile(file),'哈希字段未能对应文件：%s',key);
    assert(strcmp(digest(file),hash.(key)),'源码与运行哈希不同：%s',key);
end
end

function code = digest(path)
% 读取文件的 SHA256 校验值。
file = fopen(path,'rb'); assert(file >= 0); clean = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8'); md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes); code = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
