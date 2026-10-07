function checks = m8_pverify(name)
% 独立验收更严格容差的数据，不重新积分。
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
old = fullfile(root,'experiments','gnn_predictor_20261006');
addpath(old,'-begin'); addpath(root,'-begin'); addpath(folder,'-begin');
if nargin == 0, name = 'precision'; end
source = fullfile(folder,name); out = fullfile(folder,'verification');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'precision.log')); guard = onCleanup(@() diary('off'));
data = load(fullfile(source,'all.mat'));
assert(numel(data.results) == 4 && height(data.stats) == 4);
assert(istable(data.hash));
for k = 1:height(data.hash)
    assert(strcmp(digest(fullfile(root,data.hash.file(k))),data.hash.code(k)), ...
        '严格复核源码发生变化：%s',data.hash.file(k));
end
csv = readtable(fullfile(source,'metrics.csv'),'TextType','string');
comparecsv(csv,data.stats);
c = data.c; p = example2_problem();
shared = load(fullfile(folder,'data','reference.mat'),'ref','c');
for key = ["g0","t","start","finish","tail","delta","step","gamma","eps"]
    assert(isequaln(c.(key),shared.c.(key)),'严格复核公共配置不同：%s',key);
end
assert(isequaln(data.ref,shared.ref));
assert(c.rtol == 1e-11 && c.atol == 1e-13);
v = struct('r1',1,'r2',1,'a',1,'p',0.5,'q',0.5,'lambda1',1,'lambda2',1);
assert(isequaln(data.v,v));
t = unique([c.t;(c.start:0.00025:c.finish).']);
for peak = c.peak, t = unique([t;(peak-0.001:5e-7:peak+0.001).']); end
fine = pp_reference(t,p,c.delta);
assert(max(fine.error) < 1e-9 && max(fine.derivative) < 1e-9);
models = ["M8","M8A","M8A","VFCR-S"];
labels = ["M8 h1e-5","M8a h1e-5","M8a h1e-6","VFCR-S"];
paths = {fullfile(folder,'data','run_03.mat'),fullfile(folder,'data','run_04.mat'), ...
    fullfile(folder,'data','run_09.mat'),fullfile(old,'vfcr','run_02.mat')};
rows = table(); sensitivity = table();
for k = 1:4
    r = data.results{k};
    fprintf('严格数据验收 %d/4 %s\n',k,labels(k));
    stored = load(fullfile(source,sprintf('run_%02d.mat',k)));
    assert(isequaln(stored.hash,data.hash) && isequaln(stored.v,v));
    assert(isequaln(rmfield(stored.r,'sol'),rmfield(r,'sol')));
    assert(isequaln(stored.r.sol.x,r.sol.x) && isequaln(stored.r.sol.y,r.sol.y));
    expected = c;
    if k == 2 || k == 3, expected.correct = "adaptive"; end
    if k == 3, expected.h = 1e-6; end
    assert(isequaln(r.c,expected));
    assert(r.label == labels(k) && r.model == models(k) && r.group == "refinement");
    assert(r.complete && r.status == "complete" && r.last >= c.finish-1e-10);
    assert(isequal(r.t,c.t) && isequal(size(r.g),[numel(c.t),7]));
    assert(all(isfinite(r.g),'all') && all(diff(r.raw.t) > 0));
    assert(isequal(r.g(1,:).',c.g0) && isequal(r.raw.g(1,:).',c.g0));
    sampled = deval(r.sol,r.t).'; sampled = sampled(:,1:7);
    interpolation = max(vecnorm(sampled-r.g,2,2));
    assert(interpolation < 1e-12);
    if k == 4
        assert(size(r.y,2) == 14 && isequal(r.g,r.y(:,1:7)));
        assert(all(r.y(1,8:14) == 0));
        assert(isnan(data.stats.h(k)) && data.stats.correct(k) == "none");
    else
        assert(data.stats.h(k) == r.c.h && data.stats.correct(k) == r.c.correct);
    end
    a = measure(r.t,r.g,data.ref,p,c.delta);
    metric = max([max(abs(a.xd-r.metrics.xd)),max(abs(a.x-r.metrics.x)), ...
        max(abs(a.g-r.metrics.g)),max(abs(a.residual-r.metrics.residual)), ...
        max(abs(a.equality-r.metrics.equality)),max(abs(a.inequality-r.metrics.inequality))]);
    assert(metric < 1e-10);
    values = numbers(r.t,a,c.tail);
    expected = [data.stats.xd(k),data.stats.x(k),data.stats.g(k), ...
        data.stats.gpeak(k),data.stats.residual(k)];
    weighted = max(abs(values-expected));
    assert(weighted < 1e-12);
    dense = deval(r.sol,t).'; dense = dense(:,1:7);
    b = measure(t,dense,fine,p,c.delta);
    refined = numbers(t,b,c.tail);
    change = abs(refined-values)./max(values,realmin);
    row = table(labels(k),r.calls,interpolation,metric,weighted,values(1),values(2), ...
        values(3),values(4),values(5),refined(1),refined(2),refined(3),refined(4),refined(5), ...
        max(change(1:3)),max(change(4:5)), ...
        'VariableNames',{'model','calls','interpolation','metric','weighted','xd','x','g', ...
        'gpeak','residual','refinedxd','refinedx','refinedg','refinedgpeak', ...
        'refinedresidual','rmsrelative','peakrelative'});
    rows = [rows;row];
    prior = load(paths{k}); lower = prior.r;
    assert(lower.complete && lower.c.rtol == 1e-9 && lower.c.atol == 1e-11);
    for key = ["g0","t","start","finish","tail","delta","step","gamma","eps"]
        assert(isequaln(lower.c.(key),r.c.(key)),'较松容差的对照配置不同：%s',key);
    end
    assert(isequal(lower.t,r.t));
    if k == 4, assert(isequaln(prior.v,v)); end
    if k < 4, assert(isequaln(rmfield(lower.c,{'rtol','atol'}),rmfield(r.c,{'rtol','atol'}))); end
    d = measure(lower.t,lower.g,data.ref,p,c.delta);
    original = numbers(lower.t,d,c.tail);
    xdiff = vecnorm(lower.g(:,1:2)-r.g(:,1:2),2,2);
    gdiff = vecnorm(lower.g-r.g,2,2);
    relative = abs(values-original)./max(values,realmin);
    row = table(labels(k),original(1),values(1),relative(1),original(3),values(3),relative(3), ...
        original(4),values(4),relative(4),original(5),values(5),relative(5), ...
        weight(r.t,xdiff,c.tail),weight(r.t,gdiff,c.tail), ...
        'VariableNames',{'model','strictxd','tighterxd','xdrelative','strictg','tighterg', ...
        'grelative','strictgpeak','tightergpeak','peakrelative','strictresidual', ...
        'tighterresidual','residualrelative','xdifference','gdifference'});
    sensitivity = [sensitivity;row];
end
checks.passed = true; checks.hash = true; checks.same = true;
checks.points = numel(c.t); checks.refined = numel(t);
checks.reference = max(fine.error); checks.derivative = max(fine.derivative);
writetable(rows,fullfile(out,'tighter.csv'));
writetable(sensitivity,fullfile(out,'sensitivity.csv'));
save(fullfile(out,'tighter.mat'),'checks','rows','sensitivity','-v7.3');
disp(checks); disp(rows); disp(sensitivity);
end

function a = measure(t,g,ref,p,delta)
% 使用独立残差表达式复算误差。
a.xd = vecnorm(g(:,1:2)-ref.xd,2,2); a.x = vecnorm(g(:,1:2)-ref.x,2,2);
a.g = vecnorm(g-ref.g,2,2);
a.residual = zeros(numel(t),1); a.equality = a.residual; a.inequality = a.residual;
for j = 1:numel(t)
    time = t(j); state = g(j,:).';
    G = p.G(time); G = (G+G.')/2; P = p.P(time); Q = p.Q(time);
    slack = p.v(time)-Q*state(1:2); mu = state(4:7);
    sigma = hypot(hypot(slack,mu),sqrt(delta)); joint = slack+mu;
    phi = joint-sigma; take = joint >= 0;
    phi(take) = (2*slack(take).*mu(take)-delta)./(joint(take)+sigma(take));
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
        a = string(a); b = string(b); a(ismissing(a)) = ""; b(ismissing(b)) = "";
        assert(isequal(a,b));
    end
end
end

function code = digest(path)
% 读取文件的 SHA256 校验值。
file = fopen(path,'rb'); assert(file >= 0); guard = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8'); md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes); code = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
