function precision(name)
% 用更严格容差区分模型误差与积分误差。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
prior = fullfile(root,'experiments','gnn_predictor_20261006');
smooth = fullfile(root,'experiments','gnn_v2_20261002');
addpath(smooth,'-begin');
addpath(root,'-begin');
addpath(prior,'-begin');
addpath(folder,'-begin');
assert(strcmpi(which('pn_system'),fullfile(folder,'pn_system.m')));
assert(strcmpi(which('my_system'),fullfile(root,'my_system.m')));
assert(strcmpi(which('fsmooth'),fullfile(smooth,'fsmooth.m')));
if nargin == 0, name = 'precision'; end
out = fullfile(folder,name);
assert(~isfile(fullfile(out,'all.mat')),'已有复核结果，请使用新目录名。');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'run.log'));
guard = onCleanup(@() diary('off'));
files = ["pn_system.m","pn_config.m","pn_run.m","pn_metrics.m", ...
    "m8_system.m","m8a_system.m","precision.m"];
paths = "experiments/gnn_m8_20261007/"+files;
paths = [paths,"experiments/gnn_predictor_20261006/pp_parts.m", ...
    "experiments/gnn_predictor_20261006/pp_config.m", ...
    "experiments/gnn_predictor_20261006/pp_reference.m", ...
    "experiments/gnn_predictor_20261006/pp_vfcr.m", ...
    "experiments/gnn_v2_20261002/fsmooth.m","my_system.m","example2_problem.m"];
codes = strings(size(paths));
for k = 1:numel(paths), codes(k) = digest(fullfile(root,paths(k))); end
hash = table(paths(:),codes(:),'VariableNames',{'file','code'});
c = pn_config("M8");
c.rtol = 1e-11;
c.atol = 1e-13;
p = example2_problem();
source = load(fullfile(prior,'data','reference.mat'),'ref','c');
assert(isequal(c.t,source.ref.t) && isequal(c.g0,source.c.g0) && c.delta == source.c.delta);
ref = source.ref;
v = struct('r1',1,'r2',1,'a',1,'p',0.5,'q',0.5,'lambda1',1,'lambda2',1);
save(fullfile(out,'reference.mat'),'ref','c','v','hash','-v7.3');
labels = ["M8 h1e-5","M8a h1e-5","M8a h1e-6","VFCR-S"];
models = ["M8","M8A","M8A","VFCR-S"];
results = cell(4,1);
stats = table();
fprintf('MATLAB %s；严格复核容差 %.0e/%.0e；输出点数 %d\n', ...
    version,c.rtol,c.atol,numel(c.t));
for k = 1:4
    ci = c;
    if k == 2 || k == 3, ci.correct = "adaptive"; end
    if k == 3, ci.h = 1e-6; end
    fprintf('开始严格复核 %d/4 %s\n',k,labels(k));
    if k < 4
        r = pn_run(models(k),p,ci);
    else
        r = simulate(p,ci,v);
    end
    r.label = labels(k);
    r.group = "refinement";
    detail = table();
    if r.complete
        if k < 4
            [r.metrics,detail] = pn_metrics(r,ref,p);
        else
            r.metrics = measure(r,ref,p);
        end
    else
        r.metrics = struct();
    end
    row = summary(r);
    stats = [stats;row];
    results{k} = r;
    save(fullfile(out,sprintf('run_%02d.mat',k)),'r','detail','v','hash','-v7.3');
    writetable(stats,fullfile(out,'metrics.csv'));
    fprintf('完成 %s complete=%d last=%.9g calls=%d xd=%.9g gpeak=%.9g residual=%.9g\n', ...
        labels(k),r.complete,r.last,r.calls,row.xd,row.gpeak,row.residual);
end
for k = 1:height(hash)
    assert(strcmp(digest(fullfile(root,hash.file(k))),hash.code(k)), ...
        '运行过程中源码改变：%s',hash.file(k));
end
save(fullfile(out,'all.mat'),'results','stats','ref','c','v','hash','-v7.3');
disp(stats);
fprintf('严格复核已保存；本次仅检查精度，不能仅凭同一容差宣布模型排名。\n');
end

function r = simulate(p,c,v)
% 保留原 VFCR 的系数导数和积分状态。
if isfield(p,'noise'), p = rmfield(p,'noise'); end
p.phi = @(x) fsmooth(x,v,c.eps);
y0 = [c.g0;zeros(7,1)];
calls = 0;
clock = tic;
times = c.start;
states = y0.';
opt = odeset('RelTol',c.rtol,'AbsTol',c.atol,'MaxStep',c.step,'OutputFcn',@output);
r = struct('model',"VFCR-S",'c',c,'status',"running",'reason',"", ...
    'identifier',"",'warning',"",'log',"",'complete',false);
sol = [];
lastwarn('');
try
    log = evalc('sol = ode15s(@rhs,[c.start,c.finish],y0,opt);');
    r.log = string(log);
    r.sol = sol;
    r.raw.t = sol.x.';
    r.raw.y = sol.y.';
    r.raw.g = r.raw.y(:,1:7);
    r.last = sol.x(end);
    r.complete = r.last >= c.finish-1e-10 && all(isfinite(sol.y),'all');
    if r.complete
        r.t = c.t;
        r.y = deval(sol,c.t).';
        r.g = r.y(:,1:7);
        r.complete = all(isfinite(r.y),'all');
        r.status = "complete";
        if ~r.complete, r.status = "nonfinite"; end
    else
        r.status = "solver_failure";
    end
catch err
    r.status = "exception";
    r.reason = string(err.message);
    r.identifier = string(err.identifier);
    r.raw.t = times;
    r.raw.y = states;
    r.raw.g = states(:,1:7);
    r.last = times(end);
    if any(r.identifier == ["Precision:Calls","Precision:Time"]), r.status = "budget"; end
end
[msg,id] = lastwarn;
r.warning = string(id)+" "+string(msg);
r.calls = calls;
r.parts = calls;
r.seconds = toc(clock);

    function dy = rhs(t,y)
        calls = calls+1;
        if calls > c.calls, error('Precision:Calls','达到函数调用上限。'); end
        if mod(calls,32) == 0 && toc(clock) > c.seconds
            error('Precision:Time','达到本组运行时间上限。');
        end
        dy = my_system(t,y,7,p,v.r1,v.r2,v.lambda1,v.lambda2,v.a,v.p,v.q,c.delta);
        if any(~isfinite(dy)), error('Precision:Nonfinite','原 VFCR 产生非有限导数。'); end
    end

    function stop = output(t,y,flag)
        stop = 0;
        if ~isempty(flag), return, end
        take = t > times(end);
        times = [times;t(take).'];
        states = [states;y(:,take).'];
    end
end

function m = measure(r,ref,p)
% 只用同一残差和参考解评价 VFCR。
m.xd = vecnorm(r.g(:,1:2)-ref.xd,2,2);
m.x = vecnorm(r.g(:,1:2)-ref.x,2,2);
m.g = vecnorm(r.g-ref.g,2,2);
m.residual = zeros(numel(r.t),1);
m.equality = m.residual;
m.inequality = m.residual;
for k = 1:numel(r.t)
    t = r.t(k);
    g = r.g(k,:).';
    xi = pp_parts(t,g,p,r.c.delta);
    m.residual(k) = norm(xi);
    m.equality(k) = norm(p.P(t)*g(1:2)-p.u(t));
    m.inequality(k) = max([0;p.Q(t)*g(1:2)-p.v(t)]);
end
end

function row = summary(r)
% 统一保存完整和失败组的公共指标。
step = r.c.h;
correct = r.c.correct;
if r.model == "VFCR-S", step = NaN; correct = "none"; end
names = {'xd','x','g','xdpeak','gpeak','residual','residualrmse', ...
    'equality','inequality','tc4','tc6'};
values = nan(1,numel(names));
if r.complete
    m = r.metrics;
    take = r.t >= r.c.tail;
    t = r.t(take);
    rms = @(e) sqrt(trapz(t,e(take).^2)/(t(end)-t(1)));
    values = [rms(m.xd),rms(m.x),rms(m.g),max(m.xd(take)),max(m.g(take)), ...
        max(m.residual(take)),rms(m.residual),max(m.equality(take)), ...
        max(m.inequality(take)),crossing(r.t,m.residual,1e-4), ...
        crossing(r.t,m.residual,1e-6)];
end
row = cell2table([{r.group,r.label,r.model,r.c.rtol,r.c.atol,step,correct, ...
    r.complete,r.status,r.last,r.calls,r.parts,r.seconds},num2cell(values),{r.reason}], ...
    'VariableNames',[{'group','label','model','rtol','atol','h','correct', ...
    'complete','status','last','calls','parts','seconds'},names,{'reason'}]);
end

function time = crossing(t,v,b)
% 记录最后一次超阈值后的采样时刻。
id = find(v > b,1,'last');
if isempty(id), time = t(1);
elseif id == numel(t), time = Inf;
else, time = t(id+1);
end
end

function code = digest(path)
% 记录实际使用文件的校验值。
file = fopen(path,'rb');
assert(file >= 0);
guard = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8');
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes);
code = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
