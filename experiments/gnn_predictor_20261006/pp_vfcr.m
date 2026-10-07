function pp_vfcr()
% 按新模型的密集网格复跑旧 VFCR。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
old = fullfile(root,'experiments','gnn_v2_20261002');
addpath(old,'-begin'); addpath(root,'-begin'); addpath(folder,'-begin');
assert(strcmp(which('my_system'),fullfile(root,'my_system.m')));
assert(strcmp(which('fsmooth'),fullfile(old,'fsmooth.m')));
out = fullfile(folder,'vfcr');
assert(~isfile(fullfile(out,'all.mat')),'已有结果，请保留原数据后使用新目录。');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'run.log')); clean = onCleanup(@() diary('off'));
source = load(fullfile(folder,'data','all.mat'),'c','ref','hash');
c = source.c; ref = source.ref;
v = struct('r1',1,'r2',1,'a',1,'p',0.5,'q',0.5,'lambda1',1,'lambda2',1);
p = example2_problem();
hash.my_system = digest(fullfile(root,'my_system.m'));
hash.fsmooth = digest(fullfile(old,'fsmooth.m'));
hash.problem = digest(fullfile(root,'example2_problem.m'));
hash.driver = digest(fullfile(folder,'pp_vfcr.m'));
results = cell(3,1); stats = table();
for k = 1:3
    ci = c; group = "main"; model = "VFCR-S";
    if k == 2
        ci.rtol = 1e-9; ci.atol = 1e-11; group = "precision";
    elseif k == 3
        model = "VFCR";
    end
    fprintf('开始 %s %s\n',model,group);
    r = simulate(model,p,ci,v);
    [r.metrics,row] = measure(r,ref,p,group);
    results{k} = r; stats = [stats;row];
    save(fullfile(out,sprintf('run_%02d.mat',k)),'r','v','hash','-v7.3');
    writetable(stats,fullfile(out,'metrics.csv'));
    fprintf('完成 %s %s complete=%d last=%.9g calls=%d xd=%.9g gpeak=%.9g residual=%.9g\n', ...
        model,group,r.complete,r.last,r.calls,row.xd,row.gpeak,row.residual);
end
save(fullfile(out,'all.mat'),'results','stats','c','ref','v','hash','-v7.3');
disp(stats);
end

function r = simulate(model,p,c,v)
% 保存完整输出和中断前已接受的状态。
if model == "VFCR-S", p.phi = @(x) fsmooth(x,v,c.eps); end
calls = 0; clock = tic;
times = c.start; states = [c.g0;zeros(7,1)].';
opt = odeset('RelTol',c.rtol,'AbsTol',c.atol,'MaxStep',c.step,'OutputFcn',@output);
r = struct('model',model,'c',c,'status',"running",'complete',false, ...
    'reason',"",'identifier',"",'log',"");
sol = []; lastwarn('');
try
    log = evalc('sol = ode15s(@rhs,[c.start,c.finish],[c.g0;zeros(7,1)],opt);');
    r.log = string(log); r.sol = sol;
    r.raw.t = sol.x.'; r.raw.y = sol.y.';
    r.last = sol.x(end);
    r.complete = r.last >= c.finish-1e-10 && all(isfinite(sol.y),'all');
    if r.complete
        r.t = c.t; r.y = deval(sol,c.t).'; r.g = r.y(:,1:7);
        r.complete = all(isfinite(r.y),'all'); r.status = "complete";
    else
        r.status = "solver_failure";
    end
catch err
    r.status = "exception"; r.reason = string(err.message); r.identifier = string(err.identifier);
    r.raw.t = times; r.raw.y = states; r.last = times(end);
    if startsWith(r.identifier,"VFCR:"), r.status = "budget"; end
end
[msg,id] = lastwarn; r.warning = string(id)+" "+string(msg);
r.calls = calls; r.seconds = toc(clock);

    function dy = rhs(t,y)
        calls = calls+1;
        if calls > c.calls, error('VFCR:Calls','达到函数调用上限。'); end
        if mod(calls,32) == 0 && toc(clock) > c.seconds
            error('VFCR:Time','达到运行时间上限。');
        end
        dy = my_system(t,y,7,p,v.r1,v.r2,v.lambda1,v.lambda2,v.a,v.p,v.q,c.delta);
        if any(~isfinite(dy)), error('VFCR:Nonfinite','导数包含非有限数。'); end
    end

    function stop = output(t,y,flag)
        stop = 0;
        if ~isempty(flag), return, end
        take = t > times(end);
        times = [times;t(take).']; states = [states;y(:,take).'];
    end
end

function [m,row] = measure(r,ref,p,group)
% 用相同残差与时间加权口径计算性能。
values = nan(1,12); m = struct();
if r.complete
    m.xd = vecnorm(r.g(:,1:2)-ref.xd,2,2);
    m.x = vecnorm(r.g(:,1:2)-ref.x,2,2);
    m.g = vecnorm(r.g-ref.g,2,2);
    m.residual = zeros(numel(r.t),1);
    m.equality = m.residual; m.inequality = m.residual;
    for k = 1:numel(r.t)
        t = r.t(k); g = r.g(k,:).';
        xi = pp_parts(t,g,p,r.c.delta);
        m.residual(k) = norm(xi);
        m.equality(k) = norm(p.P(t)*g(1:2)-p.u(t));
        m.inequality(k) = max([0;p.Q(t)*g(1:2)-p.v(t)]);
    end
    take = r.t >= r.c.tail; t = r.t(take);
    rms = @(e) quadrature(t,e(take));
    values = [rms(m.xd),rms(m.x),rms(m.g),max(m.xd(take)),max(m.g(take)), ...
        max(m.residual(take)),max(m.equality(take)),max(m.inequality(take)), ...
        crossing(r.t,m.residual,1e-4),crossing(r.t,m.residual,1e-6), ...
        rms(m.residual),max(abs(r.y(take,8:14)),[],'all')];
end
names = {'xd','x','g','xdpeak','gpeak','residual','equality','inequality','tc4','tc6','residualrmse','integral'};
row = cell2table([{group,r.model,r.c.rtol,r.c.atol,r.complete,r.status,r.last,r.calls,r.seconds}, ...
    num2cell(values),{r.reason}], ...
    'VariableNames',[{'group','model','rtol','atol','complete','status','last','calls','seconds'},names,{'reason'}]);
end

function value = quadrature(t,e)
% 按时间段独立计算均方误差。
value = sqrt(sum(diff(t).*(e(1:end-1).^2+e(2:end).^2)/2)/(t(end)-t(1)));
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
% 保存模型与实验脚本的校验值。
file = fopen(path,'rb'); clean = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8'); md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes); code = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
