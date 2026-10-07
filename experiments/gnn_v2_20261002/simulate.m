function r = simulate(job,problem,bank,c)
% 运行平滑公式并保存中断轨迹。
l = job.l;
if ~startsWith(job.model,"VFCR")
    problem = rmfield(problem,{'dG','dh','dP','du','dQ','dv'});
    y0 = job.g0;
else
    y0 = [job.g0;zeros(l,1)];
end
if job.model == "VFCR-S"
    problem.phi = @(x) fsmooth(x,c.vfcr,job.eps);
end
hold = zeros(l,1);
if startsWith(job.model,"VFCR") && job.port == "paper"
    problem.noise = @eta;
end
calls = 0;
used = 1;
size0 = 4096;
times = zeros(size0,1);
states = zeros(size0,numel(y0));
times(1) = 0;
states(1,:) = y0.';
r = struct();
r.job = job;
r.bank = bank;
r.t = c.t;
r.y = nan(numel(c.t),numel(y0));
r.y(1,:) = y0.';
r.status = "running";
r.reason = "";
r.identifier = "";
r.log = "";
r.warnings = strings(0,2);
r.stack = struct([]);
timer = tic;
opt = odeset('RelTol',job.rtol,'AbsTol',job.atol, ...
    'MaxStep',job.step,'OutputFcn',@output);
knots = [0;c.t(end)];
sol = [];
if bank.random
    knots = bank.t;
end
for j = 1:numel(knots)-1
    a = knots(j);
    b = knots(j+1);
    if bank.random
        hold = bank.w(j,:).';
    end
    lastwarn('');
    try
        log = evalc('sol=ode15s(@rhs,[a,b],y0,opt);');
        r.log = r.log+string(log);
        [msg,id] = lastwarn;
        if ~isempty(msg)
            r.warnings(end+1,:) = [string(id),string(msg)];
        end
        last = sol.x(end);
        take = c.t >= a-1e-12 & c.t <= last+1e-12;
        query = min(max(c.t(take),a),last);
        r.y(take,:) = deval(sol,query).';
        y0 = sol.y(:,end);
        if last < b-1e-10
            r.status = "solver_failure";
            r.reason = "求解器提前停止，未完成时间区间。";
            r.identifier = string(id);
            break
        end
    catch err
        r.status = "exception";
        if any(string(err.identifier) == ["Compare:Time","Compare:Calls"])
            r.status = "budget";
        end
        r.reason = string(err.message);
        r.identifier = string(err.identifier);
        r.stack = err.stack;
        [msg,id] = lastwarn;
        if ~isempty(msg)
            r.warnings(end+1,:) = [string(id),string(msg)];
        end
        break
    end
end
r.seconds = toc(timer);
r.calls = calls;
r.raw.t = times(1:used);
r.raw.y = states(1:used,:);
r.last = r.raw.t(end);
if r.status == "running"
    r.status = "complete";
end
if r.status ~= "complete" && used > 1
    take = isnan(r.y(:,1)) & c.t <= r.last & c.t >= r.raw.t(1);
    r.y(take,:) = interp1(r.raw.t,r.raw.y,c.t(take),'linear');
end
r.interpolation = "成功轨迹用 deval，中断轨迹只在已接受点之间线性插值，不外推。";
r.complete = r.status == "complete" && all(isfinite(r.y),'all') ...
    && r.last >= c.t(end)-1e-10;
if ~r.complete && r.status == "complete"
    r.status = "nonfinite";
    r.reason = "完整输出包含非有限数，不能视为成功运行。";
end

    function dy = rhs(t,y)
        calls = calls+1;
        if calls > c.calls
            error('Compare:Calls','达到 %d 次函数调用上限。',c.calls);
        end
        if mod(calls,32) == 0 && toc(timer) >= c.seconds
            error('Compare:Time','达到 %.0f 秒运行上限。',c.seconds);
        end
        if job.model == "M1"
            dy = gnn_system(t,y,l,problem,job.gamma,1,c.delta);
        elseif job.model == "M2"
            dy = gsmooth(t,y,l,problem,job.gamma,job.p,c.delta,job.eps);
        elseif any(job.model == ["M3","M4"])
            dy = pg_system(t,y,l,problem,job.rates,job.lambda,job.eps,c.delta,job.model == "M4");
        else
            v = c.vfcr;
            dy = my_system(t,y,l,problem,v.r1,v.r2,v.lambda1, ...
                v.lambda2,v.a,v.p,v.q,c.delta);
        end
        if job.port == "state"
            dy(1:l) = dy(1:l)+eta(t);
        end
        if any(~isfinite(dy))
            error('Compare:Nonfinite','平滑公式产生了非有限状态导数。');
        end
    end

    function value = eta(t)
        switch bank.kind
            case "zero"
                value = zeros(l,1);
            case "constant"
                value = 2*ones(l,1);
            case "linear"
                value = 0.1*t*ones(l,1);
            case "cosine"
                value = cos(t)*ones(l,1);
            case "harmonic"
                value = sin(8*pi*t+2)*ones(l,1);
            otherwise
                value = hold;
        end
    end

    function stop = output(t,y,flag)
        stop = 0;
        if ~isempty(flag)
            return
        end
        for z = 1:numel(t)
            if t(z) <= times(used)
                continue
            end
            used = used+1;
            if used > size(times,1)
                times(end+size0,1) = 0;
                states(end+size0,size(states,2)) = 0;
            end
            times(used) = t(z);
            states(used,:) = y(:,z).';
        end
    end
end
