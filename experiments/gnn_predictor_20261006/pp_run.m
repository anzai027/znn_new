function r = pp_run(model,p,c)
% 独立积分并保存求解器实际接受的轨迹。
names = {'G','h','P','u','Q','v'};
q = struct();
for k = 1:numel(names), q.(names{k}) = p.(names{k}); end
calls = 0;
clock = tic;
times = c.start;
states = c.g0.';
opt = odeset('RelTol',c.rtol,'AbsTol',c.atol,'MaxStep',c.step,'OutputFcn',@output);
r = struct('model',string(model),'c',c,'status',"running",'reason',"", ...
    'identifier',"",'warning',"",'log',"",'complete',false);
lastwarn('');
sol = [];
try
    log = evalc('sol = ode15s(@rhs,[c.start,c.finish],c.g0,opt);');
    r.log = string(log);
    r.sol = sol;
    r.raw.t = sol.x.';
    r.raw.g = sol.y.';
    r.last = sol.x(end);
    r.complete = r.last >= c.finish-1e-10 && all(isfinite(sol.y),'all');
    if r.complete
        r.t = c.t;
        r.g = deval(sol,c.t).';
        r.complete = all(isfinite(r.g),'all');
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
    r.raw.g = states;
    r.last = times(end);
    if startsWith(string(err.identifier),"Predictor:"), r.status = "budget"; end
end
[msg,id] = lastwarn;
r.warning = string(id)+" "+string(msg);
r.calls = calls;
r.seconds = toc(clock);
r.parts = calls*(1+(model ~= "M3"));

    function dg = rhs(t,g)
        calls = calls+1;
        if calls > c.calls, error('Predictor:Calls','达到函数调用上限。'); end
        if mod(calls,32) == 0 && toc(clock) > c.seconds
            error('Predictor:Time','达到本组运行时间上限。');
        end
        dg = pp_system(t,g,q,c,model);
        if any(~isfinite(dg)), error('Predictor:Nonfinite','导数包含非有限数。'); end
    end

    function stop = output(t,g,flag)
        stop = 0;
        if ~isempty(flag), return, end
        take = t > times(end);
        times = [times;t(take).'];
        states = [states;g(:,take).'];
    end
end
