function r = bd_run(model,p,c)
% 分段启动因果差分并保留失败组的数据。
q = p;
if model ~= "VFCR-S"
    q = rmfield(p,{'dG','dh','dP','du','dQ','dv'});
end
y = c.g0;
if startsWith(model,"VFCR"), y = [y;zeros(7,1)]; end
calls = 0; parts = 0; timer = tic; sols = cell(3,1);
ends = [c.start,c.start+c.h,c.start+2*c.h,c.finish];
opt = odeset('RelTol',c.rtol,'AbsTol',c.atol,'MaxStep',c.step);
r = struct('model',model,'c',c,'complete',false,'status',"running", ...
    'last',c.start,'reason',"",'log',"",'sols',{sols});
try
    for j = 1:3
        sol = []; log = evalc('sol = ode15s(@rhs,ends(j:j+1),y,opt);');
        r.log = r.log+string(log); r.sols{j} = sol; r.last = sol.x(end);
        if r.last < ends(j+1)-1e-12, error('BD:Solver','积分提前结束。'); end
        y = sol.y(:,end);
    end
    r.t = c.t; r.y = zeros(numel(c.t),numel(y));
    for j = 1:3
        take = c.t >= ends(j) & c.t <= ends(j+1);
        if any(take), r.y(take,:) = deval(r.sols{j},c.t(take)).'; end
    end
    r.g = r.y(:,1:7); r.complete = all(isfinite(r.y),'all');
    r.status = "complete";
catch err
    r.status = "failed"; r.reason = string(err.identifier)+": "+string(err.message);
    r.stack = err.stack;
end
r.calls = calls; r.parts = parts; r.seconds = toc(timer);
    function dy = rhs(t,y)
        calls = calls+1;
        if calls > c.calls || (mod(calls,32)==0 && toc(timer)>c.seconds)
            error('BD:Budget','达到本组计算预算。');
        end
        [dy,a] = bd_system(t,y,q,c,model); parts = parts+a.calls;
    end
end
