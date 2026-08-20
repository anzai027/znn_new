function sol = solvevfcr(G,cfg,seed)
n = size(G,1);
l = 2*n+1;
problem = makeprob(G);
rng(seed,'twister')
x0 = rand(n,1);
x0 = x0/sum(x0);
g0 = [x0;0;zeros(n,1)];
y0 = zeros(l,1);
Y0 = [g0;y0];
opt = odeset('RelTol',cfg.rel,'AbsTol',cfg.abs,'MaxStep',cfg.step);
fun = @(t,Y) my_system(t,Y,l,problem,cfg.r1,cfg.r2, ...
    cfg.l1,cfg.l2,cfg.a,cfg.p,cfg.q,cfg.delta);
tic
[t,Y] = ode15s(fun,cfg.span,Y0,opt);
runtime = toc;
if t(end) < cfg.span(end)-1e-9
    error('VFCR:Incomplete','ODE stopped at t=%.6g.',t(end))
end
e = zeros(numel(t),1);
for i = 1:numel(t)
    [~,xi] = parts(Y(i,1:l)',G,cfg);
    e(i) = norm(xi,2);
end
x = Y(end,1:n)';
hit = find(e<=cfg.tol,1);
if isempty(hit)
    settle = NaN;
else
    settle = t(hit);
end
sol.t = t;
sol.z = Y;
sol.e = e;
sol.x = x;
sol.g0 = g0;
sol.y0 = y0;
sol.run = runtime;
sol.stop = settle;
sol.final = e(end);
sol.sum = sum(x);
sol.min = min(x);
end

function problem = makeprob(G)
n = size(G,1);
problem.G = @(t) G;
problem.dG = @(t) zeros(n);
problem.h = @(t) zeros(n,1);
problem.dh = @(t) zeros(n,1);
problem.P = @(t) ones(1,n);
problem.dP = @(t) zeros(1,n);
problem.u = @(t) 1;
problem.du = @(t) 0;
problem.Q = @(t) -eye(n);
problem.dQ = @(t) zeros(n);
problem.v = @(t) zeros(n,1);
problem.dv = @(t) zeros(n,1);
end
