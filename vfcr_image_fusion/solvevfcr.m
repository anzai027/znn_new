function sol = solvevfcr(G,cfg,seed)
n = size(G,1);
l = 2*n+1;
rng(seed,'twister')
x0 = rand(n,1);
x0 = x0/sum(x0);
g0 = [x0;0;zeros(n,1)];
y0 = zeros(l,1);
z0 = [g0;y0];
opt = odeset('RelTol',cfg.rel,'AbsTol',cfg.abs,'MaxStep',cfg.step);
tic
[t,z] = ode15s(@(t,z) model(t,z,G,cfg),cfg.span,z0,opt);
runtime = toc;
if t(end) < cfg.span(end)-1e-9
    error('VFCR:Incomplete','ODE stopped at t=%.6g.',t(end))
end
e = zeros(numel(t),1);
for i = 1:numel(t)
    [~,xi] = parts(z(i,1:l)',G,cfg);
    e(i) = norm(xi,2);
end
x = z(end,1:n)';
hit = find(e<=cfg.tol,1);
if isempty(hit)
    settle = NaN;
else
    settle = t(hit);
end
sol.t = t;
sol.z = z;
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
