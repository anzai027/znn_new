function dz = model(t,z,G,cfg)
n = size(G,1);
l = 2*n+1;
g = z(1:l);
y = z(l+1:end);
[J,xi] = parts(g,G,cfg);
rho = exp(cfg.l1*acot(t)+cfg.l2);
first = phi(xi,cfg);
second = phi(xi+cfg.r1*y,cfg);
rhs = -cfg.r1*rho*first-cfg.r2*second;
dg = J\rhs;
dy = rho*first;
dz = [dg;dy];
end
