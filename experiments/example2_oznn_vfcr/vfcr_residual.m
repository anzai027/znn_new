function xi = vfcr_residual(t, g, problem, delta)
%VFCR_RESIDUAL 计算完整的 KKT/PFB 残差。

G = problem.G(t);
h = problem.h(t);
P = problem.P(t);
u = problem.u(t);
Q = problem.Q(t);
v = problem.v(t);

n = size(G,1);
m = size(P,1);
w = size(Q,1);

x = g(1:n);
mu1 = g(n+1:n+m);
mu2 = g(n+m+1:n+m+w);

if isscalar(delta)
    delta = delta*ones(w,1);
else
    delta = delta(:);
end

omega = v - Q*x;
sigma = sqrt(omega.^2 + mu2.^2 + delta);

xi = [G*x + h + P.'*mu1 + Q.'*mu2; ...
      P*x - u; ...
      omega + mu2 - sigma];
end
