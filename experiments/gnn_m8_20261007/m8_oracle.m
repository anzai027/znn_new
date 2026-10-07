function [xi,J,ft] = m8_oracle(t,g,p,delta)
% 独立计算残差及固定状态的解析时间偏导。
G = p.G(t);
S = (G+G.')/2;
P = p.P(t);
Q = p.Q(t);
h = p.h(t);
u = p.u(t);
v = p.v(t);
n = size(S,1);
m = size(P,1);
w = size(Q,1);
g = g(:);
x = reshape(g(1:n),n,1);
lambda = reshape(g(n+1:n+m),m,1);
mu = reshape(g(n+m+1:end),w,1);
if isscalar(delta), delta = repmat(delta,w,1); else, delta = delta(:); end
omega = v(:)-Q*x;
sigma = hypot(hypot(omega,mu),sqrt(delta));
joint = omega+mu;
phi = joint-sigma;
take = joint >= 0;
phi(take) = (2*omega(take).*mu(take)-delta(take))./(joint(take)+sigma(take));
alpha = 1-omega./sigma;
beta = 1-mu./sigma;
take = omega > 0;
alpha(take) = (mu(take).^2+delta(take))./(sigma(take).*(sigma(take)+omega(take)));
take = mu > 0;
beta(take) = (omega(take).^2+delta(take))./(sigma(take).*(sigma(take)+mu(take)));
xi = [S*x+h(:)+P.'*lambda+Q.'*mu;P*x-u(:);phi];
J = [S,P.',Q.';P,zeros(m,m),zeros(m,w); ...
    -diag(alpha)*Q,zeros(w,m),diag(beta)];
if nargout < 3, return, end
dG = p.dG(t);
dS = (dG+dG.')/2;
dP = p.dP(t);
dQ = p.dQ(t);
dh = p.dh(t);
du = p.du(t);
dv = p.dv(t);
ft = [dS*x+dh(:)+dP.'*lambda+dQ.'*mu;dP*x-du(:); ...
    alpha.*(dv(:)-dQ*x)];
end
