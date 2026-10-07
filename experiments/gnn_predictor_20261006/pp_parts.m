function [xi,J,blocks] = pp_parts(t,g,problem,delta)
% 稳定计算 PFB 残差和雅可比。
G = problem.G(t);
S = (G+G.')/2;
h = problem.h(t);
P = problem.P(t);
u = problem.u(t);
Q = problem.Q(t);
v = problem.v(t);
h = h(:);
u = u(:);
v = v(:);
g = g(:);
n = size(G,1);
m = size(P,1);
w = size(Q,1);
blocks = [n,m,w];
if size(G,2) ~= n || size(P,2) ~= n || size(Q,2) ~= n || ...
        numel(h) ~= n || numel(u) ~= m || numel(v) ~= w || numel(g) ~= n+m+w
    error('PP:Size','题目矩阵和状态长度不匹配。');
end
if ~isnumeric(delta) || ~isreal(delta) || ...
        any(~isfinite(delta(:))) || any(delta(:) <= 0)
    error('PP:Delta','delta 必须为有限正数。');
end
if isscalar(delta)
    delta = delta*ones(w,1);
elseif numel(delta) == w
    delta = delta(:);
else
    error('PP:DeltaSize','delta 必须为标量或不等式个数长度的向量。');
end
x = reshape(g(1:n),n,1);
lambda = reshape(g(n+1:n+m),m,1);
mu = reshape(g(n+m+1:end),w,1);
omega = v-Q*x;
sigma = hypot(hypot(omega,mu),sqrt(delta));
joint = omega+mu;
phi = joint-sigma;
take = joint >= 0;
phi(take) = (2*omega(take).*mu(take)-delta(take))./(joint(take)+sigma(take));
alpha = 1-omega./sigma;
beta = 1-mu./sigma;
take = omega > 0;
alpha(take) = (hypot(mu(take),sqrt(delta(take)))./sigma(take)).^2 ...
    ./(1+omega(take)./sigma(take));
take = mu > 0;
beta(take) = (hypot(omega(take),sqrt(delta(take)))./sigma(take)).^2 ...
    ./(1+mu(take)./sigma(take));
xi = [S*x+h+P.'*lambda+Q.'*mu;P*x-u;phi];
J = [S,P.',Q.';P,zeros(m,m),zeros(m,w); ...
    -diag(alpha)*Q,zeros(w,m),diag(beta)];
if ~isreal(xi) || ~isreal(J) || any(~isfinite(xi)) || any(~isfinite(J),'all')
    error('PP:Nonfinite','残差或雅可比包含非有限实数。');
end
end
