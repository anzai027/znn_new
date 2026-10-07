function dg = gsmooth(t, g, l, problem, gamma, p, delta, eps)
% 使用平滑分母计算状态导数。
% 状态 g=[x;mu1;mu2]，参数需满足 gamma>0、0<p<2、delta>0。
% 固定 delta>0 时求解的是扰动后的 KKT/PFB 方程。
% G 可以不对称，二次目标使用 G 的对称部分。
% 若要保证目标函数凸，G 的对称部分应半正定。
% 平滑后不能直接沿用原有限时间结论。

if ~isnumeric(gamma) || ~isscalar(gamma) || ~isreal(gamma) || ...
        ~isfinite(gamma) || gamma <= 0
    error('GNN:Gamma', 'gamma 必须是正数。');
end

if ~isnumeric(p) || ~isscalar(p) || ~isreal(p) || ...
        ~isfinite(p) || p <= 0 || p >= 2
    error('GNN:Exponent', 'p 必须满足 0 < p < 2。');
end

if ~isnumeric(delta) || ~isreal(delta) || ...
        any(~isfinite(delta(:))) || any(delta(:) <= 0)
    error('GNN:Delta', 'delta 的每个分量都必须是正数。');
end

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

n = size(G, 1);
m = size(P, 1);
w = size(Q, 1);

if ~isscalar(l) || l ~= n+m+w || numel(g) ~= l
    error('GNN:State', 'l 和 g 的长度必须等于 n+m+w。');
end

if isscalar(delta)
    delta = delta*ones(w, 1);
elseif numel(delta) == w
    delta = delta(:);
else
    error('GNN:DeltaSize', 'delta 必须是标量或长度为 w 的向量。');
end

x = g(1:n);
mu1 = g(n+1:n+m);
mu2 = g(n+m+1:l);

omega = v-Q*x;
sigma = sqrt(omega.^2+mu2.^2+delta);

xi = [S*x+h+P.'*mu1+Q.'*mu2; ...
      P*x-u; ...
      omega+mu2-sigma];

k1 = diag(omega./sigma);
k2 = diag(mu2./sigma);
I = eye(w);

J = [S, P.', Q.'; ...
     P, zeros(m,m), zeros(m,w); ...
     (k1-I)*Q, zeros(w,m), I-k2];

s = J.'*xi;
r = norm(s, 2);

if r == 0
    if any(xi ~= 0)
        error('GNN:StationaryResidual', '梯度为零，但 KKT/PFB 残差不为零。');
    end
    % 零点取零速度，这只是 M1 式的一种数值取值。
    dg = zeros(l, 1);
    return
end

dg = -gamma*s/hypot(r,eps)^p;
end
