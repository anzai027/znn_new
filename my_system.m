function [g_dot, z_dot] = my_system( ...
    t, g, z, problem, ...
    r1, r2, lambda1, lambda2, a, p, q, delta)

%% 读取具体题目的矩阵
G  = problem.G(t);
dG = problem.dG(t);

h  = problem.h(t);
dh = problem.dh(t);

P  = problem.P(t);
dP = problem.dP(t);

u  = problem.u(t);
du = problem.du(t);

Q  = problem.Q(t);
dQ = problem.dQ(t);

v  = problem.v(t);
dv = problem.dv(t);

%% 每一个下面的0矩阵要和上面矩阵的列数相同
n = size(G, 1);
m = size(P, 1);
w = size(Q, 1);

%% 从 g 中取出变量
x   = g(1:n);
mu1 = g(n+1:n+m);
mu2 = g(n+m+1:n+m+w);

%% omega 和 sigma
omega = v - Q*x;

if isscalar(delta)
    delta = delta * ones(w,1);
else
    delta = delta(:);
end

sigma = sqrt(omega.^2 + mu2.^2 + delta);

%% k1 和 k2
k1 = diag(omega ./ sigma);
k2 = diag(mu2 ./ sigma);

%% 零矩阵和单位矩阵
Iw = eye(w);

Zmm = zeros(m,m);
Zmw = zeros(m,w);
Zwm = zeros(w,m);
Zww = zeros(w,w);

%% H(t)
H = [
     G,   P.',          Q.';
     P,   Zmm,          Zmw;
    -Q,   Zwm,          Iw
];

%% vartheta(t)
theta = [
     h;
    -u;
     v - sigma
];

%% J(t)
J = [
    G,                P.',    Q.';
    P,                Zmm,    Zmw;
    (k1-Iw)*Q,        Zwm,    Iw-k2
];

%% M(t)
M = [
    dG,               dP.',   dQ.';
    dP,               Zmm,    Zmw;
    (k1-Iw)*dQ,       Zwm,    Zww
];

%% varsigma(t)
varsigma = [
     dh;
    -du;
     dv - k1*dv
];

%% rho(t)
rho = exp(lambda1 * acot(t) + lambda2);

%% 误差
xi = H*g + theta;

%% 激活函数
phi_xi = Phi(xi, a, p, q);

%% 积分状态导数
z_dot = rho * phi_xi;

%% 公式（20）
RHS = -M*g ...
      - varsigma ...
      - r1*rho*phi_xi ...
      - r2*Phi(xi + r1*z, a, p, q);

g_dot = J \ RHS;

end


function y = Phi(x, a, p, q)

y = a ...
    .* exp(abs(x).^q) ...
    .* abs(x).^p ...
    .* sign(x);

end