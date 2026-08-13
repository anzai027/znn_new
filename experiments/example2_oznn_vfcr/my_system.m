function Y_dot = my_system( ...
    t, Y, l, problem, ...
    r1, r2, lambda1, lambda2, a, p, q, delta)

g = Y(1:l);
z = Y(l+1:2*l);  %% z是误差积分部分

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

%% 从 g 中取出变量，不用取mu1的原因
x   = g(1:n);
mu2 = g(n+m+1:n+m+w);

%% omega 和 sigma
omega = v - Q*x;

if isscalar(delta)
    delta = delta * ones(w,1);
else
    delta = delta(:);  %%直接到这一行
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
phi_xi = Active(xi, problem, a, p, q);

%% 积分状态导数
z_dot = rho * phi_xi;

%% 外部噪声
if isfield(problem, 'noise')
    noise = problem.noise(t);
    noise = noise(:);
else
    noise = zeros(l,1);    %直接到这一行
end

if isscalar(noise)
    noise = noise * ones(l,1);
end

%% 公式（20）和公式（47）
RHS = -M*g ...
      - varsigma ...
      - r1*rho*phi_xi ...
      - r2*Active(xi + r1*z, problem, a, p, q) ...
      + noise;

g_dot = J \ RHS;  %RHS是式子右边的全部东西，公式20

Y_dot = [g_dot; z_dot];

end


function y = Active(x, problem, a, p, q)

if isfield(problem, 'phi')
    y = problem.phi(x);
else
    y = Phi(x, a, p, q);  %直接到这行
end

end


function y = Phi(x, a, p, q)

y = a ...
    .* exp(abs(x).^q) ...
    .* abs(x).^p ...
    .* sign(x);

end
