function gdot = oznn_system(t, g, l, problem, lambda, delta)
%OZNN_SYSTEM 实现论文公式（14）代入公式（15）后的模型。

if numel(g) ~= l
    error('OZNN:WrongSize', 'g 的长度必须是 %d。', l);
end

G = problem.G(t);
dG = problem.dG(t);
h = problem.h(t);
dh = problem.dh(t);
P = problem.P(t);
dP = problem.dP(t);
u = problem.u(t);
du = problem.du(t);
Q = problem.Q(t);
dQ = problem.dQ(t);
v = problem.v(t);
dv = problem.dv(t);

n = size(G,1);
m = size(P,1);
w = size(Q,1);

x = g(1:n);
mu2 = g(n+m+1:n+m+w);

if isscalar(delta)
    delta = delta*ones(w,1);
else
    delta = delta(:);
end

omega = v - Q*x;
sigma = sqrt(omega.^2 + mu2.^2 + delta);
k1 = diag(omega./sigma);
k2 = diag(mu2./sigma);

Iw = eye(w);
Zmm = zeros(m,m);
Zmw = zeros(m,w);
Zwm = zeros(w,m);
Zww = zeros(w,w);

H = [G, P.', Q.'; ...
     P, Zmm, Zmw; ...
    -Q, Zwm, Iw];

theta = [h; -u; v-sigma];

J = [G, P.', Q.'; ...
     P, Zmm, Zmw; ...
     (k1-Iw)*Q, Zwm, Iw-k2];

M = [dG, dP.', dQ.'; ...
     dP, Zmm, Zmw; ...
     (k1-Iw)*dQ, Zwm, Zww];

varsigma = [dh; -du; dv-k1*dv];
xi = H*g + theta;

if isfield(problem,'noise')
    noise = problem.noise(t);
    noise = noise(:);
else
    noise = zeros(l,1);
end

if isscalar(noise)
    noise = noise*ones(l,1);
end

if numel(noise) ~= l
    error('OZNN:WrongNoise', 'noise 的长度必须是 %d。', l);
end

rhs = -M*g - varsigma - lambda*xi + noise;
gdot = J \ rhs;
end
