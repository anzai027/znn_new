function [xi,J] = parts(t,g,problem,delta)
% 计算任意维数的共同残差和雅可比矩阵。
G = problem.G(t); G = (G+G.')/2;
h = problem.h(t); P = problem.P(t); u = problem.u(t);
Q = problem.Q(t); v = problem.v(t);
n = size(G,1); m = size(P,1); w = size(Q,1);
x = g(1:n); mu1 = g(n+1:n+m); mu2 = g(n+m+1:n+m+w);
omega = v-Q*x;
sigma = sqrt(omega.^2+mu2.^2+delta);
xi = [G*x+h+P.'*mu1+Q.'*mu2; P*x-u; omega+mu2-sigma];
k1 = diag(omega./sigma); k2 = diag(mu2./sigma);
J = [G,P.',Q.'; P,zeros(m),zeros(m,w); ...
    (k1-eye(w))*Q,zeros(w,m),eye(w)-k2];
end
