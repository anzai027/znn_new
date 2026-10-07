function [xi,J] = parts(t,g,problem,delta)
% 计算共同的残差和雅可比矩阵。
G = problem.G(t);
G = (G+G.')/2;
h = problem.h(t);
P = problem.P(t);
u = problem.u(t);
Q = problem.Q(t);
v = problem.v(t);
x = g(1:2);
mu1 = g(3);
mu2 = g(4:7);
omega = v-Q*x;
sigma = sqrt(omega.^2+mu2.^2+delta);
xi = [G*x+h+P.'*mu1+Q.'*mu2;P*x-u;omega+mu2-sigma];
k1 = diag(omega./sigma);
k2 = diag(mu2./sigma);
J = [G,P.',Q.';P,0,zeros(1,4); ...
    (k1-eye(4))*Q,zeros(4,1),eye(4)-k2];
end
