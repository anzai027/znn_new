function xi = vfcr_residual(t, g, problem, delta)

%读取题目系数
G = problem.G(t);
h = problem.h(t);
P = problem.P(t);
u = problem.u(t);
Q = problem.Q(t);
v = problem.v(t);

%确定维数
n = size(G,1);
m = size(P,1);
w = size(Q,1);

%拆开g
x   = g(1:n); %2维
mu1 = g(n+1:n+m); %1维
mu2 = g(n+m+1:n+m+w); %4维

%把δ变成4维
if isscalar(delta) %如果传进来的δ是标量，自动拓展成四维列向量
    delta = delta*ones(w,1);
else %如果传进来的δ已经是向量，统一整理成列向量
    delta = delta(:);
end

%计算PFB残差
omega = v - Q*x;
sigma = sqrt(omega.^2 + mu2.^2 + delta);

xi = [G*x + h + P.'*mu1 + Q.'*mu2; ...
      P*x - u; ...
      omega + mu2 - sigma];
     %驻点残差，2维
     %等式约束残差，1维
     %不等式互补残差，4维
    %三类最优条件基本同时满足
end
