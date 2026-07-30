function problem = example2_problem()
%EXAMPLE2_PROBLEM  论文 Example 2 的时变二次规划数据。
%
% 与 Example 1 共用：
%   G(t), dG(t), h(t), dh(t), P(t), dP(t), u(t), du(t), Q(t), dQ(t)
%
% 与 Example 1 的核心差别：
%   Example 1 用 v = 1e8*ones(4,1)，近似表示没有有限上下界；
%   Example 2 用 v = 1.2*ones(4,1)，真正限制 -1.2 <= x1,x2 <= 1.2。

I2 = eye(2);
Q0 = [I2; -I2];

%% 目标函数 0.5*x'*G(t)*x + h(t)'*x
problem.G = @(t) [ ...
    0.25*sin(t) + 1, 0.5*cos(t); ...
    0.5*cos(t),     0.25*sin(t) + 1];

problem.dG = @(t) [ ...
    0.25*cos(t), -0.5*sin(t); ...
   -0.5*sin(t),  0.25*cos(t)];

problem.h  = @(t) [sin(3*t); cos(3*t)];
problem.dh = @(t) [3*cos(3*t); -3*sin(3*t)];

%% 等式约束 P(t)*x(t) = u(t)
problem.P  = @(t) [sin(4*t), -cos(4*t)];
problem.dP = @(t) [4*cos(4*t), 4*sin(4*t)];
problem.u  = @(t) cos(2*t);
problem.du = @(t) -2*sin(2*t);

%% 不等式约束 Q*x <= v，即 -1.2 <= x1,x2 <= 1.2
problem.Q  = @(t) Q0; %#ok<NASGU>
problem.dQ = @(t) zeros(4,2); %#ok<NASGU>
problem.v  = @(t) 1.2*ones(4,1); %#ok<NASGU>
problem.dv = @(t) zeros(4,1); %#ok<NASGU>
end
