function problem = example1_problem(phi_eps)
%EXAMPLE1_PROBLEM  VFCR-ZNN论文示例1的数据。
%
% phi_eps = 0      ：公式（18）的严格NFTAF，计算可能很慢。
% phi_eps = 1e-4   ：用于快速绘图的平滑数值近似。

if nargin < 1
    phi_eps = 1e-4;
end

I2 = eye(2);
Q0 = [I2; -I2];

problem.G  = @(t) [0.25*sin(t) + 1, 0.5*cos(t); ...
                   0.5*cos(t),       0.25*sin(t) + 1];
problem.dG = @(t) [0.25*cos(t), -0.5*sin(t); ...
                  -0.5*sin(t),  0.25*cos(t)];

problem.h  = @(t) [sin(3*t); cos(3*t)];
problem.dh = @(t) [3*cos(3*t); -3*sin(3*t)];

% sin(4t)x1 - cos(4t)x2 = cos(2t)
problem.P  = @(t) [sin(4*t), -cos(4*t)];
problem.dP = @(t) [4*cos(4*t), 4*sin(4*t)];
problem.u  = @(t) cos(2*t);
problem.du = @(t) -2*sin(2*t);

% 论文用正负无穷表示无有限边界，并在实验中将无穷替换为1e8。
problem.Q  = @(t) Q0; %#ok<NASGU>
problem.dQ = @(t) zeros(4,2); %#ok<NASGU>
problem.v  = @(t) 1e8*ones(4,1); %#ok<NASGU>
problem.dv = @(t) zeros(4,1); %#ok<NASGU>

problem.phi_eps = phi_eps;

% 此上限用于保护隐式ODE求解器被拒绝的试探步，80高于[0,1000]初值通常达到的指数。
problem.exp_cap = 80;
end
