%% run_example1_fig2.m
% 复现论文 Example 1 的 Fig. 2
% (a) x1(t), (b) x2(t), (c) Error

clear all
close all
clc
rng(42);  % 固定随机种子，使结果可重复

%% 定义 Example 1 的 problem 结构体（严格按照论文）
% TVQP 问题:
% min  0.5*(0.25*sin(t)+1)*x1^2 + 0.5*(0.25*sin(t)+1)*x2^2 
%       + 0.5*cos(t)*x1*x2 + sin(3t)*x1 + cos(3t)*x2
% s.t. sin(4t)*x1 + cos(4t)*x2 = cos(2t)
%      -1e8 <= x1 <= 1e8, -1e8 <= x2 <= 1e8

problem.G = @(t) [0.25*sin(t)+1, 0.5*cos(t);
                  0.5*cos(t), 0.25*sin(t)+1];
              
problem.dG = @(t) [0.25*cos(t), -0.5*sin(t);
                   -0.5*sin(t), 0.25*cos(t)];

problem.h = @(t) [sin(3*t); cos(3*t)];
problem.dh = @(t) [3*cos(3*t); -3*sin(3*t)];

% 论文中 P(t) = [sin(4t), cos(4t)]
problem.P = @(t) [sin(4*t), cos(4*t)];
problem.dP = @(t) [4*cos(4*t), -4*sin(4*t)];

problem.u = @(t) cos(2*t);
problem.du = @(t) -2*sin(2*t);

% 不等式约束: Q(t) = [I2; -I2], v(t) = [1e8; 1e8; 1e8; 1e8]
I2 = eye(2);
problem.Q = @(t) [I2; -I2];
problem.dQ = @(t) zeros(4,2);

INF = 1e8;
problem.v = @(t) [INF; INF; INF; INF];
problem.dv = @(t) zeros(4,1);

%% 状态维度
n = 2;  % x 的维度
m = 1;  % 等式约束个数
w = 4;  % 不等式约束个数
l = n + m + w;  % 总状态维度

%% VFCR-ZNN 模型参数（严格按照论文 Section V）
r1 = 1;
r2 = 1;
lambda1 = 1;
lambda2 = 1;
a = 1;
p = 0.5;
q = 0.5;
delta = 0.0001;

fprintf('VFCR-ZNN 模型参数:\n');
fprintf('  r1 = %d, r2 = %d\n', r1, r2);
fprintf('  p = %.1f, q = %.1f, a = %d\n', p, q, a);
fprintf('  lambda1 = %d, lambda2 = %d\n', lambda1, lambda2);
fprintf('  delta = %g\n\n', delta);

%% 生成多组随机初始值（严格按照论文要求）
% x(0) 从 [0, 1] 随机选取
% mu1(0) 从 [0, 10] 随机选取  
% mu2(0) 从 [0, 100] 和 [0, 1000] 随机选取

num_trials = 5;  % 用5条曲线来展示

g0_list = cell(num_trials, 1);
init_descriptions = cell(num_trials, 1);

for i = 1:num_trials
    % x0: [0, 1]
    x0 = rand(2, 1);
    
    % mu1_0: [0, 10]
    mu1_0 = 10 * rand();
    
    % mu2_0: 交替从 [0, 100] 和 [0, 1000] 选取
    if mod(i, 2) == 1
        mu2_0 = 100 * rand(4, 1);
        init_descriptions{i} = sprintf('mu2∈[0,100]');
    else
        mu2_0 = 1000 * rand(4, 1);
        init_descriptions{i} = sprintf('mu2∈[0,1000]');
    end
    
    g0 = [x0; mu1_0; mu2_0];
    g0_list{i} = g0;
    
    fprintf('初始值 %d: x0=[%.3f, %.3f], mu1=%.3f, %s\n', ...
            i, x0(1), x0(2), mu1_0, init_descriptions{i});
end

%% 定义颜色和线型
colors = {'b', 'r', 'g', 'm', 'k'};
line_styles = {'-', '-', '-', '-', '--'};

%% 求解时间
tspan = [0, 5];
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

%% 存储所有结果
T_all = cell(num_trials, 1);
x1_all = cell(num_trials, 1);
x2_all = cell(num_trials, 1);
error_all = cell(num_trials, 1);

%% 对每组初始值进行求解
fprintf('\n开始求解...\n');
for i = 1:num_trials
    g0 = g0_list{i};
    z0 = zeros(l, 1);
    Y0 = [g0; z0];
    
    fprintf('  求解初始值组 %d...\n', i);
    
    % 求解 ODE
    [T, Y] = ode15s(@(t,Y) znn_ode(t, Y, l, problem, r1, r2, lambda1, lambda2, a, p, q, delta), ...
                     tspan, Y0, options);
    
    % 提取结果
    g_sol = Y(:, 1:l);
    x1_sol = g_sol(:, 1);
    x2_sol = g_sol(:, 2);
    
    % 计算误差 ||H*g + theta||
    error = zeros(length(T), 1);
    for k = 1:length(T)
        t = T(k);
        gk = g_sol(k, :)';
        
        Gk = problem.G(t);
        Pk = problem.P(t);
        Qk = problem.Q(t);
        hk = problem.h(t);
        uk = problem.u(t);
        vk = problem.v(t);
        
        xk = gk(1:n);
        mu1k = gk(n+1:n+m);
        mu2k = gk(n+m+1:end);
        
        omegak = vk - Qk*xk;
        sigma_k = sqrt(omegak.^2 + mu2k.^2 + delta);
        
        Hk = [Gk, Pk', Qk';
              Pk, zeros(m,m), zeros(m,w);
              -Qk, zeros(w,m), eye(w)];
        
        theta_k = [hk; -uk; vk - sigma_k];
        
        error(k) = norm(Hk*gk + theta_k);
    end
    
    % 存储结果
    T_all{i} = T;
    x1_all{i} = x1_sol;
    x2_all{i} = x2_sol;
    error_all{i} = error;
    
    % 计算收敛时间（误差降到 1e-6 以下）
    idx = find(error < 1e-6, 1);
    if ~isempty(idx)
        fprintf('    收敛时间 ≈ %.4f s\n', T(idx));
    end
end

%% 绘制 Fig. 2 (a) x1(t)
figure('Position', [100, 100, 1200, 350]);

subplot(1,3,1);
hold on;
for i = 1:num_trials
    plot(T_all{i}, x1_all{i}, colors{i}, 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('x_1(t)');
title('(a) x_1(t)');
grid on;
legend(arrayfun(@(i) sprintf('Init %d', i), 1:num_trials, 'UniformOutput', false), ...
       'Location', 'best');
xlim([0, 5]);

%% 绘制 Fig. 2 (b) x2(t)
subplot(1,3,2);
hold on;
for i = 1:num_trials
    plot(T_all{i}, x2_all{i}, colors{i}, 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('x_2(t)');
title('(b) x_2(t)');
grid on;
legend(arrayfun(@(i) sprintf('Init %d', i), 1:num_trials, 'UniformOutput', false), ...
       'Location', 'best');
xlim([0, 5]);

%% 绘制 Fig. 2 (c) Error (对数坐标)
subplot(1,3,3);
hold on;
for i = 1:num_trials
    semilogy(T_all{i}, error_all{i}, colors{i}, 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('||H(t)g(t) + \vartheta(t)||');
title('(c) Error');
grid on;
legend(arrayfun(@(i) sprintf('Init %d', i), 1:num_trials, 'UniformOutput', false), ...
       'Location', 'southwest');
xlim([0, 5]);
ylim([1e-10, 1e5]);

% 标记理论收敛时间 T1 ≈ 2.7s
xline(2.7, '--k', 'T_1 ≈ 2.7s', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');

sgtitle('Fig. 2: VFCR-ZNN模型在Example 1上的求解结果（含不同初始值）');

%% 额外：检查等式约束是否满足
figure('Position', [100, 100, 800, 400]);

constraint_error = zeros(length(T_all{1}), 1);
for k = 1:length(T_all{1})
    t = T_all{1}(k);
    xk = [x1_all{1}(k); x2_all{1}(k)];
    Pk = problem.P(t);
    uk = problem.u(t);
    constraint_error(k) = abs(Pk*xk - uk);
end

plot(T_all{1}, constraint_error, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('|P(t)x(t) - u(t)|');
title('等式约束误差');
grid on;
set(gca, 'YScale', 'log');

fprintf('\n等式约束最终误差: %.2e\n', constraint_error(end));

%% ============ 辅助函数 ============

function Y_dot = znn_ode(t, Y, l, problem, r1, r2, lambda1, lambda2, a, p, q, delta)
    g = Y(1:l);
    z = Y(l+1:2*l);
    
    [g_dot, z_dot] = my_system(t, g, z, problem, r1, r2, lambda1, lambda2, a, p, q, delta);
    
    Y_dot = [g_dot; z_dot];
end