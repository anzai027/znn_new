clearvars
close all
clc

%% 论文参数
r1 = 1;
r2 = 1;
p = 0.5;
q = 0.5;
a = 1;
lambda1 = 1;
lambda2 = 1;
delta = 1e-4;

% 严格NFTAF在零点不满足利普希茨条件，可能使自适应ODE求解器采用极小步长。
% 默认使用1e-4以便高效复现曲线，仅在测试原始公式时设为0。
phi_eps = 0;
problem = example1_problem(phi_eps);

%% 维数和求解器
t0 = 0;
n = size(problem.G(t0),1);
m = size(problem.P(t0),1);
w = size(problem.Q(t0),1);
l = n + m + w;             % 示例1中l=7

t_plot = linspace(0,10,2001).';
initial_ranges = [1,10,100,1000];
case_count = numel(initial_ranges);

% 论文未公布随机种子，固定种子可确保结果可重复，但初始瞬态不必与论文完全相同。
rng(1,'twister')

ode_options = odeset( ...
    'RelTol',1e-5, ...
    'AbsTol',1e-7, ...
    'MaxStep',0.02);

X = cell(case_count,1);
error_norm = zeros(numel(t_plot),case_count);
elapsed = zeros(case_count,1);
settling_time = nan(case_count,1);

for k = 1:case_count
    g0 = initial_ranges(k)*rand(l,1);
    z0 = zeros(l,1);       % 历史积分从零开始
    Y0 = [g0;z0];          % ode15s接收一个14维状态

    ode_fun = @(t,Y) my_system( ...
        t, Y, l, problem, ...
        r1, r2, lambda1, lambda2, a, p, q, delta);

    tic
    [T,Y] = ode15s(ode_fun, t_plot, Y0, ode_options);
    elapsed(k) = toc;

    X{k} = Y(:,1:n);

    for j = 1:numel(T)
        xi = vfcr_residual(T(j),Y(j,1:l).',problem,delta);
        error_norm(j,k) = norm(xi,2);
    end

    % 数值收敛时间取残差最后一次高于1e-4之后的首个采样时刻。
    last_bad = find(error_norm(:,k) > 1e-4,1,'last');
    if isempty(last_bad)
        settling_time(k) = T(1);
    elseif last_bad < numel(T)
        settling_time(k) = T(last_bad + 1);
    end

    fprintf(['range [0,%4d]: runtime = %8.3f s, ', ...
             'settling estimate = %8.4f s, final ||xi|| = %.3e\n'], ...
        initial_ranges(k),elapsed(k),settling_time(k),error_norm(end,k));
end

%% 可选独立检查：求解等式约束KKT系统
x_theory = zeros(numel(t_plot),n);
for j = 1:numel(t_plot)
    t = t_plot(j);
    G = problem.G(t);
    P = problem.P(t);
    kkt = [G,P.';P,0];
    reference = kkt\[-problem.h(t);problem.u(t)];
    x_theory(j,:) = reference(1:n).';
end

for k = 1:case_count
    after_two_seconds = t_plot >= 2;
    tracking_error = vecnorm( ...
        X{k}(after_two_seconds,:) - x_theory(after_two_seconds,:),2,2);
    fprintf('range [0,%4d]: max ||x-x*|| after 2 s = %.3e\n', ...
        initial_ranges(k),max(tracking_error));
end

%% 复现图2的三个子图
labels = {'g_1(0)','g_2(0)','g_3(0)','g_4(0)'};
colors = [0.00,0.25,1.00; ...
          0.90,0.25,0.10; ...
          0.95,0.60,0.05; ...
          1.00,0.00,0.00];
styles = {'--',':','-.','-'};

figure('Color','w','Position',[80,120,1450,430])

subplot(1,3,1)
hold on
for k = 1:case_count
    plot(t_plot,X{k}(:,1),styles{k},'Color',colors(k,:), ...
        'LineWidth',1.2,'DisplayName',labels{k});
end
xlim([0,10])
ylim([-10,10])
xlabel('t (s)')
ylabel('x_1(t)')
legend('Location','northeast')
box on

subplot(1,3,2)
hold on
for k = 1:case_count
    plot(t_plot,X{k}(:,2),styles{k},'Color',colors(k,:), ...
        'LineWidth',1.2,'DisplayName',labels{k});
end
xlim([0,10])
ylim([-10,10])
xlabel('t (s)')
ylabel('x_2(t)')
legend('Location','northeast')
box on

subplot(1,3,3)
hold on
for k = 1:case_count
    plot(t_plot,error_norm(:,k),styles{k},'Color',colors(k,:), ...
        'LineWidth',1.2,'DisplayName',labels{k});
end
xlim([0,10])
ylim([0,100])
xlabel('t (s)')
ylabel('||\xi(t)||_2')
legend('Location','northeast')
box on

% 插图对应论文中放大的早期残差。
axes('Position',[0.745,0.56,0.115,0.23])
hold on
for k = 1:case_count
    plot(t_plot,error_norm(:,k),styles{k},'Color',colors(k,:), ...
        'LineWidth',1.0);
end
xlim([0,2])
ylim([0,0.1])
box on

sgtitle('VFCR-ZNN equation (20): Example 1')
