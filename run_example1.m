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

problem = example1_problem();

%% 维数和求解器
t0 = 0;
n = size(problem.G(t0),1);
m = size(problem.P(t0),1);
w = size(problem.Q(t0),1);
l = n + m + w;             % 示例1中l=7

t_plot = linspace(0,10,2001).'; %0-10秒遵循论文图，采样2001个点
initial_ranges = [1,10,100,1000];
case_count = numel(initial_ranges);

% 论文未公布随机种子，固定种子1确保结果可重复，但初始瞬态不必与论文完全相同
rng(1,'twister')

ode_options = odeset( ... % ODE求解器和容差
    'RelTol',1e-5, ...
    'AbsTol',1e-7, ...
    'MaxStep',0.02);

X = cell(case_count,1);
error_norm = zeros(numel(t_plot),case_count);

for k = 1:case_count %循环运行四次，每次采用一种初值范围

    % 1.从论文给出的范围随机生成7维g(0)
    g0 = initial_ranges(k)*rand(l,1);

    % 2.积分状态从零开始，组成14维ODE初值
    z0 = zeros(l,1);
    Y0 = [g0;z0];

    % 3.把额外参数固定进去，使函数符合fun(t,Y)形式
    ode_fun = @(t,Y) my_system( ...
        t,Y,l,problem, ...
        r1,r2,lambda1,lambda2,a,p,q,delta);

    % 4.求解式(20)，得到整个时间段的状态轨迹
    [T,Y] = ode15s(ode_fun,t_plot,Y0,ode_options);

    % 5.Y的前两列是x1(t)、x2(t)
    X{k} = Y(:,1:n); %用元胞数组储存完整的x1 x2轨迹

    % 6.在每个采样时刻计算残差范数
    for j = 1:numel(T)
        g_current = Y(j,1:l).';
        xi = vfcr_residual( ...
            T(j),g_current,problem,delta);

        error_norm(j,k) = norm(xi,2); %2001×1列向量
    end

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
