clearvars
close all
clc

r1 = 1;
r2 = 1;
p = 0.5;
q = 0.5;
a = 1;
lambda1 = 1;
lambda2 = 1;
delta = 1e-4;

alpha = 4;
beta = 0.4;
activation_names = {'BSAF','SBPAF','NSBPAF','NFTAF'};
activation_count = numel(activation_names);

problem_base = example2_problem();
t0 = 0;
n = size(problem_base.G(t0),1);
m = size(problem_base.P(t0),1);
w = size(problem_base.Q(t0),1);
l = n + m + w;               

t_plot = linspace(0,10,2001).';
rng(1,'twister')

g0 = rand(l,1);                    % 论文：g(0)从[0,1]随机选取
z0 = zeros(l,1);                   % 历史积分从零开始
Y0 = [g0;z0];                      % 14维完整ODE初值

ode_options = odeset( ...
    'RelTol',1e-5, ...
    'AbsTol',1e-7, ...
    'MaxStep',0.02);

% 第1维：2001个时刻；第2维：四种激活函数；第3维：两种噪声情形
% noise_case=1为无噪声，noise_case=2为余弦噪声
error_norm = zeros(numel(t_plot),activation_count,2);

%% 分别计算Fig.3(a)无噪声与Fig.3(b)余弦噪声
for noise_case = 1:2
    for k = 1:activation_count
        problem = problem_base;
        activation_name = activation_names{k};

        % 每次循环只选择一个激活函数
        problem.phi = @(e) vfcr_phi( ...
            e,activation_name,a,p,q,alpha,beta);

        if noise_case == 1
            % Fig.3(a)：无噪声
            problem.noise = @(t) zeros(l,1); %#ok<NASGU>
        else
            % Fig.3(b)：论文写作zeta_i(t)=cos(t)，即7个分量相同
            problem.noise = @(t) cos(t)*ones(l,1);
        end

        ode_fun = @(t,Y) my_system( ...
            t,Y,l,problem, ...
            r1,r2,lambda1,lambda2,a,p,q,delta);

        [T,Y] = ode15s(ode_fun,t_plot,Y0,ode_options);

        % 每个采样时刻重新计算7维xi，再取二范数得到一个误差标量
        for j = 1:numel(T)
            g_current = Y(j,1:l).';
            xi = vfcr_residual(T(j),g_current,problem,delta);
            error_norm(j,k,noise_case) = norm(xi,2);
        end

        fprintf('%s, noise_case=%d, final ||xi||_2=%.3e\n', ...
            activation_name,noise_case,error_norm(end,k,noise_case));
    end
end

%% 6. 绘制论文Fig.3对应的两个子图
colors = [ ...
    0.25,0.50,0.70; ...
    0.90,0.45,0.15; ...
    0.95,0.65,0.10; ...
    0.90,0.10,0.10];
styles = {'-',':','--','-'};

figure('Color','w','Position',[120,160,1050,420])

subplot(1,2,1)
hold on
for k = 1:activation_count
    plot(t_plot,error_norm(:,k,1),styles{k}, ...
        'Color',colors(k,:),'LineWidth',1.2, ...
        'DisplayName',activation_names{k});
end
xlim([0,10])
ylim([0,2])
xlabel('t (s)')
ylabel('||\xi(t)||_2')
title('(a) 无噪声')
legend('Location','northeast')
box on

subplot(1,2,2)
hold on
for k = 1:activation_count
    plot(t_plot,error_norm(:,k,2),styles{k}, ...
        'Color',colors(k,:),'LineWidth',1.2, ...
        'DisplayName',activation_names{k});
end
xlim([0,10])
ylim([0,2])
xlabel('t (s)')
ylabel('||\xi(t)||_2')
title('(b) 余弦噪声 \zeta_i(t)=cos(t)')
legend('Location','northeast')
box on

% 论文Fig.3(a)中的局部放大图。
axes('Position',[0.245,0.54,0.16,0.23])
hold on
for k = 1:activation_count
    plot(t_plot,error_norm(:,k,1),styles{k}, ...
        'Color',colors(k,:),'LineWidth',1.0);
end
xlim([1,2])
ylim([0,0.1])
box on

sgtitle('VFCR-ZNN Example 2：四种激活函数误差对比')
savefig('example2_fig3_results.fig')
