clearvars
close all
clc

folder = fileparts(mfilename('fullpath'));
cd(folder)
clear my_system oznn_system example2_problem vfcr_phi vfcr_residual

file = which('my_system');
need = fullfile(folder,'my_system.m');
assert(strcmpi(file,need),'调用的不是实验文件夹内的 my_system.m。')

r1 = 1;
r2 = 1;
p = 0.5;
q = 0.5;
a = 1;
lambda1 = 1;
lambda2 = 1;
lambda = 10;
delta = 1e-4;
alpha = 4;
beta = 0.4;

problem = example2_problem();
t0 = 0;
n = size(problem.G(t0),1);
m = size(problem.P(t0),1);
w = size(problem.Q(t0),1);
l = n + m + w;

problem.phi = @(x) vfcr_phi(x,'NFTAF',a,p,q,alpha,beta);
problem.noise = @(t) zeros(l,1);

t = linspace(0,10,2001).';
rng(1,'twister')
g0 = rand(l,1);
z0 = zeros(l,1);
y0 = [g0;z0];

opt = odeset('RelTol',1e-5,'AbsTol',1e-7,'MaxStep',0.02);
fv = @(t,y) my_system( ...
    t,y,l,problem,r1,r2,lambda1,lambda2,a,p,q,delta);   %VFCR
fo = @(t,g) oznn_system(t,g,l,problem,lambda,delta);    %OZNN

tic
[tv,yv] = ode15s(fv,t,y0,opt);
timev = toc;

tic
[to,yo] = ode15s(fo,t,g0,opt);
timeo = toc;

ev = zeros(numel(t),1);
eo = zeros(numel(t),1);
cv = zeros(numel(t),1);
co = zeros(numel(t),1);

for k = 1:numel(t)
    gv = yv(k,1:l).';
    go = yo(k,:).';
    ev(k) = norm(vfcr_residual(tv(k),gv,problem,delta),2);
    eo(k) = norm(vfcr_residual(to(k),go,problem,delta),2);

    G = problem.G(t(k));   %%纵轴
    h = problem.h(t(k));
    xv = gv(1:n);
    xo = go(1:n);
    cv(k) = 0.5*xv.'*G*xv + h.'*xv;
    co(k) = 0.5*xo.'*G*xo + h.'*xo;
end

tail = t >= 5;
xv = yv(:,1:n);
xo = yo(:,1:n);
gap = vecnorm(xv(tail,:)-xo(tail,:),2,2);
dif = abs(cv(tail)-co(tail));
errv = max(ev(tail));
erro = max(eo(tail));
errx = max(gap);
errf = max(dif);
overlap = errf < 1e-3 && errv < 1e-3 ...
    && erro < 1e-3 && errx < 1e-3;

fig = figure('Color','w','Position',[120,120,800,480]);
plot(to,co,'--','Color',[0.15,0.45,0.85],'LineWidth',1.5)
hold on
plot(tv,cv,'-','Color',[0.90,0.20,0.15],'LineWidth',1.5)
xlim([0,10])
xlabel('t (s)')
ylabel('f(x(t),t)')
title('Example 2：目标函数值对比')
legend('OZNN','VFCR-ZNN','Location','northeast')
grid on
box on

savefig(fig,'example2_oznn_vfcr.fig')
exportgraphics(fig,'example2_oznn_vfcr.png','Resolution',200)

fprintf('OZNN 用时：%.6f s\n',timeo)
fprintf('VFCR-ZNN 用时：%.6f s\n',timev)
fprintf('t>=5 时 OZNN 最大残差：%.3e\n',erro)
fprintf('t>=5 时 VFCR-ZNN 最大残差：%.3e\n',errv)
fprintf('t>=5 时两者最大 x 差：%.3e\n',errx)
fprintf('t>=5 时两者最大目标函数差：%.3e\n',errf)
fprintf('后半段重合检查：%d\n',overlap)

save all
