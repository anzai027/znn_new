function speed()
% 独立复核参考根的局部速度尖峰。
folder = fileparts(mfilename('fullpath'));
source = fullfile(fileparts(fileparts(folder)),'gnn_v2_20261002');
root = fileparts(fileparts(fileparts(folder)));
addpath(source,'-begin'); addpath(root,'-begin');
p = example2_problem(); delta = 1e-4;
t = (3.46105:1e-7:3.46115).';
values = zeros(numel(t),1);
for k = 1:numel(t), values(k) = rate(t(k),p,delta); end
[~,id] = max(values);
lo = t(max(1,id-2)); hi = t(min(numel(t),id+2));
opt = optimset('TolX',1e-13,'Display','off');
peak = fminbnd(@(t) -rate(t,p,delta),lo,hi,opt);
[value,g,dy,bdy,xi,J,ft,slack,mu,H] = rate(peak,p,delta);
h = [1e-5;1e-6;1e-7;1e-8;1e-9];
diffs = zeros(numel(h),7); errors = zeros(numel(h),1);
for k = 1:numel(h)
    ref = reference([peak-h(k);peak+h(k)],p,delta);
    diffs(k,:) = (ref.g(2,:)-ref.g(1,:))/(2*h(k));
    errors(k) = norm(diffs(k,:).'-bdy)/norm(bdy);
end
fd = table(h,vecnorm(diffs,2,2),errors,'VariableNames',{'h','speed','relative'});
writetable(fd,fullfile(folder,'speed_diff.csv'));
block = [norm(dy(1:2)),abs(dy(3)),norm(dy(4:7))];
summary = table(peak,value,norm(bdy),norm(dy-bdy)/norm(bdy),norm(xi), ...
    norm(parts(peak,g,p,delta)),min(svd(J)),block(1),block(2),block(3), ...
    'VariableNames',{'time','full','barrier','relative','stable','original','sigma','x','mu1','mu2'});
writetable(summary,fullfile(folder,'speed_summary.csv'));
details = table((1:4).',slack,mu,dy(4:7), ...
    'VariableNames',{'constraint','slack','mu','velocity'});
writetable(details,fullfile(folder,'speed_dual.csv'));
ref = reference(t,p,delta);
save(fullfile(folder,'speed.mat'),'peak','value','g','dy','bdy','xi','J','ft', ...
    'slack','mu','H','h','diffs','errors','summary','details','t','values','ref');
fig = figure('Visible','off','Color','w','Position',[80,80,1200,760]);
subplot(2,2,1); plot(t,values,'LineWidth',1.4); xline(peak,':k');
xlabel('t'); ylabel('Root speed'); grid on; title('Independent barrier derivative');
subplot(2,2,2); plot(t,ref.g(:,3:7),'LineWidth',1.2);
xlabel('t'); ylabel('Dual variables'); grid on; legend('mu1','mu2_1','mu2_2','mu2_3','mu2_4','Location','best');
subplot(2,2,3); plot(t,ref.xd,'LineWidth',1.2);
xlabel('t'); ylabel('x'); grid on; legend('x1','x2','Location','best');
subplot(2,2,4); loglog(h,errors,'o-','LineWidth',1.2);
xlabel('Central difference step'); ylabel('Relative error to barrier derivative'); grid on;
exportgraphics(fig,fullfile(folder,'speed.png'),'Resolution',150);
savefig(fig,fullfile(folder,'speed.fig')); close(fig);
disp(summary); disp(fd); disp(details);
fprintf('P = [%.12g %.12g], u = %.12g\n',p.P(peak),p.u(peak));
fprintf('x = [%.12g %.12g], mu1 = %.12g\n',g(1:3));
fprintf('Full equation relative residual %.12g\n',norm(J*dy+ft)/norm(ft));
end

function [value,g,dy,bdy,xi,J,ft,slack,mu,H] = rate(t,p,delta)
% 障碍方程与完整 PFB 方程分别求导。
ref = reference(t,p,delta); g = ref.g.';
G = p.G(t); G = (G+G.')/2;
dG = p.dG(t); dG = (dG+dG.')/2;
P = p.P(t); Q = p.Q(t); x = g(1:2); mu1 = g(3); mu = g(4:7);
slack = p.v(t)-Q*x;
sigma = sqrt(slack.^2+mu.^2+delta);
a = (mu.^2+delta)./(sigma.*(sigma+slack));
b = (slack.^2+delta)./(sigma.*(sigma+mu));
phi = (2*slack.*mu-delta)./(slack+mu+sigma);
xi = [G*x+p.h(t)+P.'*mu1+Q.'*mu; P*x-p.u(t); phi];
J = [G,P.',Q.';P,0,zeros(1,4);-diag(a)*Q,zeros(4,1),diag(b)];
omega = p.dv(t)-p.dQ(t)*x;
first = dG*x+p.dh(t)+p.dP(t).'*mu1+p.dQ(t).'*mu;
second = p.dP(t)*x-p.du(t);
ft = [first;second;a.*omega];
dy = -J\ft;
D = diag(delta./(2*slack.^2));
H = G+Q.'*D*Q;
B = [H,P.';P,0];
part = -B\[first-Q.'*D*omega;second];
bdy = [part;-D*(omega-Q*part(1:2))];
value = norm(bdy);
end
