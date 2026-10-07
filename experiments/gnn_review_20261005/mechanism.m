function mechanism()
% 用独立诊断区分跟踪误差和数值误差。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
v2 = fullfile(root,'experiments','gnn_v2_20261002');
prior = fullfile(root,'experiments','gnn_vfcr_20261002');
before = pwd; guard = onCleanup(@() cd(before)); cd(root);
addpath(v2,'-begin'); addpath(root,'-begin');
out = fullfile(folder,'mechanism'); if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'run.log')); clean = onCleanup(@() diary('off'));
c = config(); p = example2_problem();
t = (0:0.0005:10).'; ref = reference(t,p,c.delta);
r = velocities(ref,p,c.delta);
[~,id] = max(r.speed);
fine = (max(0,t(id)-0.02):0.00001:min(10,t(id)+0.02)).';
rf = reference(fine,p,c.delta); vf = velocities(rf,p,c.delta);
coarse = reference(c.t,p,c.delta); vc = velocities(coarse,p,c.delta);
rows = repmat(struct(),0,1);
for k = 1:3
    refs = {coarse,ref,rf}; speeds = {vc,r,vf}; steps = [0.005,0.0005,0.00001];
    a = refs{k}; b = speeds{k};
    diff = (a.g(3:end,:)-a.g(1:end-2,:))/(2*steps(k));
    [peak,id] = max(b.speed);
    row = struct('step',steps(k),'implicit',peak,'time',a.t(id), ...
        'difference',max(vecnorm(diff,2,2)),'x',max(b.x), ...
        'mu1',max(b.mu1),'mu2',max(b.mu2),'sigma',min(b.sigma));
    if isempty(rows), rows = row; else, rows(end+1) = row; end
end
speed = struct2table(rows); writetable(speed,fullfile(out,'speed.csv'));
disp(speed);
[~,id] = max(vf.speed); time = rf.t(id); g = rf.g(id,:).';
h = 1e-6; check = reference([time-h;time+h],p,c.delta);
diff = (check.g(2,:)-check.g(1,:)).'/(2*h);
validation = norm(diff-vf.velocity(id,:).')/max(1,norm(diff));
fprintf('峰值速度的独立差分相对误差 %.6g\n',validation);
lambda = [1e-6,1e-4,1e-2]; damping = table();
for k = 1:numel(lambda)
    values = [r.sigma;vf.sigma];
    ratio = values.^2./(values.^2+lambda(k));
    row = table(lambda(k),min(ratio),median(ratio), ...
        'VariableNames',{'lambda','minimum','median'});
    damping = [damping;row];
end
writetable(damping,fullfile(out,'damping.csv')); disp(damping);
data = load(fullfile(root,'experiments','gnn_vfcr_smooth_20261002','data','run_010.mat'),'r');
list = jobs(c,prior);
base = list(find([list.group] == "main" & [list.model] == "M3" & [list.noise] == "zero",1));
q = p; time = 6.31;
names = {'G','h','P','u','Q','v'};
for k = 1:numel(names)
    name = names{k}; value = p.(name)(time); q.(name) = @(t) value;
    q.(['d',name]) = @(t) zeros(size(value));
end
c.t = (0:0.005:10).';
frozen = reference(c.t,q,c.delta);
point = find(abs(data.r.t-time) < 1e-10,1);
base.g0 = data.r.y(point,1:7).';
models = ["M2","M3"]; results = cell(2,1);
for k = 1:2
    j = base; j.model = models(k); j.id = k; j.group = "frozen"; j.reuse = "";
    bank = noises("zero",1,c); a = simulate(j,q,bank,c);
    a.metrics = metrics(a,frozen,q,c); results{k} = a;
    save(fullfile(out,sprintf('frozen_%d.mat',k)),'a','-v7.3');
    fprintf('冻结题目 %s %s residual=%.6g full=%.6g calls=%d\n', ...
        j.model,a.status,a.metrics.last,a.metrics.full,a.calls);
end
stats = collect(results); writetable(stats,fullfile(out,'frozen.csv'));
variants = ["slow","normal","fast","gain1000","gain5000"];
rates = [0.1,1,2,1,1]; gains = [200,200,200,1000,5000];
dynamic = cell(numel(variants),1);
base = list(find([list.group] == "main" & [list.model] == "M3" & [list.noise] == "zero" & [list.eps] > 0,1));
for k = 1:numel(variants)
    rate = rates(k); ci = config(); ci.t = ci.t/rate; ci.tail = ci.tail/rate;
    ci.step = ci.step/rate;
    q = rescale(p,rate); j = base; j.id = k; j.group = variants(k);
    j.gamma = gains(k); j.rates = gains(k); j.step = ci.step; j.reuse = "";
    target = reference(ci.t,q,ci.delta); bank = noises("zero",1,ci);
    a = simulate(j,q,bank,ci); a.metrics = metrics(a,target,q,ci);
    dynamic{k} = a;
    save(fullfile(out,variants(k)+".mat"),'a','target','rate','-v7.3');
    fprintf('时变诊断 %s %s RMSE=%.6g full=%.6g residual=%.6g calls=%d\n', ...
        variants(k),a.status,a.metrics.rmsed,a.metrics.full,a.metrics.residual,a.calls);
end
variation = collect(dynamic); writetable(variation,fullfile(out,'variation.csv'));
geometry = table();
for model = ["M3","M4"]
    id = find([list.group] == "main" & [list.model] == model & [list.noise] == "zero" & [list.eps] > 0,1);
    a = load(fullfile(v2,'data',sprintf('run_%03d.mat',list(id).id)),'r'); a = a.r;
    for k = 1:numel(a.t)
        t = a.t(k); g = a.y(k,1:7).';
        [xi,J] = parts(t,g,p,c.delta);
        [dg,d] = pg_system(t,g,7,p,a.job.rates,a.job.lambda,a.job.eps,c.delta,model == "M4");
        e = g-coarse.g(k,:).';
        value = dot(J.'*xi,dg);
        angle = dot(e,-dg)/max(realmin,norm(e)*norm(dg));
        before = dot(e,d)/max(realmin,norm(e)*norm(d));
        row = table(model,t,value,angle,before,'VariableNames',{'model','t','slope','actual','direction'});
        geometry = [geometry;row];
    end
end
writetable(geometry,fullfile(out,'geometry.csv'));
fprintf('M4 采样反馈使残差上升的点数 %d\n',sum(geometry.model == "M4" & geometry.slope > 1e-10));
save(fullfile(out,'analysis.mat'),'ref','r','rf','vf','coarse','vc','speed', ...
    'damping','validation','results','stats','dynamic','variation','geometry','-v7.3');
fig = figure('Visible','off','Position',[80,80,1100,800]);
tiledlayout(3,1);
nexttile; plot(ref.t,r.speed,'LineWidth',1.1); yline(100,'--'); yline(200,'--');
ylabel('Root speed'); legend('root','100','200','Location','best'); grid on;
nexttile; semilogy(ref.t,r.sigma); ylabel('Smallest singular value'); grid on;
nexttile; hold on;
for k = 1:2, semilogy(results{k}.t,results{k}.metrics.e); end
set(gca,'YScale','log'); xlabel('Time (s)'); ylabel('Frozen residual');
legend('M2','M3','Location','best'); grid on;
exportgraphics(fig,fullfile(out,'diagnosis.png'),'Resolution',160);
savefig(fig,fullfile(out,'diagnosis.fig')); close(fig);
end

function q = rescale(p,rate)
% 只改变题目运动速度。
names = {'G','h','P','u','Q','v'}; q = p;
for k = 1:numel(names)
    name = names{k}; fun = p.(name); q.(name) = @(t) fun(rate*t);
    fun = p.(['d',name]); q.(['d',name]) = @(t) rate*fun(rate*t);
end
end

function r = velocities(ref,p,delta)
% 隐式微分只用于参考轨迹诊断。
count = numel(ref.t); r.speed = zeros(count,1); r.sigma = zeros(count,1);
r.x = zeros(count,1); r.mu1 = zeros(count,1); r.mu2 = zeros(count,1);
r.velocity = zeros(count,7);
for k = 1:count
    t = ref.t(k); g = ref.g(k,:).';
    [~,J] = parts(t,g,p,delta);
    x = g(1:2); mu1 = g(3); mu2 = g(4:7);
    omega = p.v(t)-p.Q(t)*x;
    sigma = sqrt(omega.^2+mu2.^2+delta);
    dG = p.dG(t); dG = (dG+dG.')/2;
    ft = [dG*x+p.dh(t)+p.dP(t).'*mu1+p.dQ(t).'*mu2; ...
        p.dP(t)*x-p.du(t); (1-omega./sigma).*(p.dv(t)-p.dQ(t)*x)];
    dg = -J\ft; r.velocity(k,:) = dg.';
    r.speed(k) = norm(dg); r.x(k) = norm(dg(1:2));
    r.mu1(k) = norm(dg(3)); r.mu2(k) = norm(dg(4:7));
    r.sigma(k) = min(svd(J));
end
end
