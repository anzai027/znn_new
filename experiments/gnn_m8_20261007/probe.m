function probe()
% 在参考根处分开检查预测阻尼和时间差分。
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
old = fullfile(root,'experiments','gnn_predictor_20261006');
before = pwd; guard = onCleanup(@() cd(before));
cd(root); addpath(root,'-begin'); addpath(old,'-begin'); addpath(folder,'-begin');
for name = {"pn_config","pn_system"}
    assert(strcmpi(which(name{1}),fullfile(folder,name{1}+".m")));
end
for name = {"pp_parts","pp_reference","pp_config"}
    assert(strcmpi(which(name{1}),fullfile(old,name{1}+".m")));
end
assert(strcmpi(which('example2_problem'),fullfile(root,'example2_problem.m')));
out = fullfile(folder,'diagnostics');
if ~isfolder(out), mkdir(out); end
assert(~isfile(fullfile(out,'prediction.mat')),'已有诊断结果，请先保留原输出。');
diary(fullfile(out,'prediction.log')); clean = onCleanup(@() diary('off'));
files = [fullfile(folder,{'pn_system.m','pn_config.m','probe.m'}), ...
    fullfile(old,{'pp_parts.m','pp_reference.m','pp_config.m'}), ...
    {fullfile(root,'example2_problem.m')}];
meta.files = files; meta.hash = cellfun(@digest,files,'UniformOutput',false);
c = pn_config("M8"); p = example2_problem();
t = [0;1;3.46;3.4611091472136;6.31;9.7442944543932;10];
h = [1e-5,1e-6]; damping = [1e-10,1e-9,1e-8];
ref = pp_reference(t,p,c.delta);
names = {'G','h','P','u','Q','v'}; q = struct();
for k = 1:numel(names), q.(names{k}) = p.(names{k}); end
rows = repmat(struct(),0,1); results = cell(numel(t),numel(h),4);
for i = 1:numel(t)
    time = t(i); g = ref.g(i,:).'; velocity = ref.velocity(i,:).';
    [ft,blocks] = derivative(time,g,p,c.delta); n = blocks(1); m = blocks(2); l = numel(g);
    for j = 1:numel(h)
        for k = 1:4
            model = "M8";
            if k == 4, model = "M8A"; end
            ci = pn_config(model); ci.h = h(j);
            if k < 4, ci.predict = damping(k); end
            [dg,a] = pn_system(time,g,q,ci,model);
            if model == "M8A"
                exact = -(a.J\ft); retention = 1;
            else
                exact = [a.J;sqrt(a.lp)*eye(l)]\[-ft;zeros(l,1)];
                retention = a.sigma^2/(a.sigma^2+a.lp);
            end
            loss = norm(exact-velocity); difference = norm(a.dp-exact); total = norm(a.dp-velocity);
            assert(total <= loss+difference+1e-10*max(1,total));
            assert(total >= abs(loss-difference)-1e-10*max(1,total));
            speed = [norm(velocity(1:n)),norm(velocity(n+1:n+m)),norm(velocity(n+m+1:end))];
            prediction = [norm(a.dp(1:n)),norm(a.dp(n+1:n+m)),norm(a.dp(n+m+1:end))];
            sv = svd(a.J);
            row = struct('time',time,'h',ci.h,'model',model,'correct',string(ci.correct), ...
                'lc',a.lc,'lp',a.lp,'sigma',a.sigma,'rcond',a.condition, ...
                'condition',max(sv)/min(sv),'retention',retention,'rootres',ref.error(i), ...
                'identity',norm(a.J*velocity+ft),'ftrel',norm(a.ft-ft)/max(norm(ft),realmin), ...
                'damping',loss,'difference',difference,'total',total, ...
                'relative',total/max(norm(velocity),realmin),'balance',a.balance, ...
                'mismatch',norm(a.J*a.dp+ft),'reference',norm(velocity),'prediction',norm(a.dp), ...
                'refx',speed(1),'refmu1',speed(2),'refmu2',speed(3), ...
                'predx',prediction(1),'predmu1',prediction(2),'predmu2',prediction(3), ...
                'feedback',norm(a.feedback),'executed',norm(dg-velocity),'forward',a.boundary);
            if isempty(rows), rows = row; else, rows(end+1,1) = row; end
            results{i,j,k} = struct('a',a,'exact',exact,'ft',ft,'velocity',velocity,'c',ci,'row',row);
        end
    end
end
stats = struct2table(rows); assert(height(stats) == 56);
for k = 1:numel(files), assert(strcmp(digest(files{k}),meta.hash{k}),'模型或诊断源码改变。'); end
meta.matlab = string(version); meta.date = string(datetime('now','TimeZone','Asia/Shanghai'));
meta.scope = "解析导数仅用于独立根速度诊断，模型只得到六个题目系数字段。";
save(fullfile(out,'prediction.mat'),'c','t','h','damping','ref','stats','results','meta','-v7.3');
writetable(stats,fullfile(out,'prediction.csv'),'Encoding','UTF-8');
fprintf('完成 %d 个根处诊断点，参考最大残差 %.12g。\n',height(stats),max(ref.error));
take = abs(stats.time-3.4611091472136) < 1e-12;
disp(stats(take,{'model','h','lc','lp','retention','damping','difference','total','reference','prediction'}));
end

function [ft,blocks] = derivative(t,g,p,delta)
% 在固定状态下独立计算残差时间偏导。
G = p.G(t); dG = p.dG(t); dG = (dG+dG.')/2;
P = p.P(t); Q = p.Q(t);
n = size(G,1); m = size(P,1); w = size(Q,1); blocks = [n,m,w];
x = g(1:n); mu1 = g(n+1:n+m); mu = g(n+m+1:end);
omega = p.v(t)-Q*x;
if isscalar(delta), delta = delta*ones(w,1); else, delta = delta(:); end
sigma = hypot(hypot(omega,mu),sqrt(delta));
alpha = (mu.^2+delta)./(sigma.*(sigma+omega));
ft = [dG*x+p.dh(t)+p.dP(t).'*mu1+p.dQ(t).'*mu; ...
    p.dP(t)*x-p.du(t);alpha.*(p.dv(t)-p.dQ(t)*x)];
end

function hash = digest(file)
% 保存参与计算的源码校验值。
fid = fopen(file,'rb'); assert(fid >= 0); clean = onCleanup(@() fclose(fid));
bytes = fread(fid,Inf,'*uint8'); md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes); value = typecast(md.digest(),'uint8');
hash = lower(reshape(dec2hex(value,2).',1,[]));
end
