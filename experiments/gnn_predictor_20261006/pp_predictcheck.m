function pp_predictcheck()
% 在独立参考根处分拆预测误差。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
before = pwd; guard = onCleanup(@() cd(before));
cd(root); addpath(root,'-begin'); addpath(folder,'-begin');
for name = {"pp_config","pp_system","pp_reference","pp_parts"}
    assert(strcmpi(which(name{1}),fullfile(folder,name{1}+".m")));
end
assert(strcmpi(which('example2_problem'),fullfile(root,'example2_problem.m')));
out = fullfile(folder,'diagnostics');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'prediction.log')); clean = onCleanup(@() diary('off'));
c = pp_config(); p = example2_problem();
t = [0;1;3.46;3.4611091472136;6.31;9.7442944543932;10];
h = [1e-4,1e-5,1e-6]; models = ["M5","M6","M7"];
ref = pp_reference(t,p,c.delta);
names = {'G','h','P','u','Q','v'}; q = struct();
for k = 1:numel(names), q.(names{k}) = p.(names{k}); end
rows = repmat(struct(),0,1); results = cell(numel(t),numel(h),numel(models));
for i = 1:numel(t)
    time = t(i); g = ref.g(i,:).'; velocity = ref.velocity(i,:).';
    [ft,blocks] = derivative(time,g,p,c.delta);
    n = blocks(1); m = blocks(2); l = numel(g);
    for j = 1:numel(h)
        ci = c; ci.h = h(j);
        for k = 1:numel(models)
            model = models(k);
            [dg,a] = pp_system(time,g,q,ci,model);
            scale = ones(l,1);
            if model == "M7"
                scale = [c.scale(1)*ones(n,1);c.scale(2)*ones(m,1); ...
                    c.scale(3)*ones(blocks(3),1)];
            end
            D = diag(scale); B = a.J*D;
            exact = D*([B;sqrt(a.lambda)*eye(l)]\[-ft;zeros(l,1)]);
            sigma = min(svd(a.J)); sigmab = min(svd(B));
            damping = norm(exact-velocity); difference = norm(a.dp-exact);
            total = norm(a.dp-velocity);
            assert(total <= damping+difference+1e-10*max(1,total));
            assert(total >= abs(damping-difference)-1e-10*max(1,total));
            assert(norm(a.B-B,'fro') <= 1e-12*max(1,norm(B,'fro')));
            speed = [norm(velocity(1:n)),norm(velocity(n+1:n+m)), ...
                norm(velocity(n+m+1:end))];
            prediction = [norm(a.dp(1:n)),norm(a.dp(n+1:n+m)), ...
                norm(a.dp(n+m+1:end))];
            row = struct('time',time,'h',ci.h,'model',model,'sigma',sigma, ...
                'sigmab',sigmab,'lambda',a.lambda, ...
                'retention',sigmab^2/(sigmab^2+a.lambda), ...
                'rootres',ref.error(i),'identity',norm(a.J*velocity+ft), ...
                'ftrel',norm(a.ft-ft)/max(norm(ft),realmin), ...
                'damping',damping,'difference',difference,'total',total, ...
                'relative',total/max(norm(velocity),realmin), ...
                'reference',norm(velocity),'prediction',norm(a.dp), ...
                'refx',speed(1),'refmu1',speed(2),'refmu2',speed(3), ...
                'predx',prediction(1),'predmu1',prediction(2),'predmu2',prediction(3), ...
                'alpha',a.alpha,'executed',norm(dg-velocity),'forward',a.boundary);
            if isempty(rows), rows = row; else, rows(end+1,1) = row; end
            results{i,j,k} = struct('a',a,'exact',exact,'ft',ft,'velocity',velocity,'row',row);
        end
    end
end
stats = struct2table(rows);
writetable(stats,fullfile(out,'prediction.csv'),'Encoding','UTF-8');
meta.matlab = string(version); meta.date = string(datetime('now','TimeZone','Asia/Shanghai'));
meta.scope = "解析系数导数只用于独立诊断，模型仅接收六个原始系数函数。";
save(fullfile(out,'prediction.mat'),'c','t','h','models','ref','stats','results','meta','-v7.3');
fprintf('完成 %d 个根处预测诊断点，参考最大残差 %.12g。\n',height(stats),max(ref.error));
take = abs(stats.time-3.4611091472136) < 1e-12;
disp(stats(take,{'model','h','lambda','retention','ftrel','damping','difference','total','reference','prediction','alpha'}));
end

function [ft,blocks] = derivative(t,g,p,delta)
% 固定当前状态后独立计算残差的时间偏导。
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
