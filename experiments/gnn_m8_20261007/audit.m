function audit()
% 独立核对分离阻尼和直接预测的公式。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
prior = fullfile(root,'experiments','gnn_predictor_20261006');
addpath(root,'-begin');
addpath(prior,'-begin');
addpath(folder,'-begin');
assert(strcmpi(which('pn_system'),fullfile(folder,'pn_system.m')));
assert(strcmpi(which('audit'),fullfile(folder,'audit.m')));
rng(20261007,'twister');
base = pn_config("M8");
deltas = [1e-8,1e-6,1e-4,1e-2,0.1];
rows = cell(0,1);
for trial = 1:40
    n = 2+mod(trial,2);
    m = mod(floor((trial-1)/2),2);
    w = 4*mod(floor((trial-1)/4),2);
    A = randn(n);
    S = A.'*A+eye(n);
    K = randn(n);
    G = S+K-K.';
    E = randn(n)/20;
    h = randn(n,1);
    H = randn(n,1)/5;
    P = randn(m,n);
    R = randn(m,n)/10;
    u = randn(m,1);
    U = randn(m,1)/5;
    Q = randn(w,n);
    T = randn(w,n)/10;
    g = [randn(n+m,1);0.2+0.3*rand(w,1)];
    v = Q*g(1:n)+1+rand(w,1);
    W = randn(w,1)/5;
    p = struct('G',@(t) G+t*E,'h',@(t) h+t*H,'P',@(t) P+t*R, ...
        'u',@(t) u+t*U,'Q',@(t) Q+t*T,'v',@(t) v+t*W, ...
        'dG',@(t) E,'dh',@(t) H,'dP',@(t) R,'du',@(t) U, ...
        'dQ',@(t) T,'dv',@(t) W);
    q = masked(p);
    c = base;
    c.delta = deltas(1+floor((trial-1)/8));
    if w > 0 && mod(trial,2) == 0
        c.delta = c.delta*(1+(1:w).'/w);
    end
    l = numel(g);
    for t = [0,c.h/2,c.h,0.3]
        [xi,J,trueft] = m8_oracle(t,g,p,c.delta);
        [actual,actualJ,blocks] = pp_parts(t,g,q,c.delta);
        D = zeros(l);
        for j = 1:l
            step = zeros(l,1);
            step(j) = 1e-6;
            D(:,j) = (m8_oracle(t,g+step,p,c.delta)- ...
                m8_oracle(t,g-step,p,c.delta))/(2e-6);
        end
        if t < c.start+c.h
            ft = (m8_oracle(t+c.h,g,p,c.delta)-xi)/c.h;
        else
            ft = (xi-m8_oracle(t-c.h,g,p,c.delta))/c.h;
        end
        [L,F,V] = svd(J);
        values = diag(F);
        for model = ["M8","M8A"]
            ci = c;
            ci.correct = "fixed";
            lc = ci.lambda;
            lp = ci.predict;
            if model == "M8A"
                ci.correct = "adaptive";
                lc = max(ci.floor,ci.kappa*min(values)^2);
                lp = 0;
            end
            [dg,a] = pn_system(t,g,q,ci,model);
            dc = (J.'*J+lc*eye(l))\(J.'*xi);
            dcs = V*((values./(values.^2+lc)).*(L.'*xi));
            if model == "M8"
                dp = -(J.'*J+lp*eye(l))\(J.'*ft);
                dps = -V*((values./(values.^2+lp)).*(L.'*ft));
            else
                dp = -J\ft;
                dps = -V*((L.'*ft)./values);
            end
            expected = dp-ci.gamma*dc/hypot(norm(dc),ci.eps);
            item = struct('trial',trial,'time',t,'model',model, ...
                'jacobian',norm(J-D,'fro')/max(1,norm(J,'fro')), ...
                'residual',norm(actual-xi)/max(1,norm(xi)), ...
                'derivative',norm(actualJ-J,'fro')/max(1,norm(J,'fro')), ...
                'correction',norm(a.dc-dc)/max(1,norm(dc)), ...
                'prediction',norm(a.dp-dp)/max(1,norm(dp)), ...
                'svdc',norm(a.dc-dcs)/max(1,norm(dcs)), ...
                'svdp',norm(a.dp-dps)/max(1,norm(dps)), ...
                'formula',norm(dg-expected)/max(1,norm(expected)), ...
                'fixed',norm(a.ft-ft)/max(1,norm(ft)), ...
                'drift',norm(ft-trueft)/max(1,norm(trueft)), ...
                'balance',norm(J*a.dp+a.ft)/max(1,norm(a.ft)));
            assert(abs(a.lc-lc) <= 1e-12*max(lc,realmin) && a.lp == lp);
            assert(a.boundary == (t < c.start+c.h));
            assert(isequal(blocks,[n,m,w]) && isequal(a.blocks,blocks));
            assert(abs(a.balance-norm(J*a.dp+a.ft)) < 1e-12);
            assert(abs(a.condition-rcond(J)) < 1e-12);
            rows{end+1,1} = item;
        end
        [first,metadata] = pn_system(t,g,q,c,"M8");
        pn_system(0.3,g,q,c,"M8");
        pn_system(0,g,q,c,"M8");
        [last,again] = pn_system(t,g,q,c,"M8");
        assert(isequal(first,last) && isequal(metadata,again));
    end
end
stats = struct2table(vertcat(rows{:}));
assert(max(stats.jacobian) < 1e-6 && max(stats.derivative) < 1e-12);
assert(max(stats.residual) < 1e-12 && max(stats.fixed) < 1e-12);
assert(max(stats.correction) < 1e-8 && max(stats.prediction) < 1e-8);
assert(max(stats.svdc) < 1e-8 && max(stats.svdp) < 1e-8);
assert(max(stats.formula) < 1e-8 && max(stats.drift) < 1e-4);

p = struct('G',@(t) 1,'h',@(t) -t,'P',@(t) zeros(0,1), ...
    'u',@(t) zeros(0,1),'Q',@(t) zeros(0,1),'v',@(t) zeros(0,1));
q = masked(p);
c = base;
[zero.m8,a] = pn_system(0.3,0.3,q,c,"M8");
ca = pn_config("M8A");
[zero.direct,b] = pn_system(0.3,0.3,q,ca,"M8A");
zero.expected = 1/(1+c.predict);
assert(abs(zero.m8-zero.expected) < 1e-9 && abs(zero.direct-1) < 1e-9);
assert(a.dc == 0 && b.dc == 0 && abs(a.ft+1) < 1e-9);
assert(b.lp == 0 && abs(b.lc-0.01) < 1e-12);
[~,a] = pn_system(0,0,q,c,"M8");
assert(a.boundary && abs(a.ft+1) < 1e-12);
[~,a] = pn_system(0.3,0.5,q,c,"M8");
ci = c;
ci.lambda = 10*c.lambda;
[~,b] = pn_system(0.3,0.5,q,ci,"M8");
split.predictor = norm(a.dp-b.dp);
split.correction = norm(a.dc-b.dc);
assert(split.predictor == 0 && split.correction > 1e-10);
ci = c;
ci.predict = 100*c.predict;
[~,b] = pn_system(0.3,0.5,q,ci,"M8");
split.fixed = norm(a.dc-b.dc);
split.changed = norm(a.dp-b.dp);
assert(split.fixed == 0 && split.changed > 1e-10);
ci = ca;
ci.correct = "fixed";
[~,b] = pn_system(0.3,0.5,q,ci,"M8A");
assert(norm(a.dc-b.dc) == 0 && b.lc == c.lambda && b.lp == 0);

p = struct('G',@(t) 1,'h',@(t) 0,'P',@(t) zeros(0,1), ...
    'u',@(t) zeros(0,1),'Q',@(t) 1,'v',@(t) 1e8);
[xi,J] = pp_parts(0,[0;1e-9],p,1e-4);
pfb.positive = xi(end);
pfb.alpha = -J(end,1);
assert(abs(xi(end)-9.995e-10) < 1e-23);
assert(abs(pfb.alpha-5e-21) < 1e-32);
p.v = @(t) -1e8;
[xi,J] = pp_parts(0,[0;-1e8],p,1e-4);
pfb.negative = xi(end);
assert(all(isfinite(xi)) && all(isfinite(J),'all'));

p = example2_problem();
q = masked(p);
times = unique([0,0.17,1,2.3,5,6.31,10, ...
    base.peak(1)+[-1e-3,-1e-4,-1e-5,0,1e-5,1e-4,1e-3], ...
    base.peak(2)+[-1e-3,-1e-4,-1e-5,0,1e-5,1e-4,1e-3]]).';
ref = pp_reference(times,p,base.delta);
rows = cell(0,1);
for k = 1:numel(times)
    t = times(k);
    g = ref.g(k,:).';
    velocity = ref.velocity(k,:).';
    [xi,J,trueft] = m8_oracle(t,g,p,base.delta);
    oracle = -J\trueft;
    for step = [1e-4,1e-5,1e-6]
        for model = ["M8","M8A"]
            c = pn_config(model);
            c.h = step;
            [~,a] = pn_system(t,g,q,c,model);
            if model == "M8"
                [L,F,V] = svd(J);
                values = diag(F);
                analytic = -V*((values./(values.^2+a.lp)).*(L.'*trueft));
            else
                analytic = oracle;
            end
            item = struct('time',t,'step',step,'model',model, ...
                'sigma',min(svd(J)),'lc',a.lc,'lp',a.lp, ...
                'root',norm(xi),'speed',norm(velocity),'predictor',norm(a.dp), ...
                'identity',norm(J*velocity+trueft)/max(1,norm(trueft)), ...
                'oracle',norm(oracle-velocity),'damping',norm(analytic-oracle), ...
                'difference',norm(a.dp-analytic),'total',norm(a.dp-velocity), ...
                'balance',a.balance);
            rows{end+1,1} = item;
        end
    end
end
peaks = struct2table(vertcat(rows{:}));
assert(max(peaks.identity) < 1e-9 && max(peaks.root) < 1e-9);
writetable(stats,fullfile(folder,'audit.csv'));
writetable(peaks,fullfile(folder,'prediction.csv'));
save(fullfile(folder,'audit.mat'),'stats','peaks','zero','split','pfb','base');
fprintf('Formula cases=%d, Jacobian=%g, independent derivative=%g\n', ...
    height(stats),max(stats.jacobian),max(stats.derivative));
fprintf('Correction=%g, predictor=%g, combined formula=%g\n', ...
    max(stats.correction),max(stats.prediction),max(stats.formula));
fprintf('SVD correction=%g, SVD predictor=%g, fixed-g difference=%g\n', ...
    max(stats.svdc),max(stats.svdp),max(stats.fixed));
fprintf('Difference versus analytic partial=%g; direct normalized balance=%g\n', ...
    max(stats.drift),max(stats.balance(stats.model == "M8A")));
disp(zero); disp(split); disp(pfb);
disp(peaks(peaks.step == 1e-5 & ismember(peaks.time,base.peak),:));
fprintf('All assertions passed; no ODE was run; derivative fields were forbidden.\n');
end

function q = masked(p)
% 保留六个系数并禁止读取解析导数。
names = {'G','h','P','u','Q','v'};
q = struct();
for k = 1:numel(names), q.(names{k}) = p.(names{k}); end
for name = ["dG","dh","dP","du","dQ","dv"]
    q.(name) = @(t) forbidden();
end
end

function value = forbidden()
% 阻止模型调用解析导数。
error('Audit:Derivative','模型读取了解析系数导数。');
value = [];
end
