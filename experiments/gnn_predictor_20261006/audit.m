function audit()
% 独立核对四个模型的公式和数值入口。
folder = fileparts(mfilename('fullpath'));
before = pwd;
guard = onCleanup(@() cd(before));
cd(folder);
rng(20261006,'twister');
c = struct('gamma',200,'eps',0.01,'delta',1e-4,'lambda',1e-6, ...
    'h',1e-5,'floor',1e-10,'kappa',0.05,'scale',[1,10,10], ...
    'limit',[100,5000,5000],'start',0);
rows = cell(0,1);
for trial = 1:20
    n = 3;
    m = mod(trial,2);
    w = 4*(mod(trial,3) ~= 0);
    A = randn(n);
    S = A.'*A+eye(n);
    K = randn(n);
    E = randn(n)/20;
    E = (E+E.')/2;
    G = S+K-K.';
    h = randn(n,1);
    H = randn(n,1)/5;
    P = randn(m,n);
    R = randn(m,n)/10;
    u = randn(m,1);
    U = randn(m,1)/5;
    Q = randn(w,n);
    T = randn(w,n)/10;
    v = randn(w,1);
    W = randn(w,1)/5;
    problem = struct('G',@(t) G+t*E,'h',@(t) h+t*H, ...
        'P',@(t) P+t*R,'u',@(t) u+t*U,'Q',@(t) Q+t*T, ...
        'v',@(t) v+t*W);
    for name = ["dG","dh","dP","du","dQ","dv"]
        problem.(name) = @(t) forbidden();
    end
    l = n+m+w;
    g = randn(l,1);
    for t = [0,c.h/2,c.h,0.3]
        [xi,J] = pp_parts(t,g,problem,c.delta);
        direct = original(t,g,problem,c.delta);
        D = zeros(l);
        for j = 1:l
            step = zeros(l,1);
            step(j) = 1e-6;
            D(:,j) = (pp_parts(t,g+step,problem,c.delta)- ...
                pp_parts(t,g-step,problem,c.delta))/(2e-6);
        end
        jacobian = norm(J-D,'fro')/max(1,norm(J,'fro'));
        residual = norm(xi-direct)/max(1,norm(direct));
        x = g(1:n);
        mu = g(n+m+1:end);
        omega = problem.v(t)-problem.Q(t)*x;
        sigma = sqrt(omega.^2+mu.^2+c.delta);
        analytic = [E*x+H+R.'*g(n+1:n+m)+T.'*mu; ...
            R*x-U;(1-omega./sigma).*(W-T*x)];
        if t >= c.start+c.h
            fixed = (xi-pp_parts(t-c.h,g,problem,c.delta))/c.h;
            boundary = false;
        else
            fixed = (pp_parts(t+c.h,g,problem,c.delta)-xi)/c.h;
            boundary = true;
        end
        for model = ["M3","M5","M6","M7"]
            [dg,a] = pp_system(t,g,problem,c,model);
            scale = ones(l,1);
            expected = c.lambda;
            ft = fixed;
            if model == "M3"
                ft = zeros(l,1);
            end
            if model == "M7"
                scale = [c.scale(1)*ones(n,1);c.scale(2)*ones(m,1); ...
                    c.scale(3)*ones(w,1)];
            end
            B = J*diag(scale);
            if any(model == ["M6","M7"])
                expected = max(c.floor,c.kappa*min(svd(B))^2);
            end
            matrix = B.'*B+expected*eye(l);
            dc = matrix\(B.'*xi);
            dp = -matrix\(B.'*ft);
            raw = scale.*(dp-c.gamma*dc/hypot(norm(dc),c.eps));
            speed = [norm(raw(1:n)),norm(raw(n+1:n+m)),norm(raw(n+m+1:end))];
            alpha = 1;
            if model == "M7"
                alpha = min([1,c.limit./speed]);
            end
            formula = norm(dg-alpha*raw)/max(1,norm(alpha*raw));
            predictor = norm(a.dp-scale.*dp)/max(1,norm(scale.*dp));
            correction = norm(a.dc-scale.*dc)/max(1,norm(scale.*dc));
            difference = norm(a.ft-ft)/max(1,norm(ft));
            damping = abs(a.lambda-expected)/expected;
            drift = norm(fixed-analytic)/max(1,norm(analytic));
            assert(a.boundary == (boundary && model ~= "M3"));
            assert(isequal(a.blocks,[n,m,w]));
            assert(norm(a.B-B,'fro') < 1e-12);
            assert(abs(a.sigma-min(svd(J))) < 1e-12);
            item = struct('trial',trial,'time',t,'model',model, ...
                'jacobian',jacobian,'residual',residual,'formula',formula, ...
                'predictor',predictor,'correction',correction, ...
                'difference',difference,'damping',damping,'drift',drift);
            rows{end+1,1} = item;
        end
        [a,b] = m5_system(t,g,problem,c);
        [e,f] = pp_system(t,g,problem,c,"M5");
        assert(isequal(a,e) && isequal(b,f));
        [a,b] = m6_system(t,g,problem,c);
        [e,f] = pp_system(t,g,problem,c,"M6");
        assert(isequal(a,e) && isequal(b,f));
        [a,b] = m7_system(t,g,problem,c);
        [e,f] = pp_system(t,g,problem,c,"M7");
        assert(isequal(a,e) && isequal(b,f));
    end
end
stats = struct2table(vertcat(rows{:}));
assert(max(stats.jacobian) < 1e-6);
assert(max(stats.residual) < 1e-12);
assert(max(stats.formula) < 1e-8);
assert(max(stats.predictor) < 1e-8);
assert(max(stats.correction) < 1e-8);
assert(max(stats.difference) < 1e-12);
assert(max(stats.damping) < 1e-12);
assert(max(stats.drift) < 1e-4);

problem = struct('G',@(t) 1,'h',@(t) -t,'P',@(t) zeros(0,1), ...
    'u',@(t) zeros(0,1),'Q',@(t) zeros(0,1),'v',@(t) zeros(0,1));
[zero.m3,~] = pp_system(0.3,0.3,problem,c,"M3");
[zero.m5,a] = m5_system(0.3,0.3,problem,c);
zero.expected = 1/(1+c.lambda);
assert(zero.m3 == 0 && abs(zero.m5-zero.expected) < 1e-9);
assert(a.dc == 0 && abs(a.ft+1) < 1e-10);
[~,a] = m5_system(c.h/2,c.h/2,problem,c);
assert(a.boundary && abs(a.ft+1) < 1e-10);

problem = struct('G',@(t) 1,'h',@(t) 0,'P',@(t) zeros(0,1), ...
    'u',@(t) zeros(0,1),'Q',@(t) 0,'v',@(t) 1e8);
[xi,J] = pp_parts(0,[0;1e-9],problem,c.delta);
pfb.positive = xi(end);
pfb.derivative = J(end,end);
assert(abs(xi(end)-9.995e-10) < 1e-23);
problem.v = @(t) -1e8;
[xi,J] = pp_parts(0,[0;-1e8],problem,c.delta);
pfb.negative = xi(end);
assert(all(isfinite(xi)) && all(isfinite(J),'all'));

problem = struct('G',@(t) [2,0.2;0.2,1],'h',@(t) [0.4;-0.3], ...
    'P',@(t) [1,-1],'u',@(t) 0.2,'Q',@(t) [eye(2);-eye(2)], ...
    'v',@(t) 1.2*ones(4,1));
g = [0.2;-0.1;0.3;0.1;0.2;0.3;0.4];
d = c;
d.limit = [0.01,0.02,0.03];
[dg,a] = m7_system(0.3,g,problem,d);
limit.alpha = a.alpha;
limit.speed = [norm(dg(1:2)),norm(dg(3)),norm(dg(4:7))];
limit.direction = norm(dg-a.alpha*a.raw);
limit.descent = a.xi.'*a.J*dg;
assert(a.alpha < 1 && limit.direction == 0);
assert(all(limit.speed <= d.limit*(1+1e-12)));
assert(limit.descent < 0);
assert(norm(a.ft) == 0 && norm(a.dp) == 0);

writetable(stats,fullfile(folder,'audit.csv'));
save(fullfile(folder,'audit.mat'),'stats','zero','pfb','limit','c');
fprintf('cases=%d, Jacobian=%g, residual=%g, matrix formula=%g\n', ...
    height(stats),max(stats.jacobian),max(stats.residual),max(stats.formula));
fprintf('predictor=%g, correction=%g, fixed-g difference=%g, damping=%g\n', ...
    max(stats.predictor),max(stats.correction),max(stats.difference),max(stats.damping));
fprintf('time-difference versus analytic partial=%g\n',max(stats.drift));
disp(zero);
disp(pfb);
disp(limit);
fprintf('All assertions passed; derivative fields were forbidden.\n');
end

function xi = original(t,g,problem,delta)
% 用原平方根公式独立计算残差。
G = problem.G(t);
G = (G+G.')/2;
P = problem.P(t);
Q = problem.Q(t);
n = size(G,1);
m = size(P,1);
x = g(1:n);
mu = g(n+m+1:end);
omega = problem.v(t)-Q*x;
xi = [G*x+problem.h(t)+P.'*g(n+1:n+m)+Q.'*mu; ...
    P*x-problem.u(t);omega+mu-sqrt(omega.^2+mu.^2+delta)];
end

function value = forbidden()
% 禁止读取系数导数字段。
error('Audit:Derivative','模型读取了禁止使用的导数字段。');
value = [];
end
