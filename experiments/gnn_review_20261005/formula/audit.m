function audit()
% 独立检查残差、雅可比和状态导数。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(folder)));
addpath(root,'-begin');
addpath(fullfile(root,'experiments','gnn_compare_20261002'),'-end');
assert(strcmpi(which('gnn_system'),fullfile(root,'gnn_system.m')));
assert(strcmpi(which('ftcgnn_system'),fullfile(root,'experiments','gnn_compare_20261002','ftcgnn_system.m')));
rng(20261005,'twister');
rows = cell(0,1);
for n = [2,3]
    for m = [0,1]
        for w = [0,4]
            B = randn(n);
            S = B.'*B+eye(n);
            K = randn(n);
            G = S+K-K.';
            P = randn(m,n);
            Q = randn(w,n);
            h = randn(n,1);
            u = randn(m,1);
            v = 2+rand(w,1);
            problem = struct('G',@(t) G,'h',@(t) h,'P',@(t) P, ...
                'u',@(t) u,'Q',@(t) Q,'v',@(t) v);
            plain = problem;
            plain.G = @(t) S;
            l = n+m+w;
            for delta = [1e-8,1e-4,1e-2]
                for trial = 1:5
                    g = randn(l,1);
                    if w && mod(trial,2) == 0
                        d = delta*(1+(1:w).'/w);
                    else
                        d = delta;
                    end
                    [xi,J] = parts(g,problem,d,false,false);
                    D = zeros(l);
                    for j = 1:l
                        step = zeros(l,1);
                        step(j) = 1e-6;
                        D(:,j) = (parts(g+step,problem,d,false,false)- ...
                            parts(g-step,problem,d,false,false))/(2e-6);
                    end
                    for p = [0.5,1,1.5]
                        s = J.'*xi;
                        expected = -37*s/norm(s)^p;
                        actual = gnn_system(0,g,l,problem,37,p,d);
                        symmetric = gnn_system(0,g,l,plain,37,p,d);
                        item = struct('n',n,'m',m,'w',w,'delta',delta, ...
                            'vector',~isscalar(d),'trial',trial,'p',p, ...
                            'jacobian',norm(J-D,'fro')/max(1,norm(J,'fro')), ...
                            'formula',norm(actual-expected)/max(1,norm(expected)), ...
                            'symmetry',norm(actual-symmetric)/max(1,norm(actual)));
                        rows{end+1,1} = item;
                    end
                end
            end
        end
    end
end
gnn = struct2table(vertcat(rows{:}));
assert(max(gnn.jacobian) < 1e-6);
assert(max(gnn.formula) < 1e-12);
assert(max(gnn.symmetry) < 1e-12);
writetable(gnn,fullfile(folder,'gnn.csv'));

rows = cell(0,1);
for trial = 1:30
    A = randn(4,3);
    b = A*randn(3,1);
    x = randn(3,1);
    for p = [0.5,1,1.5]
        for port = 0:2
            problem = struct('A',@(t) A,'b',@(t) b);
            noise = zeros(3,1);
            if port == 1
                problem.noise = @(t) 0.7;
                noise(:) = 0.7;
            elseif port == 2
                noise = [0.1;-0.2;0.3];
                problem.noise = @(t) noise;
            end
            s = A.'*(A*x-b);
            expected = -53*s/norm(s)^p+noise;
            actual = ftcgnn_system(0,x,3,problem,53,p);
            item = struct('trial',trial,'p',p,'port',port, ...
                'formula',norm(actual-expected)/max(1,norm(expected)));
            rows{end+1,1} = item;
        end
    end
end
ftc = struct2table(vertcat(rows{:}));
assert(max(ftc.formula) < 1e-12);
writetable(ftc,fullfile(folder,'ftc.csv'));

problem = struct('G',@(t) eye(2),'h',@(t) zeros(2,1), ...
    'P',@(t) zeros(0,2),'u',@(t) zeros(0,1), ...
    'Q',@(t) zeros(0,2),'v',@(t) zeros(0,1));
zero.gnn = gnn_system(0,zeros(2,1),2,problem,37,1,1e-4);
assert(all(zero.gnn == 0));
problem.G = @(t) zeros(2);
problem.h = @(t) [1;0];
zero.stationary = "";
try
    gnn_system(0,zeros(2,1),2,problem,37,1,1e-4);
catch err
    zero.stationary = string(err.identifier);
end
assert(zero.stationary == "GNN:StationaryResidual");
problem = struct('A',@(t) eye(2),'b',@(t) zeros(2,1),'noise',@(t) [0.2;-0.3]);
zero.ftc = ftcgnn_system(0,zeros(2,1),2,problem,53,1);
assert(isequal(zero.ftc,[0.2;-0.3]));
problem = struct('A',@(t) [1;0],'b',@(t) [0;1]);
zero.inconsistent = ftcgnn_system(0,0,1,problem,53,1);
zero.residual = norm(problem.A(0)*0-problem.b(0));

omega = 1e8;
delta = 1e-4;
mu = [0;5e-13;1e-9;1e-8;2e-8;0.5];
sigma = sqrt(omega^2+mu.^2+delta);
direct = omega+mu-sigma;
stable = (2*omega*mu-delta)./(omega+mu+sigma);
alpha = 1-omega./sigma;
correct = (mu.^2+delta)./(sigma.*(sigma+omega));
pfb = table(mu,direct,stable,alpha,correct);
writetable(pfb,fullfile(folder,'pfb.csv'));

rows = cell(0,1);
weak = cell(0,1);
for ex = 1:2
    if ex == 1
        problem = example1_problem();
    else
        problem = example2_problem();
    end
    low = Inf;
    for t = 0:0.01:10
        clean = current(t,problem);
        g = reference(clean,1e-4);
        [xi,J] = parts(g,clean,1e-4,false,false);
        [eta,L] = parts(g,clean,1e-4,true,true);
        [~,D,V] = svd(J);
        last = D(end,end);
        item = struct('example',ex,'time',t,'sigma',last, ...
            'condition',D(1,1)/last,'direct',norm(xi),'stable',norm(eta), ...
            'gradient',norm(J.'*xi),'gradientstable',norm(L.'*eta), ...
            'difference',norm(J-L,'fro'));
        rows{end+1,1} = item;
        if last < low
            low = last;
            state = g;
            local = clean;
            stamp = t;
            direction = V(:,end);
        end
    end
    [~,J] = parts(state,local,1e-4,false,false);
    [~,L] = parts(state,local,1e-4,true,true);
    for step = [1e-3,1e-4,1e-6,1e-8,1e-10]
        direct = (parts(state+step*direction,local,1e-4,false,false)- ...
            parts(state-step*direction,local,1e-4,false,false))/(2*step);
        stable = (parts(state+step*direction,local,1e-4,true,true)- ...
            parts(state-step*direction,local,1e-4,true,true))/(2*step);
        item = struct('example',ex,'time',stamp,'step',step,'sigma',low, ...
            'predicted',norm(J*direction),'direct',norm(direct), ...
            'stable',norm(stable),'error',norm(direct-J*direction), ...
            'errorstable',norm(stable-L*direction),'difference',norm((J-L)*direction));
        weak{end+1,1} = item;
    end
end
roots = struct2table(vertcat(rows{:}));
weak = struct2table(vertcat(weak{:}));
writetable(roots,fullfile(folder,'roots.csv'));
writetable(weak,fullfile(folder,'weak.csv'));
rows = cell(0,1);
for ex = 1:2
    if ex == 1
        problem = example1_problem();
    else
        problem = example2_problem();
    end
    for t = [0,0.17,1,2.3,5,10]
        clean = current(t,problem);
        g = reference(clean,1e-4)+randn(7,1)/10;
        [~,J] = parts(g,clean,1e-4,false,false);
        [~,L] = parts(g,clean,1e-4,true,true);
        D = zeros(7);
        E = zeros(7);
        for j = 1:7
            step = zeros(7,1);
            step(j) = 1e-6;
            D(:,j) = (parts(g+step,clean,1e-4,false,false)- ...
                parts(g-step,clean,1e-4,false,false))/(2e-6);
            E(:,j) = (parts(g+step,clean,1e-4,true,true)- ...
                parts(g-step,clean,1e-4,true,true))/(2e-6);
        end
        item = struct('example',ex,'time',t, ...
            'direct',norm(J-D,'fro')/norm(J,'fro'), ...
            'stable',norm(L-E,'fro')/norm(L,'fro'));
        rows{end+1,1} = item;
    end
end
moving = struct2table(vertcat(rows{:}));
writetable(moving,fullfile(folder,'moving.csv'));
assert(max(moving.stable) < 1e-6);
rows = cell(0,1);
problem = current(0,example2_problem());
base = reference(problem,1e-4);
for amplitude = [1,10,100,1e4,1e8]
    g = base;
    g(4:7) = amplitude;
    [~,J] = parts(g,problem,1e-4,false,false);
    [~,L] = parts(g,problem,1e-4,true,true);
    item = struct('amplitude',amplitude,'direct',min(svd(J)), ...
        'stable',min(svd(L)));
    rows{end+1,1} = item;
end
limits = struct2table(vertcat(rows{:}));
writetable(limits,fullfile(folder,'global.csv'));
save(fullfile(folder,'audit.mat'),'gnn','ftc','zero','pfb','roots','weak','moving','limits');
fprintf('GNN cases=%d, Jacobian relative error=%g, formula=%g, symmetry=%g\n', ...
    height(gnn),max(gnn.jacobian),max(gnn.formula),max(gnn.symmetry));
fprintf('FTCGNN cases=%d, formula=%g\n',height(ftc),max(ftc.formula));
disp(zero);
disp(pfb);
disp(weak);
disp(moving);
disp(limits);
for ex = 1:2
    row = roots(roots.example == ex,:);
    fprintf('Example %d: min sigma=%g, max condition=%g, max direct residual=%g, stable=%g\n', ...
        ex,min(row.sigma),max(row.condition),max(row.direct),max(row.stable));
end
end

function [xi,J] = parts(g,problem,delta,stable,accurate)
% 按定义计算独立公式。
G = problem.G(0);
S = (G+G.')/2;
h = problem.h(0);
P = problem.P(0);
u = problem.u(0);
Q = problem.Q(0);
v = problem.v(0);
n = size(G,1);
m = size(P,1);
w = size(Q,1);
x = g(1:n);
mu = g(n+m+1:end);
omega = v-Q*x;
if isscalar(delta)
    delta = delta*ones(w,1);
end
sigma = sqrt(omega.^2+mu.^2+delta);
phi = omega+mu-sigma;
if stable
    phi = (2*omega.*mu-delta)./(omega+mu+sigma);
end
xi = [S*x+h+P.'*g(n+1:n+m)+Q.'*mu;P*x-u;phi];
alpha = 1-omega./sigma;
beta = 1-mu./sigma;
if accurate
    take = omega > 0;
    alpha(take) = (mu(take).^2+delta(take))./(sigma(take).*(sigma(take)+omega(take)));
    take = mu > 0;
    beta(take) = (omega(take).^2+delta(take))./(sigma(take).*(sigma(take)+mu(take)));
end
J = [S,P.',Q.';P,zeros(m,m),zeros(m,w); ...
    -diag(alpha)*Q,zeros(w,m),diag(beta)];
end

function clean = current(t,problem)
% 固定题目在当前时刻的数值。
G = problem.G(t);
h = problem.h(t);
P = problem.P(t);
u = problem.u(t);
Q = problem.Q(t);
v = problem.v(t);
clean = struct('G',@(s) G,'h',@(s) h,'P',@(s) P, ...
    'u',@(s) u,'Q',@(s) Q,'v',@(s) v);
end

function g = reference(problem,delta)
% 用等式的一维方向求障碍驻点。
G = problem.G(0);
h = problem.h(0);
P = problem.P(0);
u = problem.u(0);
Q = problem.Q(0);
v = problem.v(0);
base = P.'*u/(P*P.');
dir = [-P(2);P(1)]/norm(P);
a = Q*dir;
b = v-Q*base;
lo = max(b(a < 0)./a(a < 0));
hi = min(b(a > 0)./a(a > 0));
for trial = 1:100
    mid = lo+(hi-lo)/2;
    if mid == lo || mid == hi
        break
    end
    x = base+dir*mid;
    slack = v-Q*x;
    slope = dir.'*(G*x+h)+delta/2*sum(a./slack);
    if slope > 0
        hi = mid;
    else
        lo = mid;
    end
end
x = base+dir*(lo+(hi-lo)/2);
mu = delta./(2*(v-Q*x));
lambda = -P*(G*x+h+Q.'*mu)/(P*P.');
g = [x;lambda;mu];
end
