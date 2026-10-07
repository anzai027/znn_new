function m = metrics(r,ref,problem,c)
% 使用同一口径评价完整和中断的运行。
m = curves(r.t,r.y(:,1:c.l),ref,problem,c.delta);
m.tc = [NaN,NaN];
m.rmse = NaN;
m.rmsed = NaN;
m.peak = NaN;
m.peakd = NaN;
m.residual = NaN;
m.equality = NaN;
m.inequality = NaN;
m.bias = max(vecnorm(ref.x-ref.xd,2,2));
if r.complete
    for j = 1:2
        bad = find(m.e > c.threshold(j) | ~isfinite(m.e),1,'last');
        if isempty(bad)
            k = 1;
        else
            k = bad+1;
        end
        m.tc(j) = Inf;
        if k <= numel(r.t) && r.t(end)-r.t(k) >= c.dwell-1e-12
            m.tc(j) = r.t(k);
        end
    end
    tail = r.t >= c.tail;
    m.rmse = sqrt(mean(m.ex(tail).^2));
    m.rmsed = sqrt(mean(m.ed(tail).^2));
    m.peak = max(m.ex(tail));
    m.peakd = max(m.ed(tail));
    m.residual = max(m.e(tail));
    m.equality = max(m.eq(tail));
    m.inequality = max(m.bound(tail));
end
take = unique(round(linspace(1,numel(r.raw.t),min(512,numel(r.raw.t)))));
raw = reference(r.raw.t(take),problem,c.delta);
m.raw = curves(raw.t,r.raw.y(take,1:c.l),raw,problem,c.delta);
m.last = m.raw.e(end);
m.lastx = m.raw.ex(end);
m.lastd = m.raw.ed(end);
m.lastgrad = m.raw.grad(end);
m.min = min(m.raw.sigma);
step = diff(r.raw.t);
m.steps = numel(step);
m.smallest = NaN;
if ~isempty(step)
    m.smallest = min(step);
end
m.motion = motion(r,c);
end

function m = curves(t,g,ref,problem,delta)
m.t = t;
m.x = g(:,1:2);
names = {'e','ex','ed','eg','grad','sigma','eq','bound','dual','comp','f','gap'};
for k = 1:numel(names)
    m.(names{k}) = nan(numel(t),1);
end
for k = 1:numel(t)
    if any(~isfinite(g(k,:)))
        continue
    end
    v = g(k,:).';
    x = v(1:2);
    [xi,J] = parts(t(k),v,problem,delta);
    m.e(k) = norm(xi);
    m.ex(k) = norm(x-ref.x(k,:).');
    m.ed(k) = norm(x-ref.xd(k,:).');
    m.eg(k) = norm(v-ref.g(k,:).');
    m.grad(k) = norm(J.'*xi);
    m.sigma(k) = min(svd(J));
    m.eq(k) = norm(problem.P(t(k))*x-problem.u(t(k)));
    slack = problem.v(t(k))-problem.Q(t(k))*x;
    m.bound(k) = max([0;-slack]);
    m.dual(k) = max([0;-v(4:7)]);
    m.comp(k) = norm(slack.*v(4:7));
    G = problem.G(t(k));
    h = problem.h(t(k));
    m.f(k) = 0.5*x.'*G*x+h.'*x;
    m.gap(k) = m.f(k)-ref.f(k);
end
end
