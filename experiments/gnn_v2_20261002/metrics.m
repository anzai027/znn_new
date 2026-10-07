function m = metrics(r,ref,problem,c)
% 按旧实验口径记录误差和计算成本。
count = numel(r.t); l = r.job.l;
names = {'e','ex','ed','eg','grad','sigma','eq','bound','speed','direction','condition','angle'};
for k = 1:numel(names), m.(names{k}) = nan(count,1); end
for k = 1:count
    g = r.y(k,1:l).';
    if any(~isfinite(g)), continue, end
    [xi,J] = parts(r.t(k),g,problem,c.delta);
    s = J.'*xi; sv = svd(J); error = g-ref.g(k,:).';
    m.e(k) = norm(xi); m.ex(k) = norm(g(1:2)-ref.x(k,:).');
    m.ed(k) = norm(g(1:2)-ref.xd(k,:).'); m.eg(k) = norm(error);
    m.grad(k) = norm(s); m.sigma(k) = min(sv); m.condition(k) = max(sv)/min(sv);
    m.eq(k) = norm(problem.P(r.t(k))*g(1:2)-problem.u(r.t(k)));
    m.bound(k) = max([0;problem.Q(r.t(k))*g(1:2)-problem.v(r.t(k))]);
    if any(r.job.model == ["M3","M4"])
        p = rmfield(problem,{'dG','dh','dP','du','dQ','dv'});
        [dy,d] = pg_system(r.t(k),g,l,p,r.job.rates,r.job.lambda, ...
            r.job.eps,c.delta,r.job.model == "M4");
        m.direction(k) = norm(d);
    elseif r.job.model == "M2"
        dy = -r.job.gamma*s/hypot(norm(s),r.job.eps);
        d = s; m.direction(k) = norm(d);
    elseif r.job.model == "M1"
        dy = zeros(l,1);
        if norm(s) > 0, dy = -r.job.gamma*s/norm(s); end
        d = s; m.direction(k) = norm(d);
    else
        dy = zeros(l,1); d = zeros(l,1);
        v = c.vfcr;
        p = problem;
        if r.job.model == "VFCR-S", p.phi = @(x) fsmooth(x,v,r.job.eps); end
        out = my_system(r.t(k),r.y(k,:).',l,p,v.r1,v.r2,v.lambda1,v.lambda2,v.a,v.p,v.q,c.delta);
        dy = out(1:l);
    end
    m.speed(k) = norm(dy);
    if norm(error)*norm(d) > 0
        m.angle(k) = dot(error,d)/(norm(error)*norm(d));
    end
end
m.tc = [NaN,NaN];
fields = {'rmse','rmsed','peak','peakd','residual','equality','inequality','full'};
for k = 1:numel(fields), m.(fields{k}) = NaN; end
if r.complete
    tail = r.t >= c.tail;
    m.rmse = sqrt(mean(m.ex(tail).^2)); m.rmsed = sqrt(mean(m.ed(tail).^2));
    m.full = sqrt(mean(m.eg(tail).^2)); m.peak = max(m.ex(tail)); m.peakd = max(m.ed(tail));
    m.residual = max(m.e(tail)); m.equality = max(m.eq(tail)); m.inequality = max(m.bound(tail));
    for j = 1:2
        bad = find(m.e > c.threshold(j) | ~isfinite(m.e),1,'last');
        k = 1; if ~isempty(bad), k = bad+1; end
        m.tc(j) = Inf;
        if k <= count && r.t(end)-r.t(k) >= c.dwell-1e-12, m.tc(j) = r.t(k); end
    end
end
m.last = norm(parts(r.last,r.raw.y(end,1:l).',problem,c.delta));
m.motion = motion(r,c);
end
