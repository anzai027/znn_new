function r = reference(t,problem,delta)
% 用一维问题独立计算两种参考解。
t = t(:);
n = numel(t);
r.t = t;
r.x = zeros(n,2);
r.xd = zeros(n,2);
r.g = zeros(n,7);
r.f = zeros(n,1);
r.fd = zeros(n,1);
r.eq = zeros(n,1);
r.bound = zeros(n,1);
r.station = zeros(n,1);
r.pfb = zeros(n,1);
r.complement = zeros(n,1);
r.raw = zeros(n,1);
for k = 1:n
    G = problem.G(t(k));
    G = (G+G.')/2;
    h = problem.h(t(k));
    P = problem.P(t(k));
    u = problem.u(t(k));
    Q = problem.Q(t(k));
    v = problem.v(t(k));
    base = P.'*(u/(P*P.'));
    dir = [-P(2);P(1)]/norm(P);
    a = Q*dir;
    b = v-Q*base;
    lo = -Inf;
    hi = Inf;
    for j = 1:numel(a)
        if a(j) > 0
            hi = min(hi,b(j)/a(j));
        elseif a(j) < 0
            lo = max(lo,b(j)/a(j));
        else
            assert(b(j) >= 0,'参考问题不可行。');
        end
    end
    assert(isfinite(lo) && isfinite(hi) && lo < hi,'参考区间无效。');
    alpha = -(dir.'*(G*base+h))/(dir.'*G*dir);
    alpha = min(max(alpha,lo),hi);
    x = base+dir*alpha;
    for j = 1:80
        mid = lo+(hi-lo)/2;
        if mid == lo || mid == hi
            break
        end
        xd = base+dir*mid;
        slack = v-Q*xd;
        assert(all(slack > 0),'障碍参考解不在内部。');
        slope = dir.'*(G*xd+h)+(delta/2)*sum(a./slack);
        if slope > 0
            hi = mid;
        else
            lo = mid;
        end
    end
    xd = base+dir*(lo+(hi-lo)/2);
    slack = v-Q*xd;
    mu2 = delta./(2*slack);
    mu1 = -(P*(G*xd+h+Q.'*mu2))/(P*P.');
    g = [xd;mu1;mu2];
    sigma = sqrt(slack.^2+mu2.^2+delta);
    stable = (2*slack.*mu2-delta)./(slack+mu2+sigma);
    r.x(k,:) = x.';
    r.xd(k,:) = xd.';
    r.g(k,:) = g.';
    r.f(k) = 0.5*x.'*G*x+h.'*x;
    r.fd(k) = 0.5*xd.'*G*xd+h.'*xd;
    r.eq(k) = max(abs([P*x-u;P*xd-u]));
    r.bound(k) = max([0;Q*x-v;Q*xd-v]);
    r.station(k) = norm(G*xd+h+P.'*mu1+Q.'*mu2);
    r.pfb(k) = norm(stable);
    r.complement(k) = max(abs(slack.*mu2-delta/2));
    r.raw(k) = norm(parts(t(k),g,problem,delta));
end
end
