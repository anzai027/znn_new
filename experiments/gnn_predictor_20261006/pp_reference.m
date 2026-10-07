function ref = pp_reference(t,p,delta)
% 独立计算二维约束题目的解与参考速度。
t = t(:);
assert(isreal(t) && all(isfinite(t)),'参考时间必须为有限实数。');
assert(isreal(delta) && all(isfinite(delta(:))) && all(delta(:) > 0), ...
    'delta 必须为正有限实数。');
if isscalar(delta), delta = delta*ones(4,1); else, delta = delta(:); end
assert(numel(delta) == 4,'delta 必须为标量或四维向量。');
count = numel(t);
ref.t = t;
ref.delta = delta;
ref.x = zeros(count,2);
ref.xd = zeros(count,2);
ref.g = zeros(count,7);
ref.velocity = zeros(count,7);
ref.speed = zeros(count,1);
ref.sigma = zeros(count,1);
ref.bias = zeros(count,1);
ref.error = zeros(count,1);
ref.f = zeros(count,1);
ref.fd = zeros(count,1);
ref.slack = zeros(count,4);
ref.check = zeros(count,5);
ref.derivative = zeros(count,1);
ref.identity = zeros(count,1);
ref.columns = ["true equality","delta equality","stationarity", ...
    "inequality","PFB"];
if count == 0, return, end
assert(isequal(size(p.G(t(1))),[2,2]) && isequal(size(p.P(t(1))),[1,2]) ...
    && isequal(size(p.Q(t(1))),[4,2]),'参考题目必须为二维、一条等式和四条不等式。');
names = {'dG','dh','dP','du','dQ','dv'};
assert(all(isfield(p,names)),'解析导数只用于参考轨迹核验。');
for k = 1:count
    time = t(k);
    G = p.G(time); G = (G+G.')/2;
    h = p.h(time); h = h(:);
    P = p.P(time); u = p.u(time);
    Q = p.Q(time); v = p.v(time); v = v(:);
    normp = norm(P);
    assert(normp > 0,'等式矩阵必须满行秩。');
    base = P.'*(u/(normp^2));
    dir = [-P(2);P(1)]/normp;
    a = Q*dir; b = v-Q*base;
    lo = -Inf; hi = Inf;
    for j = 1:4
        if a(j) > 0, hi = min(hi,b(j)/a(j));
        elseif a(j) < 0, lo = max(lo,b(j)/a(j));
        else, assert(b(j) > 0,'题目必须存在严格可行点。');
        end
    end
    assert(isfinite(lo) && isfinite(hi) && lo < hi,'可行区间必须有内部。');
    curvature = dir.'*G*dir;
    slope = dir.'*(G*base+h);
    assert(curvature > 0,'目标沿可行方向必须严格凸。');
    alpha = -slope/curvature;
    x = base+dir*min(max(alpha,lo),hi);
    for j = 1:80
        mid = lo+(hi-lo)/2;
        if mid == lo || mid == hi, break, end
        slack = b-a*mid;
        assert(all(slack > 0),'障碍解的松弛必须为正。');
        value = curvature*mid+slope+sum((delta/2).*a./slack);
        if value > 0, hi = mid; else, lo = mid; end
    end
    xd = base+dir*(lo+(hi-lo)/2);
    slack = v-Q*xd;
    assert(all(slack > 0),'扰动解必须严格可行。');
    mu = delta./(2*slack);
    force = G*xd+h+Q.'*mu;
    mu1 = -(P*force)/(normp^2);
    g = [xd;mu1;mu];
    sigma = hypot(hypot(slack,mu),sqrt(delta));
    a = (mu.^2+delta)./(sigma.*(sigma+slack));
    b = (slack.^2+delta)./(sigma.*(sigma+mu));
    phi = (2*slack.*mu-delta)./(slack+mu+sigma);
    xi = [force+P.'*mu1;P*xd-u;phi];
    J = [G,P.',Q.';P,0,zeros(1,4);-a.*Q,zeros(4,1),diag(b)];
    dG = p.dG(time); dG = (dG+dG.')/2;
    omega = p.dv(time)-p.dQ(time)*xd;
    first = dG*xd+p.dh(time)+p.dP(time).'*mu1+p.dQ(time).'*mu;
    second = p.dP(time)*xd-p.du(time);
    ft = [first;second;a.*omega];
    D = delta./(2*slack.^2);
    H = G+Q.'*(D.*Q);
    part = -[H,P.';P,0]\[first-Q.'*(D.*omega);second];
    velocity = [part;-D.*(omega-Q*part(1:2))];
    full = -J\ft;
    ref.x(k,:) = x.';
    ref.xd(k,:) = xd.';
    ref.g(k,:) = g.';
    ref.velocity(k,:) = velocity.';
    ref.speed(k) = norm(velocity);
    ref.sigma(k) = min(svd(J));
    ref.bias(k) = norm(xd-x);
    ref.error(k) = norm(xi);
    ref.f(k) = 0.5*x.'*G*x+h.'*x;
    ref.fd(k) = 0.5*xd.'*G*xd+h.'*xd;
    ref.slack(k,:) = slack.';
    ref.check(k,:) = [abs(P*x-u),abs(P*xd-u),norm(xi(1:2)), ...
        max([0;Q*x-v;Q*xd-v]),norm(phi)];
    ref.derivative(k) = norm(full-velocity)/max(1,norm(velocity));
    ref.identity(k) = norm(J*velocity+ft)/max(1,norm(ft));
end
end
