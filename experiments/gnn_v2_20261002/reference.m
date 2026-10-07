function ref = reference(t,problem,delta)
% 独立计算二维等式题目的真解和扰动解。
t = t(:); count = numel(t);
w = size(problem.Q(t(1)),1); l = 3+w;
ref.t = t; ref.x = zeros(count,2); ref.xd = zeros(count,2);
ref.g = zeros(count,l); ref.f = zeros(count,1); ref.fd = zeros(count,1);
ref.check = zeros(count,3);
for k = 1:count
    G = problem.G(t(k)); G = (G+G.')/2;
    h = problem.h(t(k)); P = problem.P(t(k)); u = problem.u(t(k));
    Q = problem.Q(t(k)); v = problem.v(t(k));
    base = P.'*(u/(P*P.'));
    dir = [-P(2);P(1)]/norm(P);
    alpha = -(dir.'*(G*base+h))/(dir.'*G*dir);
    if w == 0
        x = base+dir*alpha; xd = x; mu2 = zeros(0,1);
    else
        a = Q*dir; b = v-Q*base; lo = -Inf; hi = Inf;
        for j = 1:w
            if a(j) > 0, hi = min(hi,b(j)/a(j));
            elseif a(j) < 0, lo = max(lo,b(j)/a(j));
            else, assert(b(j) >= 0);
            end
        end
        assert(isfinite(lo) && isfinite(hi) && lo < hi);
        x = base+dir*min(max(alpha,lo),hi);
        for j = 1:80
            mid = lo+(hi-lo)/2;
            if mid == lo || mid == hi, break, end
            xd = base+dir*mid; slack = v-Q*xd;
            assert(all(slack > 0));
            slope = dir.'*(G*xd+h)+(delta/2)*sum(a./slack);
            if slope > 0, hi = mid; else, lo = mid; end
        end
        xd = base+dir*(lo+(hi-lo)/2);
        mu2 = delta./(2*(v-Q*xd));
    end
    mu1 = -(P*(G*xd+h+Q.'*mu2))/(P*P.');
    ref.x(k,:) = x.'; ref.xd(k,:) = xd.'; ref.g(k,:) = [xd;mu1;mu2].';
    ref.f(k) = 0.5*x.'*G*x+h.'*x; ref.fd(k) = 0.5*xd.'*G*xd+h.'*xd;
    ref.check(k,:) = [abs(P*xd-u),norm(G*xd+h+P.'*mu1+Q.'*mu2), ...
        max([0;Q*x-v;Q*xd-v])];
end
end
