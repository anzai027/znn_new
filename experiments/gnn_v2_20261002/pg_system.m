function [dg,d] = pg_system(t,g,l,problem,gamma,lambda,eps,delta,block)
% 使用阻尼预条件方向和可选分块速度。
assert(isscalar(lambda) && isfinite(lambda) && lambda > 0,'lambda 必须为正数。');
assert(isscalar(eps) && isfinite(eps) && eps >= 0,'eps 必须非负。');
assert(all(isfinite(gamma)) && all(gamma > 0),'gamma 必须为正数。');
assert(isscalar(delta) && delta > 0,'delta 必须为正数。');
g = g(:);
n = size(problem.G(t),1); m = size(problem.P(t),1); w = size(problem.Q(t),1);
assert(l == n+m+w && numel(g) == l,'状态维数不一致。');
[xi,J] = parts(t,g,problem,delta);
d = [J;sqrt(lambda)*eye(l)]\[xi;zeros(l,1)];
if ~block
    assert(isscalar(gamma),'统一速度需要标量 gamma。');
    den = hypot(norm(d),eps);
    dg = zeros(l,1);
    if den > 0, dg = -gamma*d/den; end
else
    assert(numel(gamma) == 3,'分块速度需要三个 gamma。');
    dg = zeros(l,1);
    ids = {1:n,n+1:n+m,n+m+1:l};
    for k = 1:3
        take = ids{k};
        den = hypot(norm(d(take)),eps);
        if den > 0, dg(take) = -gamma(k)*d(take)/den; end
    end
end
end
