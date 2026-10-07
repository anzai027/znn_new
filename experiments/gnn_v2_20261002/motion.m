function m = motion(r,c)
% 用统一或分块速度上限排查明显失真。
m = struct('tested',false,'passed',NaN,'rate',NaN,'margin',NaN);
if startsWith(r.job.model,"VFCR"), return, end
m.tested = true;
rate = r.job.gamma;
if r.job.model == "M4", rate = norm(r.job.rates); end
noise = 0;
if r.job.port == "state"
    switch r.bank.kind
        case "constant", noise = 2*sqrt(r.job.l);
        case "linear", noise = 0.1*c.t(end)*sqrt(r.job.l);
        case {"cosine","harmonic"}, noise = sqrt(r.job.l);
        case {"gamma","gaussian"}, noise = max(vecnorm(r.bank.w,2,2));
    end
end
m.rate = rate+noise;
m.margin = -Inf;
for source = 1:2
    if source == 1, t = r.raw.t; g = r.raw.y(:,1:r.job.l);
    else, t = r.t; g = r.y(:,1:r.job.l); end
    take = all(isfinite(g),2); t = t(take); g = g(take,:);
    scale = max([ones(numel(t),1),vecnorm(g,2,2),norm(r.job.g0)*ones(numel(t),1)],[],2);
    tol = 100*(r.job.atol+r.job.rtol*scale);
    margin = max(vecnorm(g-r.job.g0.',2,2)-m.rate*t-tol);
    if numel(t) > 1
        margin = max(margin,max(vecnorm(diff(g),2,2)-m.rate*diff(t)-max(tol(1:end-1),tol(2:end))));
    end
    m.margin = max(m.margin,margin);
end
m.passed = double(m.margin <= 0);
end
