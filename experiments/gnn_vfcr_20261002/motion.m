function m = motion(r,c)
% 用 p=1 的速度上限检查轨迹是否明显失真。
m = struct('tested',false,'rate',NaN,'raw',NaN,'grid',NaN, ...
    'jump',NaN,'output',NaN,'excess',NaN,'margin',NaN,'passed',NaN,'factor',100);
if r.job.model ~= "GNN" || r.job.p ~= 1
    return
end
noise = 0;
if r.job.port == "state"
    switch r.bank.kind
        case "constant", noise = 2*sqrt(c.l);
        case "linear", noise = 0.1*c.t(end)*sqrt(c.l);
        case {"cosine","harmonic"}, noise = sqrt(c.l);
        case {"gamma","gaussian"}, noise = max(vecnorm(r.bank.w,2,2));
    end
end
m.tested = true;
m.rate = r.job.gamma+noise;
[m.raw,m.jump,a] = measure(r.raw.t,r.raw.y(:,1:c.l));
[m.grid,m.output,b] = measure(r.t,r.y(:,1:c.l));
m.excess = max([m.raw,m.grid,m.jump,m.output]);
m.margin = max(a,b);
m.passed = double(m.margin <= 0);

    function [excess,jump,margin] = measure(t,g)
        take = all(isfinite(g),2);
        t = t(take);
        g = g(take,:);
        distance = vecnorm(g-r.job.g0.',2,2);
        limit = m.rate*t;
        scale = max([ones(numel(t),1),vecnorm(g,2,2), ...
            norm(r.job.g0)*ones(numel(t),1)],[],2);
        tolerance = m.factor*(r.job.atol+r.job.rtol*scale);
        excess = max(distance-limit);
        margin = max(distance-limit-tolerance);
        jump = 0;
        if numel(t) > 1
            difference = vecnorm(diff(g),2,2)-m.rate*diff(t);
            allowance = max(tolerance(1:end-1),tolerance(2:end));
            jump = max(difference);
            margin = max(margin,max(difference-allowance));
        end
    end
end
