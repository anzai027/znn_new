function [m,row] = pp_metrics(r,ref,p)
% 对不规则输出按时间积分计算误差。
c = r.c;
t = r.t;
count = numel(t);
m.residual = zeros(count,1);
m.xd = vecnorm(r.g(:,1:2)-ref.xd,2,2);
m.x = vecnorm(r.g(:,1:2)-ref.x,2,2);
m.g = vecnorm(r.g-ref.g,2,2);
m.predictor = zeros(count,1);
m.oracle = zeros(count,1);
m.sigma = zeros(count,1);
m.lambda = zeros(count,1);
m.alpha = zeros(count,1);
m.speed = zeros(count,1);
m.dp = zeros(count,7);
m.dg = zeros(count,7);
m.equality = zeros(count,1);
m.inequality = zeros(count,1);
names = {'G','h','P','u','Q','v'};
q = struct();
for k = 1:numel(names), q.(names{k}) = p.(names{k}); end
for k = 1:count
    g = r.g(k,:).';
    [dg,a] = pp_system(t(k),g,q,c,r.model);
    [~,b] = pp_system(t(k),ref.g(k,:).',q,c,r.model);
    m.residual(k) = norm(a.xi);
    m.predictor(k) = norm(a.dp-ref.velocity(k,:).');
    m.oracle(k) = norm(b.dp-ref.velocity(k,:).');
    m.sigma(k) = a.sigma;
    m.lambda(k) = a.lambda;
    m.alpha(k) = a.alpha;
    m.dp(k,:) = a.dp.';
    m.dg(k,:) = dg.';
    m.speed(k) = norm(dg);
    m.equality(k) = norm(p.P(t(k))*g(1:2)-p.u(t(k)));
    m.inequality(k) = max([0;p.Q(t(k))*g(1:2)-p.v(t(k))]);
end
take = t >= c.tail;
at = t(take);
weighted = @(v) sqrt(trapz(at,v(take).^2)/(at(end)-at(1)));
settle = @(v,b) crossing(t,v,b);
row = table(r.model,r.complete,r.calls,r.parts,r.seconds,weighted(m.xd), ...
    weighted(m.x),weighted(m.g),max(m.residual(take)),max(m.xd(take)), ...
    max(m.x(take)),max(m.g(take)),weighted(m.predictor),max(m.predictor(take)), ...
    max(m.oracle(take)),max(m.speed),min(m.sigma),min(m.lambda),max(m.lambda), ...
    min(m.alpha),sum(m.alpha < 1-1e-12),max(m.equality(take)),max(m.inequality(take)), ...
    settle(m.residual,1e-4),settle(m.residual,1e-6), ...
    'VariableNames',{'model','complete','calls','parts','seconds','xd','x','g', ...
    'residual','xdpeak','xpeak','gpeak','prediction','predictionpeak','oraclepeak', ...
    'speed','sigma','lambdamin','lambdamax','alpha','limited','equality','inequality','tc4','tc6'});
end

function time = crossing(t,v,b)
% 记录最后一次超阈值之后的采样时刻。
id = find(v > b,1,'last');
if isempty(id), time = t(1);
elseif id == numel(t), time = Inf;
else, time = t(id+1);
end
end
