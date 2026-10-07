function [m,row] = pn_metrics(r,ref,p)
% 在同一密集网格上按时间积分统计误差。
c = r.c; t = r.t; count = numel(t);
m.xd = vecnorm(r.g(:,1:2)-ref.xd,2,2);
m.x = vecnorm(r.g(:,1:2)-ref.x,2,2);
m.g = vecnorm(r.g-ref.g,2,2);
fields = {'residual','predictor','oracle','sigma','lc','lp','speed', ...
    'balance','feedback','condition','equality','inequality'};
for k = 1:numel(fields), m.(fields{k}) = zeros(count,1); end
m.dp = zeros(count,7); m.dg = m.dp;
names = {'G','h','P','u','Q','v'}; q = struct();
for k = 1:numel(names), q.(names{k}) = p.(names{k}); end
for k = 1:count
    g = r.g(k,:).';
    [dg,a] = pn_system(t(k),g,q,c,r.model);
    [~,b] = pn_system(t(k),ref.g(k,:).',q,c,r.model);
    m.residual(k) = norm(a.xi);
    m.predictor(k) = norm(a.dp-ref.velocity(k,:).');
    m.oracle(k) = norm(b.dp-ref.velocity(k,:).');
    m.sigma(k) = a.sigma; m.lc(k) = a.lc; m.lp(k) = a.lp;
    m.dp(k,:) = a.dp.'; m.dg(k,:) = dg.'; m.speed(k) = norm(dg);
    m.balance(k) = a.balance; m.feedback(k) = norm(a.feedback);
    m.condition(k) = a.condition;
    m.equality(k) = norm(p.P(t(k))*g(1:2)-p.u(t(k)));
    m.inequality(k) = max([0;p.Q(t(k))*g(1:2)-p.v(t(k))]);
end
take = t >= c.tail; at = t(take);
weight = @(v) sqrt(trapz(at,v(take).^2)/(at(end)-at(1)));
row = table(r.group,r.label,r.model,c.rtol,c.atol,c.h,c.correct, ...
    r.complete,r.calls,r.parts,r.seconds,weight(m.xd),weight(m.x),weight(m.g), ...
    max(m.residual(take)),weight(m.residual),max(m.xd(take)),max(m.g(take)), ...
    weight(m.predictor),max(m.predictor(take)),max(m.oracle(take)), ...
    max(m.speed),max(m.feedback),min(m.sigma),min(m.lc),max(m.lc), ...
    min(m.lp),max(m.lp),max(m.balance(take)),min(m.condition), ...
    max(m.equality(take)),max(m.inequality(take)), ...
    crossing(t,m.residual,1e-4),crossing(t,m.residual,1e-6), ...
    'VariableNames',{'group','label','model','rtol','atol','h','correct', ...
    'complete','calls','parts','seconds','xd','x','g','residual','residualrmse', ...
    'xdpeak','gpeak','prediction','predictionpeak','oraclepeak','speed', ...
    'feedback','sigma','lcmin','lcmax','lpmin','lpmax','balance','condition', ...
    'equality','inequality','tc4','tc6'});
end

function time = crossing(t,v,b)
% 记录最后一次超阈值之后的采样时刻。
id = find(v > b,1,'last');
if isempty(id), time = t(1);
elseif id == numel(t), time = Inf;
else, time = t(id+1);
end
end
