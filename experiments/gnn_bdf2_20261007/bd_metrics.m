function [m,row] = bd_metrics(r,ref,p)
% 分开统计扰动根跟踪误差和原题误差。
c = r.c; t = c.t; vals = nan(1,11); m = struct();
if r.complete
    m.xd = vecnorm(r.g(:,1:2)-ref.xd,2,2);
    m.x = vecnorm(r.g(:,1:2)-ref.x,2,2);
    m.g = vecnorm(r.g-ref.g,2,2);
    m.residual = zeros(numel(t),1); m.obs = m.residual;
    m.speed = m.residual; m.prediction = m.residual;
    for j = 1:numel(t)
        [dy,a] = bd_system(t(j),r.y(j,:).',p,c,r.model);
        m.residual(j) = norm(a.xi); m.obs(j) = norm(a.obs);
        m.speed(j) = norm(dy(1:7));
        m.prediction(j) = norm(a.dp-ref.velocity(j,:).');
    end
    take = t >= c.tail; at = t(take);
    rms = @(e) sqrt(trapz(at,e(take).^2)/(at(end)-at(1)));
    vals = [rms(m.xd),rms(m.x),rms(m.g),max(m.g(take)), ...
        rms(m.residual),max(m.residual(take)),rms(m.obs), ...
        rms(m.prediction),max(m.prediction(take)),max(m.speed),rms(ref.bias)];
end
row = cell2table([{r.group,r.label,r.model,c.order,c.delta,c.h,c.kind,c.amp,c.seed, ...
    c.rtol,c.atol,r.complete,r.status,r.last,r.calls,r.parts,r.seconds},num2cell(vals),{r.reason}], ...
    'VariableNames',{'group','label','model','order','delta','h','noise','amp','seed', ...
    'rtol','atol','complete','status','last','calls','parts','seconds', ...
    'xd','x','g','gpeak','residual','respeak','obs','prediction','predpeak','speed','bias','reason'});
end
