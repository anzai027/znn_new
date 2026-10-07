function bd_verify()
% 独立重算误差并检查更密网格是否改变结论。
folder = bd_setup(); p = example2_problem(); c = bd_config();
t = unique([c.t;(0:0.00025:10).']);
for peak = c.peak
    t = unique([t;(peak-0.001:5e-7:peak+0.001).'; ...
        (peak-0.0001:5e-8:peak+0.0001).']);
end
refs = cell(3,1); smooth = cell(3,1); deltas = [1e-4,1e-5,1e-6];
check = table();
cache = fullfile(folder,'reference_fine.mat');
codes = [bd_hash(which('pp_reference')),bd_hash(which('pp_parts')),bd_hash(which('bd_polish'))];
if isfile(cache)
    saved = load(cache); assert(isequal(saved.t,t) && isequal(saved.codes,codes));
    refs = saved.refs; smooth = saved.smooth; check = saved.check;
else
for j = 1:3
    refs{j} = pp_reference(t,p,deltas(j));
    ref = refs{j}; take = t>=5;
    smooth{j} = bd_polish(ref,p,deltas(j)); alt = smooth{j};
    bias = sqrt(trapz(t(take),ref.bias(take).^2)/5);
    row = table(deltas(j),max(ref.error),min(ref.sigma),max(ref.speed), ...
        bias,max(ref.derivative),max(alt.polished),max(alt.shift),max(alt.xshift), ...
        'VariableNames',{'delta','residual','sigma','speed','bias','identity', ...
        'polished','gshift','xshift'});
    check = [check;row];
end
save(cache,'refs','smooth','check','t','codes','-v7.3');
end
writetable(check,fullfile(folder,'reference_check.csv')); rows = table();
for k = 1:37
    file = fullfile(folder,'results',sprintf('run_%02d.mat',k));
    if ~isfile(file), break, end
    data = load(file,'r','row'); r = data.r;
    if ~r.complete, continue, end
    assert(isequal(r.g,r.y(:,1:7)) && isequal(r.t,r.c.t));
    ci = r.c; ref = refs{find(deltas==ci.delta)};
    y = values(r,t); old = values(r,r.t);
    difference = max(abs(old-r.y),[],'all'); assert(difference==0);
    xd = vecnorm(y(:,1:2)-ref.xd,2,2); x = vecnorm(y(:,1:2)-ref.x,2,2);
    g = vecnorm(y(:,1:7)-ref.g,2,2); take = t>=ci.tail;
    rms = @(v) sqrt(trapz(t(take),v(take).^2)/(t(find(take,1,'last'))-ci.tail));
    fine = [rms(xd),rms(x),rms(g),max(g(take))];
    alt = smooth{find(deltas==ci.delta)};
    xdp = vecnorm(y(:,1:2)-alt.xd,2,2); gp = vecnorm(y(:,1:7)-alt.g,2,2);
    source = load(fullfile(folder,'results',sprintf('ref_%g.mat',ci.delta)),'ref');
    xd0 = vecnorm(r.g(:,1:2)-source.ref.xd,2,2); g0 = vecnorm(r.g-source.ref.g,2,2);
    at = r.t >= ci.tail; z = r.t(at);
    xdcheck = sqrt(trapz(z,xd0(at).^2)/(z(end)-z(1)));
    assert(abs(xdcheck-data.row.xd)<=max(1e-20,1e-12*data.row.xd));
    assert(isequal(xd0,r.metrics.xd) && isequal(g0,r.metrics.g));
    row = table(k,difference,fine(1),fine(2),fine(3),fine(4), ...
        abs(fine(1)/data.row.xd-1),abs(fine(4)/data.row.gpeak-1),rms(xdp),max(gp(take)), ...
        'VariableNames',{'run','deval','xd','x','g','gpeak','rmsechange','peakchange','xdp','gpeakp'});
    rows = [rows;row]; fprintf('验证 %d/37\n',k);
end
writetable(rows,fullfile(folder,'verification.csv'));
complete = height(rows)==37;
save(fullfile(folder,'verification.mat'),'rows','check','complete'); disp(check);
fprintf('验证完成 %d/37\n',height(rows));
end
function y = values(r,t)
% 按实际分段积分结果读取轨迹。
ends = [r.c.start,r.c.start+r.c.h,r.c.start+2*r.c.h,r.c.finish];
y = zeros(numel(t),size(r.y,2));
for j = 1:3
    take = t>=ends(j)&t<=ends(j+1);
    if any(take), y(take,:) = deval(r.sols{j},t(take)).'; end
end
end
