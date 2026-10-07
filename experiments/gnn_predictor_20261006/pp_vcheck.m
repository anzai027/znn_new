function rows = pp_vcheck()
% 独立核对 VFCR 与新模型的统一口径。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(root,'-begin'); addpath(folder,'-begin');
out = fullfile(folder,'vfcr');
main = load(fullfile(folder,'data','all.mat'),'results','c','ref');
extra = load(fullfile(folder,'extra','all.mat'),'results');
a = load(fullfile(out,'run_01.mat')); b = load(fullfile(out,'run_02.mat'));
assert(strcmp(a.hash.my_system,digest(fullfile(root,'my_system.m'))));
assert(strcmp(a.hash.fsmooth,digest(fullfile(root,'experiments','gnn_v2_20261002','fsmooth.m'))));
assert(isequal(a.hash,b.hash));
c = main.c; p = example2_problem();
t = unique([c.t;(0:0.00025:10).']);
ref = pp_reference(t,p,c.delta);
assert(max(ref.error) < 1e-9);
list = {a.r,b.r,main.results{6},extra.results{1}};
labels = ["VFCR-S main","VFCR-S strict","M5 strict","M6 kappa=.01 strict"];
rows = table();
for k = 1:4
    r = list{k}; assert(r.complete && r.status == "complete");
    assert(isequal(r.t,c.t) && isequal(r.c.g0,c.g0) && r.c.delta == c.delta);
    assert(r.c.step == c.step && r.c.tail == c.tail);
    if k > 1, assert(r.c.rtol == 1e-9 && r.c.atol == 1e-11); end
    if k <= 2
        assert(size(r.y,2) == 14 && isequal(r.g,r.y(:,1:7)));
        assert(all(r.y(1,8:14) == 0));
    else
        assert(size(r.g,2) == 7);
    end
    dense = deval(r.sol,t).'; dense = dense(:,1:7);
    xd = vecnorm(dense(:,1:2)-ref.xd,2,2);
    x = vecnorm(dense(:,1:2)-ref.x,2,2);
    full = vecnorm(dense-ref.g,2,2);
    old = [rms(r.t,r.metrics.xd,c.tail),rms(r.t,r.metrics.x,c.tail),rms(r.t,r.metrics.g,c.tail)];
    fine = [rms(t,xd,c.tail),rms(t,x,c.tail),rms(t,full,c.tail)];
    take = r.t >= c.tail;
    row = table(labels(k),r.c.rtol,r.c.atol,r.calls,r.seconds,old(1),old(2),old(3), ...
        max(r.metrics.g(take)),max(r.metrics.residual(take)),fine(1),fine(2),fine(3), ...
        max(abs(fine-old)./max(old,realmin)), ...
        'VariableNames',{'model','rtol','atol','calls','seconds','xd','x','g','gpeak', ...
        'residual','refinedxd','refinedx','refinedg','refinement'});
    rows = [rows;row];
end
writetable(rows,fullfile(out,'comparison.csv'));
checks = struct('same',true,'complete',true,'hash',true,'points',numel(t), ...
    'reference',max(ref.error),'bias',rms(t,ref.bias,c.tail));
save(fullfile(out,'verification.mat'),'checks','rows');
disp(checks); disp(rows);
end

function value = rms(t,e,start)
% 对密集和不规则时间网格使用相同权重。
take = t >= start; t = t(take); e = e(take);
value = sqrt(trapz(t,e.^2)/(t(end)-t(1)));
end

function code = digest(path)
% 核对运行时使用的旧模型。
file = fopen(path,'rb'); clean = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8'); md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes); code = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
