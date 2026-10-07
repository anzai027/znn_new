function checks = check(c,refs)
% 验证公式实现与独立参考解。
problem = example2_problem(); g = ones(7,1)/2;
[xi,J] = parts(0.3,g,problem,c.delta);
h = 1e-6; D = zeros(7);
for k = 1:7
    e = zeros(7,1); e(k) = h;
    D(:,k) = (parts(0.3,g+e,problem,c.delta)-parts(0.3,g-e,problem,c.delta))/(2*h);
end
checks.jacobian = norm(D-J,'fro')/norm(J,'fro'); assert(checks.jacobian < 1e-7);
problem = rmfield(problem,{'dG','dh','dP','du','dQ','dv'});
[dy,d] = pg_system(0.3,g,7,problem,200,1e-6,1e-2,c.delta,false);
direct = (J.'*J+1e-6*eye(7))\(J.'*xi);
checks.direction = norm(d-direct)/max(1,norm(direct)); assert(checks.direction < 1e-8);
assert(norm(dy) <= 200*(1+1e-12));
dy = pg_system(0.3,g,7,problem,[100,200,200],1e-6,1e-2,c.delta,true);
checks.block = [norm(dy(1:2)),norm(dy(3)),norm(dy(4:7))];
assert(all(checks.block <= [100,200,200]*(1+1e-12)));
checks.derivatives = true;
dy = pg_system(0.3,g,7,problem,200,1e-6,0,c.delta,false);
assert(abs(norm(dy)-200) < 1e-9);
dy = pg_system(0.3,g,7,problem,[100,200,200],1e-6,0,c.delta,true);
assert(all([norm(dy(1:2)),norm(dy(3)),norm(dy(4:7))] <= [100,200,200]*(1+1e-12)));
checks.plain = true;
for k = 1:numel(refs)
    assert(max(refs{k}.check,[],'all') < 1e-10,'参考解检查未通过。');
end
checks.reference = true;
checks.claims = struct();
for k = 1:numel(c.t)
    % 只复查聊天指出的时间点。
    if abs(c.t(k)-6.31) < 1e-12
        old = load(fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
            'gnn_vfcr_smooth_20261002','data','run_010.mat'),'r');
        checks.claims.residual = old.r.metrics.e(k);
        checks.claims.gradient = old.r.metrics.grad(k);
        checks.claims.sigma = min(old.r.metrics.sigma);
    end
end
end
