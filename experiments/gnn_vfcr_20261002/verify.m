function v = verify(c,problems,refs)
% 检查参考解、题目符号和无导数调用。
v = struct();
for j = 1:2
    r = refs{j};
    assert(max(r.eq) < 1e-10,'参考解违反等式。');
    assert(max(r.bound) < 1e-10,'参考解违反边界。');
    assert(max(r.station) < 1e-9,'障碍参考解不满足驻点条件。');
    assert(max(r.pfb) < 1e-12,'障碍参考解不满足 PFB 条件。');
    assert(max(r.complement) < 1e-12,'障碍参考解不满足扰动互补条件。');
    assert(max(r.raw) < 1e-7,'参考解的原始残差过大。');
end
v.equality = max(cellfun(@(r) max(r.eq),refs));
v.stationarity = max(cellfun(@(r) max(r.station),refs));
v.raw = max(cellfun(@(r) max(r.raw),refs));
v.jacobian = 0;
problem = problems{2};
clean = rmfield(problem,{'dG','dh','dP','du','dQ','dv'});
g = (1:c.l).'/10;
for t = [0,0.17,1,2.3,5,10]
    x = reference(t,problem,c.delta);
    assert(abs(sin(4*t)*x.x(1)-cos(4*t)*x.x(2)-cos(2*t)) < 1e-10, ...
        '等式约束的负号有误。');
    [xi,J] = parts(t,g,problem,c.delta);
    D = zeros(c.l);
    for k = 1:c.l
        h = zeros(c.l,1);
        h(k) = 1e-6;
        D(:,k) = (parts(t,g+h,problem,c.delta)-parts(t,g-h,problem,c.delta))/(2e-6);
    end
    v.jacobian = max(v.jacobian,norm(J-D,'fro')/max(1,norm(J,'fro')));
    dg = gnn_system(t,g,c.l,clean,100,1,c.delta);
    assert(all(isfinite(dg)) && (J.'*xi).'*dg < 0,'梯度方向不正确。');
    assert(abs(norm(dg)-100) < 1e-9,'归一化速度不正确。');
end
assert(v.jacobian < 1e-5,'雅可比与独立差分不一致。');
v.derivatives = "新模型已在删除全部六个导数字段后成功调用。";
v.passed = true;
end
