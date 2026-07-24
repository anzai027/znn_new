function Y_dot = znn_ode( ...
    t, Y, l, problem, ...
    r1, r2, lambda1, lambda2, a, p, q, delta)
%ZNN_ODE  将my_system.m的两个状态合并后交给ODE求解器。

if numel(Y) ~= 2*l
    error('VFCR:WrongStateSize', ...
        'Y must contain [g; z] and therefore have length %d.', 2*l);
end

g = Y(1:l);
z = Y(l+1:2*l);

[g_dot, z_dot] = my_system( ...
    t, g, z, problem, ...
    r1, r2, lambda1, lambda2, a, p, q, delta);

Y_dot = [g_dot; z_dot];
end
