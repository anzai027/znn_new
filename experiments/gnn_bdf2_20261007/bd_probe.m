function bd_probe()
% 分开检查差分截断误差和独立样本噪声放大。
folder = bd_setup(); p = example2_problem(); c = bd_config(); rows = table();
for delta = [1e-4,1e-5,1e-6]
    ci = c; ci.delta = delta;
    ref = pp_reference([3.4611091472136;9.74429445439319],p,delta);
    for j = 1:2
        t = ref.t(j); g = ref.g(j,:).';
        for order = [1,2]
            for h = [1e-4,1e-5,1e-6]
                ci.order = order; ci.h = h;
                [~,a] = bd_system(t,g,p,ci,"M8a-BDF2");
                exact = -a.J*ref.velocity(j,:).';
                row = table(delta,t,order,h,norm(a.ft-exact), ...
                    norm(a.dp-ref.velocity(j,:).'),a.sigma,norm(ref.velocity(j,:)), ...
                    'VariableNames',{'delta','t','order','h','timeerror','velocityerror','sigma','speed'});
                rows = [rows;row];
            end
        end
    end
end
writetable(rows,fullfile(folder,'prediction.csv'));
s = RandStream('mt19937ar','Seed',29); count = 100000;
n = randn(s,7,count,3); a = n(:,:,1)-n(:,:,2);
b = (3*n(:,:,1)-4*n(:,:,2)+n(:,:,3))/2;
ref = pp_reference(9.74429445439319,p,c.delta); [~,J] = pp_parts(ref.t,ref.g.',p,c.delta);
rows2 = table();
for h = [1e-4,1e-5,1e-6]
    amp = 1e-5;
    fd = amp*a/h; bd = amp*b/h;
    expected1 = amp*sqrt(2)/h; expected2 = amp*sqrt(6.5)/h;
    v1 = J\fd; v2 = J\bd;
    row = table(h,amp,sqrt(mean(fd(:).^2)),sqrt(mean(bd(:).^2)), ...
        expected1,expected2,sqrt(mean(sum(v1.^2,1))),sqrt(mean(sum(v2.^2,1))), ...
        'VariableNames',{'h','amp','fd','bdf','expectedfd','expectedbdf','velocityfd','velocitybdf'});
    rows2 = [rows2;row];
end
assert(all(abs(rows2.fd./rows2.expectedfd-1)<0.01));
assert(all(abs(rows2.bdf./rows2.expectedbdf-1)<0.01));
writetable(rows2,fullfile(folder,'iid.csv'));
save(fullfile(folder,'probe.mat'),'rows','rows2'); disp(rows); disp(rows2);
end
