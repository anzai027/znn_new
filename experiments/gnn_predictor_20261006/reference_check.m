function checks = reference_check()
% 核验独立参考解及乘子速度尖峰。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(folder,'-begin'); addpath(root,'-begin');
p = example2_problem(); delta = 1e-4;
peak = 3.4611091472136;
t = [0;0.3;1;3.46;3.461;peak;3.46111;5;6.31;peak+2*pi;10];
ref = pp_reference(t,p,delta);
checks.feasible = max(ref.check,[],'all');
checks.residual = max(ref.error);
checks.derivative = max(ref.derivative);
checks.identity = max(ref.identity);
assert(checks.feasible < 1e-9 && checks.residual < 1e-9);
assert(checks.derivative < 1e-9 && checks.identity < 1e-9);
opt = optimset('TolX',1e-13,'Display','off');
time = fminbnd(@rate,peak-2e-7,peak+2e-7,opt);
point = pp_reference(time,p,delta);
checks.time = time;
checks.speed = point.speed;
assert(abs(checks.speed-4735.7083211) < 0.01,'参考速度峰与独立复核不一致。');
h = [1e-6;1e-7;1e-8;1e-9];
velocity = zeros(numel(h),7); relative = zeros(numel(h),1);
for k = 1:numel(h)
    pair = pp_reference([time-h(k);time+h(k)],p,delta);
    velocity(k,:) = (pair.g(2,:)-pair.g(1,:))/(2*h(k));
    relative(k) = norm(velocity(k,:)-point.velocity)/point.speed;
end
checks.difference = relative;
assert(relative(3) < 1e-6,'参考速度的中心差分核验未通过。');
checks.periodic = norm(ref.g(6,:)-ref.g(10,:));
assert(checks.periodic < 1e-9,'周期相同的参考根不一致。');
checks.vector = pp_reference(0.3,p,delta*ones(4,1)).error;
checks.empty = size(pp_reference([],p,delta).g);
assert(isequal(checks.empty,[0,7]));
check = table(h,vecnorm(velocity,2,2),relative, ...
    'VariableNames',{'h','speed','relative'});
writetable(check,fullfile(folder,'reference_difference.csv'));
save(fullfile(folder,'reference_check.mat'),'checks','ref','point','check');
disp(checks); disp(check);

    function value = rate(t)
        r = pp_reference(t,p,delta); value = -r.speed;
    end
end
