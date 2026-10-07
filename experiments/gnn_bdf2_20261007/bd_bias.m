function bd_bias(name,tol,step)
% 用自适应积分复核小 delta 的原题误差。
folder = bd_setup(); p = example2_problem(); rows = table();
if nargin < 1, name = 'bias_adaptive'; end
if nargin < 2, tol = 1e-6; end
if nargin < 3, step = 0.02; end
deltas = [1e-4,1e-5,1e-6]; ids = [3,8,11];
for k = 1:3
    delta = deltas(k); data = load(fullfile(folder,'results',sprintf('run_%02d.mat',ids(k))),'r');
    r = data.r;
    [bias,error1] = quadgk(@value,5,10,'Waypoints',5:step:10, ...
        'RelTol',tol,'AbsTol',1e-20,'MaxIntervalCount',10000);
    [original,error2] = quadgk(@actual,5,10,'Waypoints',5:step:10, ...
        'RelTol',tol,'AbsTol',1e-20,'MaxIntervalCount',10000);
    row = table(delta,sqrt(bias/5),sqrt(original/5),error1,error2, ...
        'VariableNames',{'delta','bias','x','biasbound','xbound'});
    rows = [rows;row]; fprintf('自适应误差检查 delta=%.0e x=%.9g\n',delta,row.x);
end
writetable(rows,fullfile(folder,[name,'.csv'])); save(fullfile(folder,[name,'.mat']),'rows','tol','step');
    function values = value(t)
        ref = pp_reference(t,p,delta); values = reshape(ref.bias.^2,size(t));
    end
    function values = actual(t)
        ref = pp_reference(t,p,delta); y = deval(r.sols{3},t).';
        values = reshape(sum((y(:,1:2)-ref.x).^2,2),size(t));
    end
end
