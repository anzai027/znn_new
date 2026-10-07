function [dg,a] = pn_system(t,g,problem,c,model)
% 分别计算纠错方向和固定状态的预测速度。
model = upper(string(model));
if ~isscalar(model) || ~any(model == ["M8","M8A"])
    error('PN:Model','模型必须为 M8 或 M8a。');
end
check(c.gamma,'gamma',false); check(c.eps,'eps',false);
check(c.lambda,'lambda',false); check(c.h,'h',false);
g = g(:);
[xi,J,blocks] = pp_parts(t,g,problem,c.delta);
n = numel(g);
sigma = min(svd(J));
lc = c.lambda;
mode = string(c.correct);
if mode == "adaptive"
    check(c.floor,'floor',false); check(c.kappa,'kappa',true);
    lc = max(c.floor,c.kappa*sigma^2);
elseif mode ~= "fixed"
    error('PN:Correction','纠错阻尼必须为 fixed 或 adaptive。');
end
if ~isnumeric(c.start) || ~isscalar(c.start) || ~isreal(c.start) || ...
        ~isfinite(c.start) || t < c.start
    error('PN:Time','时间不能早于 start。');
end
boundary = t < c.start+c.h;
if boundary
    future = pp_parts(t+c.h,g,problem,c.delta);
    ft = (future-xi)/c.h;
else
    prior = pp_parts(t-c.h,g,problem,c.delta);
    ft = (xi-prior)/c.h;
end
dc = [J;sqrt(lc)*eye(n)]\[xi;zeros(n,1)];
lp = 0;
if model == "M8"
    check(c.predict,'predict',false);
    lp = c.predict;
    dp = [J;sqrt(lp)*eye(n)]\[-ft;zeros(n,1)];
else
    dp = -J\ft;
end
feedback = -c.gamma*dc/hypot(norm(dc),c.eps);
dg = dp+feedback;
if any(~isfinite(dg)), error('PN:Nonfinite','导数包含非有限数。'); end
a = struct('model',model,'xi',xi,'J',J,'sigma',sigma, ...
    'lc',lc,'lp',lp,'dc',dc,'dp',dp,'ft',ft,'feedback',feedback, ...
    'boundary',boundary,'blocks',blocks,'balance',norm(J*dp+ft), ...
    'condition',rcond(J));
end

function check(value,name,zero)
% 检查参数是否为范围内的有限实数。
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 0 || (~zero && value == 0)
    error('PN:Parameter','参数 %s 的数值不正确。',name);
end
end
