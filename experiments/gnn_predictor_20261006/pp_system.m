function [dg,a] = pp_system(t,g,problem,c,model)
% 统一计算预条件反馈与固定状态时间预测。
model = upper(string(model));
if ~isscalar(model) || ~any(model == ["M3","M5","M6","M7"])
    error('PP:Model','model 必须为 M3、M5、M6 或 M7。');
end
check(c.gamma,'gamma',false);
check(c.eps,'eps',false);
check(c.lambda,'lambda',false);
g = g(:);
[xi,J,blocks] = pp_parts(t,g,problem,c.delta);
l = numel(g);
sigma = min(svd(J));
scale = ones(l,1);
if model == "M7"
    if ~isnumeric(c.scale) || numel(c.scale) ~= 3 || ~isreal(c.scale) || ...
            any(~isfinite(c.scale(:))) || any(c.scale(:) <= 0)
        error('PP:Scale','scale 必须包含三个有限正数。');
    end
    scale = [c.scale(1)*ones(blocks(1),1); ...
        c.scale(2)*ones(blocks(2),1);c.scale(3)*ones(blocks(3),1)];
end
D = diag(scale);
B = J*D;
if model == "M7"
    sigmab = min(svd(B));
else
    sigmab = sigma;
end
lambda = c.lambda;
if any(model == ["M6","M7"])
    check(c.floor,'floor',false);
    check(c.kappa,'kappa',true);
    lambda = max(c.floor,c.kappa*sigmab^2);
end
ft = zeros(l,1);
boundary = false;
if model ~= "M3"
    check(c.h,'h',false);
    if ~isnumeric(c.start) || ~isscalar(c.start) || ...
            ~isreal(c.start) || ~isfinite(c.start) || t < c.start
        error('PP:Time','时间不能早于 start。');
    end
    if t >= c.start+c.h
        prior = pp_parts(t-c.h,g,problem,c.delta);
        ft = (xi-prior)/c.h;
    else
        future = pp_parts(t+c.h,g,problem,c.delta);
        ft = (future-xi)/c.h;
        boundary = true;
    end
end
matrix = [B;sqrt(lambda)*eye(l)];
rhs = [xi,-ft;zeros(l,2)];
steps = matrix\rhs;
dc = steps(:,1);
dp = steps(:,2);
raw = D*(dp-c.gamma*dc/hypot(norm(dc),c.eps));
alpha = 1;
speed = [norm(raw(1:blocks(1))), ...
    norm(raw(blocks(1)+1:blocks(1)+blocks(2))), ...
    norm(raw(blocks(1)+blocks(2)+1:end))];
if model == "M7"
    if ~isnumeric(c.limit) || numel(c.limit) ~= 3 || ~isreal(c.limit) || ...
            any(~isfinite(c.limit(:))) || any(c.limit(:) <= 0)
        error('PP:Limit','limit 必须包含三个有限正数。');
    end
    alpha = min([1,c.limit(:).'./speed]);
end
dg = alpha*raw;
if any(~isfinite(dg))
    error('PP:Nonfinite','模型产生了非有限状态导数。');
end
a = struct('model',model,'xi',xi,'J',J,'sigma',sigma, ...
    'lambda',lambda,'dc',D*dc,'dp',D*dp,'raw',raw, ...
    'alpha',alpha,'boundary',boundary,'ft',ft,'B',B, ...
    'sigmab',sigmab,'scale',scale,'blocks',blocks,'speed',speed);
end

function check(value,name,zero)
% 检查标量参数的范围。
if ~isnumeric(value) || ~isscalar(value) || ~isreal(value) || ...
        ~isfinite(value) || value < 0 || (~zero && value == 0)
    error('PP:Parameter','参数 %s 的数值不正确。',name);
end
end
