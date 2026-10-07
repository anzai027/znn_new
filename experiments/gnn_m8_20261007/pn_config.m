function c = pn_config(model)
% 保持旧实验条件并分别设置两个阻尼。
c = pp_config();
c.kappa = 0.01;
c.correct = "fixed";
c.predict = 1e-10;
if nargin > 0 && upper(string(model)) == "M8A", c.correct = "adaptive"; end
end
