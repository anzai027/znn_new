function [dg,a] = m8a_bdf2_system(t,g,problem,c)
% 调用二阶后向差分的预测纠错模型。
if nargin < 4, c = bd_config(); end
c.order = 2;
[dg,a] = bd_system(t,g,problem,c,"M8a-BDF2");
end
