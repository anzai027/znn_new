function [dg,a] = m8_system(t,g,problem,c)
% 使用独立的小预测阻尼。
if nargin < 4, c = pn_config("M8"); end
[dg,a] = pn_system(t,g,problem,c,"M8");
end
