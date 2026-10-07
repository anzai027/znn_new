function [dg,a] = m8a_system(t,g,problem,c)
% 直接解预测方程并保留 M6 纠错。
if nargin < 4, c = pn_config("M8A"); end
[dg,a] = pn_system(t,g,problem,c,"M8A");
end
