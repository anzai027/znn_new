function [dg,a] = m6_system(t,g,problem,c)
% 使用自适应阻尼的预测模型。
[dg,a] = pp_system(t,g,problem,c,"M6");
end
