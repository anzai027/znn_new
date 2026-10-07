function [dg,a] = m5_system(t,g,problem,c)
% 使用固定阻尼的预测模型。
[dg,a] = pp_system(t,g,problem,c,"M5");
end
