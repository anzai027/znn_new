function [dg,a] = m7_system(t,g,problem,c)
% 使用变量尺度和共同限速的预测模型。
[dg,a] = pp_system(t,g,problem,c,"M7");
end
