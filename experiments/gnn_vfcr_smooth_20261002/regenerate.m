function regenerate()
% 用已有数据重新生成报告和图像。
folder = fileparts(mfilename('fullpath'));
data = load(fullfile(folder,'data','all.mat'));
conclude(data.results,data.refs,data.c,data.meta,folder);
end
