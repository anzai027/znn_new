function folder = bd_setup()
% 将本轮及历史共用函数加入搜索路径。
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
addpath(root,'-begin');
addpath(fullfile(root,'experiments','gnn_v2_20261002'),'-begin');
addpath(fullfile(root,'experiments','gnn_predictor_20261006'),'-begin');
addpath(fullfile(root,'experiments','gnn_m8_20261007'),'-begin');
addpath(folder,'-begin');
end
