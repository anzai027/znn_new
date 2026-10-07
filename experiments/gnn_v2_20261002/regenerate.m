function regenerate()
% 只用保存的数据重建图像和报告。
folder = fileparts(mfilename('fullpath'));
data = load(fullfile(folder,'data','all.mat'));
addpath(data.meta.root,'-end');
assert(strcmpi(which('my_system'),fullfile(data.meta.root,'my_system.m')));
assert(strcmpi(which('gnn_system'),fullfile(data.meta.root,'gnn_system.m')));
conclude(data.results,data.refs,data.c,data.meta,folder);
end
