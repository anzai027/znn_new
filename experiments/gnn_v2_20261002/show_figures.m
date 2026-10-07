function show_figures()
% 打开本次保存的全部 MATLAB 图像。
folder = fileparts(mfilename('fullpath'));
files = dir(fullfile(folder,'figures','*.fig'));
for k = 1:numel(files)
    openfig(fullfile(files(k).folder,files(k).name),'new','visible');
end
end
