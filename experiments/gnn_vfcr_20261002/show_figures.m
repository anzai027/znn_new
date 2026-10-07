function show_figures(name)
% 打开指定图像，省略名称时打开全部图像。
folder = fileparts(mfilename('fullpath'));
if nargin == 0
    files = dir(fullfile(folder,'figures','*.fig'));
else
    files = dir(fullfile(folder,'figures',[char(name),'.fig']));
    assert(~isempty(files),'没有找到图像：%s',name);
end
for k = 1:numel(files)
    openfig(fullfile(files(k).folder,files(k).name),'new','visible');
end
end
