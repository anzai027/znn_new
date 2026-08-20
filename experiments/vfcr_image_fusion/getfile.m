function path = getfile(root,name)
item = dir(fullfile(root,'**',name));
if isempty(item)
    error('VFCR:Data','Missing dataset file: %s',name)
end
path = fullfile(item(1).folder,item(1).name);
end
