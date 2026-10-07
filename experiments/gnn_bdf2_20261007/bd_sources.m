function bd_sources()
% 保存最终代码与共用依赖的校验值。
folder=bd_setup(); files=dir(fullfile(folder,'*.m'));
names=string({files.name}); codes=strings(size(names));
for k=1:numel(names), codes(k)=bd_hash(fullfile(folder,names(k))); end
for name=["pp_parts","pp_reference","pp_config","pn_config","fsmooth","example2_problem","my_system"]
    path=which(name); names(end+1)=string(path); codes(end+1)=bd_hash(path);
end
paths=table(names.',codes.','VariableNames',{'file','sha256'});
writetable(paths,fullfile(folder,'sources.csv')); runtime=version;
save(fullfile(folder,'sources.mat'),'paths','runtime');
end
