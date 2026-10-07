function run_predictor(name)
% 先比较四个模型再统一提高积分精度。
folder = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(folder));
addpath(root,'-begin'); addpath(folder,'-begin');
if nargin == 0, name = 'data'; end
out = fullfile(folder,name);
assert(~isfile(fullfile(out,'all.mat')),'已有结果，请传入新的输出文件夹名称。');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'run.log'));
clean = onCleanup(@() diary('off'));
c = pp_config();
p = example2_problem();
fprintf('MATLAB %s\n',version);
fprintf('输出 %d 个时间点\n',numel(c.t));
files = dir(fullfile(folder,'*.m'));
hash = struct();
for k = 1:numel(files), hash.(erase(files(k).name,'.m')) = digest(fullfile(folder,files(k).name)); end
hash.problem = digest(fullfile(root,'example2_problem.m'));
ref = pp_reference(c.t,p,c.delta);
save(fullfile(out,'reference.mat'),'ref','c','hash','-v7.3');
fprintf('参考残差最大 %.9g，参考速度峰 %.9g\n',max(ref.error),max(ref.speed));
models = ["M3","M5","M6","M7"];
results = cell(8,1);
stats = table();
for k = 1:8
    ci = c;
    group = "main";
    if k > 4
        ci.rtol = 1e-9; ci.atol = 1e-11; group = "precision";
    end
    model = models(mod(k-1,4)+1);
    fprintf('开始 %d/8 %s %s\n',k,group,model);
    r = pp_run(model,p,ci);
    if ~r.complete
        row = failure(r);
        fprintf('未完成 %s %s：%s\n',model,group,r.reason);
    else
        [r.metrics,row] = pp_metrics(r,ref,p);
    end
    row = addvars(row,group,ci.rtol,ci.atol,'Before',1,'NewVariableNames',{'group','rtol','atol'});
    stats = [stats;row];
    results{k} = r;
    save(fullfile(out,sprintf('run_%02d.mat',k)),'r','hash','-v7.3');
    writetable(stats,fullfile(out,'metrics.csv'));
    fprintf('完成 %s %s calls=%d xd=%.9g x=%.9g gpeak=%.9g residual=%.9g\n', ...
        group,model,r.calls,row.xd,row.x,row.gpeak,row.residual);
end
save(fullfile(out,'all.mat'),'results','stats','ref','c','hash','-v7.3');
comparison = table();
for k = 1:4
    a = results{k}; b = results{k+4};
    gdiff = NaN; xdiff = NaN;
    if a.complete && b.complete
        gdiff = max(vecnorm(a.g-b.g,2,2));
        xdiff = max(vecnorm(a.g(:,1:2)-b.g(:,1:2),2,2));
    end
    row = table(a.model,gdiff,xdiff, ...
        abs(stats.xd(k)-stats.xd(k+4))/max(stats.xd(k+4),realmin), ...
        abs(stats.residual(k)-stats.residual(k+4))/max(stats.residual(k+4),realmin), ...
        'VariableNames',{'model','gdiff','xdiff','rmserel','residualrel'});
    comparison = [comparison;row];
end
writetable(comparison,fullfile(out,'precision.csv'));
disp(stats);
disp(comparison);
end

function row = failure(r)
% 未完成的组保留状态且不计算完整误差。
names = {'model','complete','calls','parts','seconds','xd','x','g', ...
    'residual','xdpeak','xpeak','gpeak','prediction','predictionpeak','oraclepeak', ...
    'speed','sigma','lambdamin','lambdamax','alpha','limited','equality','inequality','tc4','tc6'};
values = [{r.model,r.complete,r.calls,r.parts,r.seconds},num2cell(nan(1,numel(names)-5))];
row = cell2table(values,'VariableNames',names);
end

function code = digest(path)
% 记录运行时源码的校验值。
file = fopen(path,'rb'); clean = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8');
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes);
value = typecast(md.digest(),'uint8');
code = lower(reshape(dec2hex(value,2).',1,[]));
end
