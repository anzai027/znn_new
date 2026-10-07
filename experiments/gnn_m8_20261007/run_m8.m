function run_m8(name)
% 先验证两个模型再做五个预测项小对照。
folder = fileparts(mfilename('fullpath')); root = fileparts(fileparts(folder));
old = fullfile(root,'experiments','gnn_predictor_20261006');
addpath(root,'-begin'); addpath(old,'-begin'); addpath(folder,'-begin');
if nargin == 0, name = 'data'; end
out = fullfile(folder,name);
assert(~isfile(fullfile(out,'all.mat')),'已有结果，请保留原数据后使用新目录。');
if ~isfolder(out), mkdir(out); end
diary(fullfile(out,'run.log')); clean = onCleanup(@() diary('off'));
c = pn_config(); p = example2_problem();
files = ["pn_system.m","pn_config.m","pn_run.m","pn_metrics.m", ...
    "m8_system.m","m8a_system.m","run_m8.m"];
paths = "experiments/gnn_m8_20261007/"+files;
paths = [paths,"experiments/gnn_predictor_20261006/pp_parts.m", ...
    "experiments/gnn_predictor_20261006/pp_config.m", ...
    "experiments/gnn_predictor_20261006/pp_reference.m","example2_problem.m"];
codes = strings(size(paths));
for k = 1:numel(paths), codes(k) = digest(fullfile(root,paths(k))); end
hash = table(paths(:),codes(:),'VariableNames',{'file','code'});
source = load(fullfile(old,'data','reference.mat'),'ref','c');
assert(isequal(c.t,source.c.t) && c.delta == source.c.delta && isequal(c.g0,source.c.g0));
ref = source.ref;
save(fullfile(out,'reference.mat'),'ref','c','hash','-v7.3');
fprintf('MATLAB %s，输出 %d 个时间点\n',version,numel(c.t));
labels = ["M8","M8a","M8","M8a","M8a fixed", ...
    "M8 lp1e-9","M8 lp1e-8","M8 h1e-6","M8a h1e-6"];
models = ["M8","M8A","M8","M8A","M8A","M8","M8","M8","M8A"];
results = cell(numel(labels),1); stats = table();
for k = 1:numel(labels)
    ci = pn_config(models(k)); group = "main";
    if k > 2, ci.rtol = 1e-9; ci.atol = 1e-11; group = "precision"; end
    if k > 4, group = "ablation"; end
    if k == 5, ci.correct = "fixed"; end
    if k == 6, ci.predict = 1e-9; end
    if k == 7, ci.predict = 1e-8; end
    if k >= 8, ci.h = 1e-6; end
    fprintf('开始 %d/%d %s %s\n',k,numel(labels),group,labels(k));
    r = pn_run(models(k),p,ci); r.label = labels(k); r.group = group;
    if ~r.complete
        save(fullfile(out,sprintf('failed_%02d.mat',k)),'r','hash','-v7.3');
        error('M8:Incomplete','%s 未完成，原因 %s',labels(k),r.reason);
    end
    [r.metrics,row] = pn_metrics(r,ref,p);
    results{k} = r; stats = [stats;row];
    save(fullfile(out,sprintf('run_%02d.mat',k)),'r','hash','-v7.3');
    writetable(stats,fullfile(out,'metrics.csv'));
    fprintf('完成 %s %s calls=%d xd=%.9g gpeak=%.9g residual=%.9g\n', ...
        group,r.label,r.calls,row.xd,row.gpeak,row.residual);
end
save(fullfile(out,'all.mat'),'results','stats','ref','c','hash','-v7.3');
comparison = baseline(root,old,stats,name);
writetable(comparison,fullfile(out,'comparison.csv'));
save(fullfile(out,'comparison.mat'),'comparison');
disp(stats); disp(comparison);
end

function rows = baseline(root,old,stats,name)
% 复用已核验的旧轨迹并保留来源和精度。
paths = ["data/run_05.mat","data/run_06.mat","data/run_07.mat", ...
    "extra/run_01.mat","data/run_08.mat","vfcr/run_02.mat"];
labels = ["M3","M5","M6 k0.05","M6 k0.01","M7","VFCR-S"];
rows = table();
for k = 1:numel(paths)
    data = load(fullfile(old,paths(k))); r = data.r;
    fields = fieldnames(data.hash);
    for j = 1:numel(fields)
        field = fields{j};
        if strcmp(field,'problem'), path = fullfile(root,'example2_problem.m');
        elseif strcmp(field,'my_system'), path = fullfile(root,'my_system.m');
        elseif strcmp(field,'fsmooth'), path = fullfile(root,'experiments','gnn_v2_20261002','fsmooth.m');
        elseif strcmp(field,'driver'), path = fullfile(old,'pp_vfcr.m');
        else, path = fullfile(old,[field,'.m']); end
        assert(strcmp(digest(path),data.hash.(field)),'旧模型源码改变：%s',field);
    end
    assert(r.complete && r.c.rtol == 1e-9 && r.c.atol == 1e-11);
    source = load(fullfile(old,'data','reference.mat'),'ref','c'); ref = source.ref;
    assert(isequal(r.t,ref.t) && isequal(r.c.g0,source.c.g0) && r.c.delta == source.c.delta);
    m = r.metrics; take = r.t >= r.c.tail; at = r.t(take);
    weight = @(v) sqrt(trapz(at,v(take).^2)/(at(end)-at(1)));
    file = "experiments/gnn_predictor_20261006/"+paths(k);
    row = table("old","precision",labels(k),r.c.rtol,r.c.atol,r.c.delta,r.calls,r.seconds, ...
        weight(m.xd),weight(m.x),weight(m.g),max(m.g(take)),max(m.residual(take)), ...
        'VariableNames',{'origin','group','label','rtol','atol','delta','calls','seconds', ...
        'xd','x','g','gpeak','residual'});
    row.file = file; rows = [rows;row];
end
for k = 1:height(stats)
    if stats.group(k) == "main", continue, end
    s = stats(k,:);
    row = table("new",s.group,s.label,s.rtol,s.atol,1e-4,s.calls,s.seconds, ...
        s.xd,s.x,s.g,s.gpeak,s.residual,'VariableNames',rows.Properties.VariableNames(1:13));
    row.file = "experiments/gnn_m8_20261007/"+string(name)+"/"+sprintf('run_%02d.mat',k);
    rows = [rows;row];
end
end

function code = digest(path)
% 记录实际使用源码的校验值。
file = fopen(path,'rb'); assert(file >= 0); clean = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8'); md = java.security.MessageDigest.getInstance('SHA-256');
md.update(bytes); code = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
