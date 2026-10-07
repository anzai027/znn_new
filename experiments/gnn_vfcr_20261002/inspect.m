function audit = inspect(folder,results,c)
% 重新打开数据和图像，检查保存结果。
files = dir(fullfile(folder,'figures','*.fig'));
audit.figures = strings(numel(files),1);
audit.sizes = zeros(numel(files),2);
audit.lines = zeros(numel(files),1);
audit.legends = zeros(numel(files),1);
audit.labels = cell(numel(files),1);
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);
    fig = openfig(file,'invisible');
    assert(isgraphics(fig),'FIG 无法打开。');
    lines = findall(fig,'Type','line');
    text = findall(fig,'Type','text');
    assert(~isempty(lines) || any(strlength(string({text.String})) > 0), ...
        'FIG 没有曲线或文字。');
    [~,name] = fileparts(file);
    labels = strings(0,1);
    if ~isempty(lines)
        labels = string({lines.DisplayName});
    end
    if any(string(name) == ["supplement_states","constraints"])
        assert(sum(contains(labels,"GNN gamma=100 p=0.5")) == 4, ...
            '补充参数图没有读取四条正确的新 GNN 曲线。');
    end
    axes = findall(fig,'Type','axes');
    for j = 1:numel(axes)
        if ~isempty(findall(axes(j),'Type','line'))
            assert(strlength(string(axes(j).XLabel.String)) > 0 ...
                && strlength(string(axes(j).YLabel.String)) > 0, ...
                '曲线图缺少坐标标签。');
        end
    end
    audit.lines(k) = numel(lines);
    audit.legends(k) = numel(findall(fig,'Type','legend'));
    audit.labels{k} = labels;
    close(fig);
    pic = imread(fullfile(folder,'figures',[name,'.png']));
    assert(~isempty(pic) && max(pic(:)) > min(pic(:)),'PNG 为空白。');
    audit.figures(k) = string(name);
    audit.sizes(k,:) = size(pic,[1,2]);
end
audit.runs = numel(results);
for k = 1:numel(results)
    r = results{k};
    file = fullfile(folder,'data',sprintf('run_%03d.mat',k));
    stored = load(file,'r');
    assert(stored.r.job.id == k,'保存的运行编号不一致。');
    assert(all(diff(r.raw.t) > 0),'接受点时间不是严格递增的。');
    assert(all(isfinite(r.raw.y),'all'),'接受点包含非有限状态。');
    assert(isequal(size(r.y,1),numel(c.t)),'输出网格长度不一致。');
    assert(isequaln(r.metrics.motion,motion(r,c)),'速度诊断与轨迹不一致。');
    if ~r.complete
        assert(all(isnan(r.metrics.tc)) && isnan(r.metrics.rmse), ...
            '失败运行不应拥有完整区间的性能指标。');
        assert(all(isnan(r.y(c.t > r.last+1e-10,:)),'all'), ...
            '中断后的未知状态不能被补造。');
    end
end
stats = collect(results);
for kind = unique(stats.noise(stats.group == "main")).'
    for seed = unique(stats.seed(stats.group == "main" & stats.noise == kind)).'
        ids = stats.id(stats.group == "main" & stats.noise == kind & stats.seed == seed);
        assert(numel(ids) == 2,'共同噪声组缺少模型。');
        assert(isequal(results{ids(1)}.bank,results{ids(2)}.bank), ...
            '两个模型的噪声样本不同。');
        assert(isequal(results{ids(1)}.job.g0,results{ids(2)}.job.g0), ...
            '两个模型的初值不同。');
    end
end
audit.noise = true;
audit.motion = sum(stats.tested & stats.motion == 0);
audit.passed = true;
save(fullfile(folder,'data','figures.mat'),'audit');
end
