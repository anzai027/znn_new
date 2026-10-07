function bd_figures()
% 用实际保存的数据绘制本轮比较图。
folder = bd_setup(); out = fullfile(folder,'figures');
if ~isfolder(out), mkdir(out); end
s = load(fullfile(folder,'results','all.mat'),'stats'); stats = s.stats;
colors = lines(5); ids = [1,3,6,7]; labels = ["M8a-FD1","M8a-BDF2","VFCR-FD2","VFCR-S"];
f = figure('Visible','off','Position',[50,50,1200,800]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for panel = 1:4
    ax = nexttile; hold(ax,'on');
    for k = 1:numel(ids)
        s = load(fullfile(folder,'results',sprintf('run_%02d.mat',ids(k))),'r'); r = s.r;
        take = r.t>=5;
        if panel > 2, take = abs(r.t-cast(9.74429445439319,'like',r.t))<=0.001; end
        if mod(panel,2)==1, v=r.metrics.xd; else, v=r.metrics.g; end
        plot(ax,r.t(take),max(v(take),1e-16),'LineWidth',1.3,'Color',colors(k,:), ...
            'DisplayName',labels(k));
    end
    set(ax,'YScale','log'); grid(ax,'on'); xlabel(ax,'Time (s)');
    if mod(panel,2)==1, ylabel(ax,'Primal tracking error'); else, ylabel(ax,'Full-state tracking error'); end
    legend(ax,'Location','best');
end
title(tl,'Zero noise: common tolerance 1e-11 / 1e-13, h = 1e-5'); finish(f,out,'zero');
f = figure('Visible','off','Position',[50,50,1200,800]); tl = tiledlayout(2,2,'TileSpacing','compact');
for panel = 1:3
    ax=nexttile; hold(ax,'on');
    for k = 1:3
        if panel==2 && k>1, continue, end
        model = ["M8a-BDF2","VFCR-FD","VFCR-S"]; take = stats.model==model(k) & stats.amp==0;
        if k==1, take=take & stats.h==1e-5; end
        if k==2, take=take & stats.order==2; end
        [d,ix]=sort(stats.delta(take),'descend'); part=stats(take,:); part=part(ix,:);
        fields=["xd","x","gpeak"]; loglog(ax,d,part.(fields(panel)),'o-','Color',colors(k,:), ...
            'LineWidth',1.5,'DisplayName',model(k));
        if panel==2
            delete(findobj(ax,'Type','line'));
            bias=readtable(fullfile(folder,'bias_tight.csv'));
            loglog(ax,bias.delta,bias.x,'o-','LineWidth',1.5,'Color',colors(1,:), ...
                'DisplayName','M8a-BDF2 (adaptive integral)');
        end
    end
    set(ax,'XScale','log','YScale','log'); xlabel(ax,'delta'); grid(ax,'on'); legend(ax,'Location','best');
    titles=["Primal tracking RMSE","Original QP RMSE","Sampled full-state peak error"]; ylabel(ax,titles(panel));
end
ax=nexttile; s=readtable(fullfile(folder,'reference_check.csv'));
yyaxis(ax,'left'); loglog(ax,s.delta,s.sigma,'o-','LineWidth',1.5); ylabel(ax,'Minimum singular value');
yyaxis(ax,'right'); loglog(ax,s.delta,s.speed,'s-','LineWidth',1.5); ylabel(ax,'Maximum root speed');
xlabel(ax,'delta'); grid(ax,'on'); title(tl,'Delta study: tracking, QP bias and conditioning'); finish(f,out,'delta');
f=figure('Visible','off','Position',[50,50,1200,800]); tl=tiledlayout(2,2,'TileSpacing','compact');
ids=[14,15,16,17]; labels=["M8a-FD1","M8a-BDF2","M8-reg","VFCR-FD2"];
for panel=1:4
    ax=nexttile; hold(ax,'on');
    for k=1:4
        s=load(fullfile(folder,'results',sprintf('run_%02d.mat',ids(k))),'r'); r=s.r;
        take=r.t>=5;
        if panel>=3, take=abs(r.t-9.74429445439319)<=0.003; end
        if mod(panel,2)==1, v=r.metrics.xd; else, v=r.metrics.residual; end
        plot(ax,r.t(take),max(v(take),1e-16),'Color',colors(k,:),'DisplayName',labels(k),'LineWidth',1);
    end
    set(ax,'YScale','log'); xlabel(ax,'Time (s)'); grid(ax,'on'); legend(ax,'Location','best');
    if mod(panel,2)==1, ylabel(ax,'Clean primal tracking error'); else, ylabel(ax,'Clean residual'); end
end
title(tl,'Residual measurement noise: 20-200 Hz, amplitude 1e-5, seed 1'); finish(f,out,'noise');
f=figure('Visible','off','Position',[50,50,1200,800]); tl=tiledlayout(2,2,'TileSpacing','compact');
for panel=1:2
    ax=nexttile; hold(ax,'on');
    for k=1:3
        models=["M8a-BDF2","M8-reg","VFCR-FD"];
        take=stats.model==models(k) & stats.amp==1e-5 & stats.seed==1 & stats.noise=="colored";
        part=stats(take,:); [h,ix]=sort(part.h); part=part(ix,:);
        fields=["xd","gpeak"]; loglog(ax,h,part.(fields(panel)),'o-','LineWidth',1.5, ...
            'Color',colors(k,:),'DisplayName',models(k));
    end
    set(ax,'XScale','log','YScale','log'); xlabel(ax,'h'); grid(ax,'on'); legend(ax,'Location','best');
    if panel==1, ylabel(ax,'Clean primal RMSE'); else, ylabel(ax,'Full-state peak error'); end
end
ax=nexttile; hold(ax,'on'); s=readtable(fullfile(folder,'iid.csv'));
loglog(ax,s.h,s.fd,'o-','DisplayName','FD1'); loglog(ax,s.h,s.bdf,'s-','DisplayName','BDF2');
set(ax,'XScale','log','YScale','log'); xlabel(ax,'h'); ylabel(ax,'Derivative noise RMS per component');
grid(ax,'on'); legend(ax,'Location','best'); title(ax,'Independent samples: algebraic test only');
ax=nexttile; hold(ax,'on');
for k=1:3
    models=["M8a-BDF2","M8-reg","VFCR-FD"];
    take=stats.model==models(k)&stats.group=="noise";
    plot(ax,stats.seed(take),stats.xd(take),'o-','Color',colors(k,:),'DisplayName',models(k),'LineWidth',1.5);
end
xlabel(ax,'Noise seed'); ylabel(ax,'Clean primal RMSE'); grid(ax,'on'); legend(ax,'Location','best');
title(tl,'Noise sensitivity: step size, independent-sample amplification and seeds'); finish(f,out,'sensitivity');
end
function finish(f,out,name)
% 保存图片和可编辑图并检查曲线数据。
savefig(f,fullfile(out,[name,'.fig'])); exportgraphics(f,fullfile(out,[name,'.png']),'Resolution',160);
lines1=findall(f,'Type','line'); values1=arrayfun(@(h) {h.XData,h.YData},lines1,'UniformOutput',false);
g=openfig(fullfile(out,[name,'.fig']),'invisible'); lines2=findall(g,'Type','line');
values2=arrayfun(@(h) {h.XData,h.YData},lines2,'UniformOutput',false);
assert(isequaln(values1,values2)); close(g); close(f);
end
