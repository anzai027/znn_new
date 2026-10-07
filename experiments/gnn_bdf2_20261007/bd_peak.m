function bd_peak()
% 在检测到的窄误差尖峰附近继续细化。
folder=bd_setup(); p=example2_problem(); cached=load(fullfile(folder,'reference_fine.mat'),'t','smooth');
ids=[3,8,11,12,13,15,17]; rows=table();
for k=ids
    s=load(fullfile(folder,'results',sprintf('run_%02d.mat',k)),'r'); r=s.r;
    d=find([1e-4,1e-5,1e-6]==r.c.delta); ref=cached.smooth{d};
    take=cached.t>=5; t=cached.t(take); y=deval(r.sols{3},t).';
    error=vecnorm(y(:,1:7)-ref.g(take,:),2,2); [coarse,at]=max(error); center=t(at);
    local=unique([center;(center-2e-6:5e-9:center+2e-6).']);
    local=local(local>=5 & local<=10); truth=pp_reference(local,p,r.c.delta);
    truth=bd_polish(truth,p,r.c.delta); y=deval(r.sols{3},local).';
    error=vecnorm(y(:,1:7)-truth.g,2,2); [fine,at]=max(error);
    row=table(k,r.c.delta,center,coarse,local(at),fine,abs(fine/coarse-1), ...
        'VariableNames',{'run','delta','center','coarse','peak','fine','change'});
    rows=[rows;row];
end
writetable(rows,fullfile(folder,'peak_check.csv')); save(fullfile(folder,'peak_check.mat'),'rows'); disp(rows);
end
