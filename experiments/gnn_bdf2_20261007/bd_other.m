function bd_other()
% 独立复算复核组与联合组的保存指标。
folder = bd_setup(); p = example2_problem(); rows = table();
for group = ["extra","joint"]
    count = 8; if group=="joint", count=2; end
    for k=1:count
        s=load(fullfile(folder,group,sprintf('run_%02d.mat',k)),'r','row'); r=s.r;
        assert(r.complete && r.last==10);
        sref=load(fullfile(folder,'results',sprintf('ref_%g.mat',r.c.delta)),'ref'); ref=sref.ref;
        xd=vecnorm(r.g(:,1:2)-ref.xd,2,2); g=vecnorm(r.g-ref.g,2,2);
        take=r.t>=5; t=r.t(take); rms=sqrt(trapz(t,xd(take).^2)/5);
        assert(isequal(xd,r.metrics.xd) && isequal(g,r.metrics.g));
        error=abs(rms-s.row.xd); assert(error<=max(1e-20,s.row.xd*1e-12));
        fresh=zeros(size(r.y)); ends=[0,r.c.h,2*r.c.h,10];
        for m=1:3
            take=r.t>=ends(m)&r.t<=ends(m+1);
            if any(take), fresh(take,:)=deval(r.sols{m},r.t(take)).'; end
        end
        trajectory=max(abs(fresh-r.y),[],'all'); assert(trajectory==0);
        residual=0;
        for m=unique(round(linspace(1,numel(r.t),500)))
            xi=pp_parts(r.t(m),r.g(m,:).',p,r.c.delta);
            residual=max(residual,abs(norm(xi)-r.metrics.residual(m)));
        end
        assert(residual==0);
        row=table(group,k,rms,max(g(r.t>=5)),error,trajectory,residual, ...
            'VariableNames',{'group','run','xd','gpeak','metricerror','trajectory','residual'});
        rows=[rows;row];
    end
end
writetable(rows,fullfile(folder,'other_check.csv')); save(fullfile(folder,'other_check.mat'),'rows');
s=load(fullfile(folder,'results','plan.mat'),'paths','hash');
for k=1:numel(s.paths), assert(s.hash(k)==bd_hash(fullfile(folder,s.paths(k)))); end
files=dir(fullfile(folder,'*.m')); names=string({files.name}); codes=strings(size(names));
for k=1:numel(names), codes(k)=bd_hash(fullfile(folder,names(k))); end
extra=["pp_parts","pp_reference","pp_config","pn_config","fsmooth","example2_problem","my_system"];
for name=extra
    path=which(name); names(end+1)=string(path); codes(end+1)=bd_hash(path);
end
paths=table(names.',codes.','VariableNames',{'file','sha256'});
writetable(paths,fullfile(folder,'sources.csv')); runtime=version;
save(fullfile(folder,'sources.mat'),'paths','runtime'); disp(rows);
end
