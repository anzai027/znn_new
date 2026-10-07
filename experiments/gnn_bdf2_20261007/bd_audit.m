function bd_audit()
% 独立核对差分、因果性和 VFCR 的等价实现。
folder = bd_setup(); p = example2_problem(); c = bd_config();
s = RandStream('mt19937ar','Seed',17); errors = zeros(100,4);
for k = 1:100
    t = 0.1+9.8*rand(s); g = randn(s,7,1); y = [g;randn(s,7,1)];
    [dy,a] = bd_system(t,g,p,c,"M8a-BDF2");
    xi = pp_parts(t,g,p,c.delta); z1 = pp_parts(t-c.h,g,p,c.delta);
    z2 = pp_parts(t-2*c.h,g,p,c.delta);
    ft = (3*xi-4*z1+z2)/(2*c.h);
    dc = (a.J.'*a.J+a.lc*eye(7))\(a.J.'*xi);
    expected = -a.J\a.ft-c.gamma*dc/hypot(norm(dc),c.eps);
    errors(k,1) = norm(ft-a.ft)/max(1,norm(ft));
    errors(k,2) = norm(dy-expected)/max(1,norm(dy));
    v = struct('a',1,'p',0.5,'q',0.5); q = p; q.phi = @(x) fsmooth(x,v,c.eps);
    actual = my_system(t,y,7,q,1,1,1,1,1,0.5,0.5,c.delta);
    ours = bd_system(t,y,p,c,"VFCR-S");
    errors(k,3) = norm(actual-ours)/max(1,norm(actual));
    h = 1e-4; exact = -Jtime(t,g,p,c.delta,a.J);
    ci = c; ci.h = h;
    [~,b] = bd_system(t,g,p,ci,"M8a-BDF2");
    ci.h = h/2; [~,d] = bd_system(t,g,p,ci,"M8a-BDF2");
    errors(k,4) = norm(b.dp-exact)/max(norm(d.dp-exact),realmin);
end
assert(max(errors(:,1))<1e-8 && max(errors(:,2))<1e-8 && max(errors(:,3))<1e-10);
q = rmfield(p,{'dG','dh','dP','du','dQ','dv'});
for t = [0,c.h/2,c.h,1.5*c.h,2*c.h,3]
    [v,a] = bd_system(t,c.g0,q,c,"M8a-BDF2");
    assert(a.stage == (t>=c.h)+(t>=2*c.h));
    assert(isequal(v,bd_system(t,c.g0,q,c,"M8a-BDF2")));
end
ci = c; ci.amp = 1e-5; ci.kind = "colored";
ci.wave.a = randn(s,7,16)/4; ci.wave.b = randn(s,7,16)/4;
t = 2.5; [n,dn] = bd_noise(t,ci);
fd = (bd_noise(t+1e-7,ci)-bd_noise(t-1e-7,ci))/2e-7;
assert(norm(dn-fd)/norm(dn)<1e-7 && isequal(n,bd_noise(t,ci)));
data = table(max(errors(:,1)),max(errors(:,2)),max(errors(:,3)), ...
    median(errors(:,4)),min(errors(:,4)),max(errors(:,4)), ...
    'VariableNames',{'bdf','formula','vfcr','ratio','ratiomin','ratiomax'});
writetable(data,fullfile(folder,'audit.csv')); save(fullfile(folder,'audit.mat'),'data','errors');
disp(data);
end
function ft = Jtime(t,g,p,delta,J)
% 用解析导数独立核对预测方向。
x = g(1:2); mu = g(4:7); w = p.v(t)-p.Q(t)*x;
s = hypot(hypot(w,mu),sqrt(delta)); alpha = 1-w./s;
ft = J\[p.dG(t)*x+p.dh(t)+p.dP(t).'*g(3)+p.dQ(t).'*mu; ...
    p.dP(t)*x-p.du(t);alpha.*(p.dv(t)-p.dQ(t)*x)];
end
