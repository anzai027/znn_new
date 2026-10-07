function [dy,a] = bd_system(t,y,p,c,model)
% 在同一当前状态上用历史残差预测移动速度。
g = y(1:7); [xi,J] = pp_parts(t,g,p,c.delta);
n = bd_noise(t,c); obs = xi+n;
ft = zeros(7,1); calls = 1; stage = 0;
if model == "VFCR-S"
    x = g(1:2); mu = g(4:7);
    w = p.v(t)-p.Q(t)*x; s = hypot(hypot(w,mu),sqrt(c.delta));
    alpha = (mu.^2+c.delta)./(s.*(s+w));
    ft = [p.dG(t)*x+p.dh(t)+p.dP(t).'*g(3)+p.dQ(t).'*mu; ...
        p.dP(t)*x-p.du(t);alpha.*(p.dv(t)-p.dQ(t)*x)];
elseif t >= c.start+c.h
    prior = pp_parts(t-c.h,g,p,c.delta)+bd_noise(t-c.h,c);
    ft = (obs-prior)/c.h; calls = 2; stage = 1;
    if c.order == 2 && t >= c.start+2*c.h
        earlier = pp_parts(t-2*c.h,g,p,c.delta)+bd_noise(t-2*c.h,c);
        ft = (3*(obs-prior)-(prior-earlier))/(2*c.h);
        calls = 3; stage = 2;
    end
end
sigma = min(svd(J)); lc = max(c.floor,c.kappa*sigma^2);
dc = [J;sqrt(lc)*eye(7)]\[obs;zeros(7,1)];
lp = 0;
if startsWith(model,"VFCR")
    v = struct('a',1,'p',0.5,'q',0.5);
    rho = exp(acot(t)+1); z = y(8:14);
    phi = fsmooth(obs,v,c.eps);
    rhs = -rho*phi-fsmooth(obs+z,v,c.eps);
    dy = [J\(-ft+rhs);rho*phi];
    dp = -J\ft; fb = J\rhs;
else
    if model == "M8-reg"
        lp = c.predict;
        dp = [J;sqrt(lp)*eye(7)]\[-ft;zeros(7,1)];
    else
        dp = -J\ft;
    end
    fb = -c.gamma*dc/hypot(norm(dc),c.eps);
    dy = dp+fb;
end
assert(all(isfinite(dy)),'动力学产生非有限数。');
a = struct('xi',xi,'obs',obs,'J',J,'ft',ft,'dp',dp,'dc',dc, ...
    'fb',fb,'sigma',sigma,'lc',lc,'lp',lp,'stage',stage,'calls',calls);
end
