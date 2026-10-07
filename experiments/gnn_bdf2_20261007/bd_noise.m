function [n,dn] = bd_noise(t,c)
% 用固定参数生成可重复的残差测量噪声。
n = zeros(7,1); dn = n;
if c.amp == 0, return, end
if c.kind == "constant"
    n = c.amp*ones(7,1);
else
    a = c.wave.a; b = c.wave.b; w = c.wave.w;
    n = c.amp*(a*cos(w*t)+b*sin(w*t));
    dn = c.amp*(-a*(w.*sin(w*t))+b*(w.*cos(w*t)));
end
end
