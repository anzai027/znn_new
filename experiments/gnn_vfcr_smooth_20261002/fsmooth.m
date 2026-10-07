function y = fsmooth(x,v,eps)
% 在零点附近用三次曲线连接原激活函数。
y = v.a.*exp(abs(x).^v.q).*abs(x).^v.p.*sign(x);
take = abs(x) < eps;
u = x(take)/eps;
k = v.p+v.q*eps^v.q;
b = v.a*exp(eps^v.q)*eps^v.p;
y(take) = b*((3-k)*u+(k-1)*u.^3)/2;
end
