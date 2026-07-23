function dY = my_system(t,X)
x = X(1);
y = X(2);
dx = x*sin(t)+sin(y);
dy = x*sin(t)*cos(y)+cos(y);
dY=[dx;dy];

end