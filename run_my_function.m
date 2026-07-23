clear all
close all
clc
x0 = [5;10];
[t,x] = ode15s('my_system',[0:0.001:5]',x0);

figure
plot(t,x(:,1));
xlabel('t')
ylabel('x')

figure
plot(t,x(:,2));
xlabel('t')
ylabel('y')
