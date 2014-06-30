clear all
close all
clc



a1 = 0;
b1 = pi;
N1 = 250;
h1=(b1-a1)/N1;
x1=a1:h1/2:b1;
y = sin(x1);
sol1=h1/6*( y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end))

x2=0:pi/500:pi;
y = sin(x2);
h2=2*(x2(end)-x2(1))/(length(x2)-1);
sol2=h2/6*( y(1) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)) + y(end))

yf = @(x) sin(x);
a = 0;
b = pi;
N = 250;
h=(b-a)/N;
x=a:h/2:b;
fx=feval(yf,x);
sol=h/6*( fx(1) + 4*sum(fx(2:2:end-1)) + 2*sum(fx(3:2:end-2)) + fx(end))

z = [(0:pi/500:pi)', sin(0:pi/500:pi)'];
simpson_mod(z)