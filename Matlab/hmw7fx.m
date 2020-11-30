function [y,df]=hmw7fx(x)

x1=x(1);
x2=x(2);
y=10*x1^4-20*x1^2*x2+10*x2^2+x1^2-2*x1+5;

dydx1=40*x1^3 - 40*x1*x2 + 2*x1 - 22;
dydx2=-20*x1^2 + 20*x2;
df=[dydx1;dydx2];