function [g,dg,h,dh]=hmw7gx(x)
x1=x(1);
x2=x(2);

g1=-2-x1;
g2=x1-.5;
g3=-.5-x2;
g4=x2-4.5;
g=[g1 g2 g3 g4];
dgdx1=[-1 1 0 0];
dgdx2=[0 0 -1 1];
dg=[dgdx1;dgdx2];
h=[];
dh=[];