function [a_opt,f_opt,x_history,y_history]=golden_section_with_x_mod01(fx_fun,gx_fun,xk,dk,Rk,delta)

x_history=[];
y_history=[];

tol=0.1;
r=(1+sqrt(5))/2;
ir=1/r;
%delta=.01;
a0=delta;
%a0=0.0000000001;
x0=xk+a0*dk;
%Rk=2;

[f0]=feval(fx_fun,x0);
[g0]=feval(gx_fun,x0);
Vt=max([0,g0]);
f0=f0+Rk*Vt;

a1=delta+delta*r;
x1=xk+a1*dk;

[f1]=feval(fx_fun,x1);
[g1]=feval(gx_fun,x1);
Vt=max([0,g1]);
f1=f1+Rk*Vt;

% [f2]=feval(fx_fun,x1);
% [g2]=feval(gx_fun,x1);
% Vt=max([0,g2]);
% f2=f2+Rk*Vt;
% 
% [f3]=feval(fx_fun,x1);
% [g3]=feval(gx_fun,x1);
% Vt=max([0,g3]);
% f3=f3+Rk*Vt;
% 
% [f4]=feval(fx_fun,x1);
% [g4]=feval(gx_fun,x1);
% Vt=max([0,g4]);
% f4=f4+Rk*Vt;
% 
% f1=max(f1,f2,f3,f4)

id=2;

if f1>f0 
    delta=-delta;
end

x_history=[x_history x0 x1];
y_history=[y_history f0 f1];

whileN=1;

while 1
    
    a2=delta*sum(r.^[0:id]);
   x2=xk+a2*dk;
    
    [f2]=feval(fx_fun,x2);
    [g2]=feval(gx_fun,x2);
    Vt=max([0,g2]);
    f2=f2+Rk*Vt;
    x_history=[x_history x2];
    y_history=[y_history f2];
    
    if f0>f1 & f1<f2  
        break;
    else
        
        id=id+1;
        a0=a1;
        a1=a2;
        f0=f1;
        f1=f2;
    end
    
    whileN=whileN+1;
end

aL=a0; aA=a1; aU=a2;
fL=f0; fA=f1; fU=f2;
Intv0=aU-aL;
aB=aL+ir*Intv0;
xB=xk+aB*dk;

[fB]=feval(fx_fun,xB);
[gB]=feval(gx_fun,xB);
Vt=max([0,gB]);
fB=fB+Rk*Vt;

x_history=[x_history xB];
y_history=[y_history fB];

while 1
    
    if fA<fB
    
        aL=aL; aU=aB; aB=aA;
        fL=fL; fU=fB; fB=fA;
        Intv1=aU-aL;
        aA=aL+(1-ir)*Intv1;
        xA=xk+aA*dk;
        [fA]=feval(fx_fun,xA);
        [gA]=feval(gx_fun,xA);
        Vt=max([0,gA]);
        fA=fA+Rk*Vt;
        x_history=[x_history xA];
        y_history=[y_history fA];
    else
        
        aL=aA;  aU=aU;  aA=aB;
        fL=fA;  fU=fU;  fA=fB;
        Intv1=aU-aL;
        aB=aL+ir*Intv1;
        xB=xk+aB*dk;
        [fB]=feval(fx_fun,xB);
        [gB]=feval(gx_fun,xB);
        Vt=max([0,gB]);
        fB=fB+Rk*Vt;
        x_history=[x_history xB];
        y_history=[y_history fB];
        
    end
    
    if abs(Intv1-Intv0)<tol
        break
    else
        Intv0=Intv1;
    end
end

a_opt=(aU+aL)/2;
f_opt=(fU+fL)/2;

