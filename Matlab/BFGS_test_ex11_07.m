function BFGS_test_ex11_07()
%close all
fx_fun='hmw7fx'
gx_fun='hmw7gx'
DVn=2;
delta=.000000000005;
xk=[-1 3.5];
Hk=eye(DVn);
epsilon=0.0001;
k=0;
Rk=2;
[fk,dfk]=feval(fx_fun,xk);
ck=dfk;
iterN=0;
iter_x_history=[xk];
iter_y_history=[fk];
figure,hold on




%[x1,x2]=meshgrid(0:2:100,0:2:100);
[x1,x2]=meshgrid(-3.5:.1:4.5,-3.5:0.1:4.5);
for i=1:length(x1(:,1))
    for j=1:length(x2(:,1))
        xi=[x1(i,j) x2(i,j)];
        [fk,dfk]=feval(fx_fun,xi);
        [gk,dgk]=feval(gx_fun,xi);
        f(i,j)=fk;
        g1(i,j)=gk(1);
        g2(i,j)=gk(2);
        g3(i,j)=gk(3);
        g4(i,j)=gk(4);
    end
end

[const,h]=contour(x1,x2,f,100);
var=max(g1)
cv1=[0.001 0.0015]*max(var)

[const1,h]=contour(x1,x2,g1,cv1,'k','LineWidth',2);

cv1=[0.001 0.0015]*max(max(g2));
[const,h]=contour(x1,x2,g2,cv1,'k','LineWidth',2);
cv1=[0.001 0.0015]*max(max(g3));
[const,h]=contour(x1,x2,g3,cv1,'k','LineWidth',2);
cv1=[0.001 0.0015]*max(max(g4));
[const,h]=contour(x1,x2,g4,cv1,'k','LineWidth',2);


while 1
    if norm(ck)<epsilon|iterN>100
        break;
    end
    
    [fk,dfk]=feval(fx_fun,xk);
    [gk,dgk,hk,dhk]=feval(gx_fun,xk);
    
    
    if k>=1
    sk=ak*dk;
    yk=cknew-ck;
    Dk=(yk*yk')/(yk'*sk);
    Ek=(ck*ck')/(ck'*dk);
    Hknew=Hk+Dk+Ek;
    H=Hknew

    end
    
    
    
    ci=dfk;
    A=dgk';
    b=-gk;
    Aeq=dhk';
    beq=-hk;
    H=eye(DVn);
    [dxp,fval,exitflag,output,lambda]=quadprog(H,ci,A,b,Aeq,beq);
    dk=dxp;
    
    rk=sum(lambda.ineqlin)+sum(abs(lambda.eqlin));
    Vk=max([0,gk]);
    
    if Vk>0
        Rk=2;
    end
    
    Rk=max(Rk,rk);
    
    
    
%     if norm(dk)>1
%          dk=dk/norm(dk);
%     end
   
    
    [ak,fak,x_history,y_history]=golden_section_with_x_mod01(fx_fun,gx_fun,xk,dk,Rk,delta);
    scatter(xk(1),xk(2),15,'MarkerEdgeColor','b','MarkerFaceColor','b')
    
     
    %dk=inv(Hk)*-ck;
    xknew=xk+ak*dk;
    [fknew,dfknew]=feval(fx_fun,xknew);
    cknew=dfknew;
    
    %%%%hk update
    
   
    if k>1
    xk=xknew;
    ck=cknew;
    Hk=Hknew;
    fk=fknew;
    end
    iterN=iterN+1;
    iter_x_history=[iter_x_history; xk];
    iter_y_history=[iter_y_history; fk];
    k=k+1
end


scatter(xk(1),xk(2),20,'MarkerEdgeColor','r','MarkerFaceColor','r')

% xk;
% 
% %[x1,x2]=meshgrid(-40:0.1:40,-40:0.1:40);
% 
% for i=1:length(x1(:,1))
%     for j=1:length(x2(:,1))
%         xi=[x1(i,j) x2(i,j)];
%         [fk,dfk]=feval(fx_fun,xi);
%         f(i,j)=fk;
%     end
% end
% [const1,h]=contour(x1,x2,f,100)
% 
%     
    