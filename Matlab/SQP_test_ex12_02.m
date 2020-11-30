function SQP_test_ex12_02()
%close all
fx_fun = 'hmw7fx'
gx_fun = 'hmw7gx'


figure, hold on

%[x1,x2]=meshgrid(0:2:100,0:2:100);
[x1,x2]=meshgrid(-2:.1:4.5,-2.0:0.1:4.5);
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

cv1=[0.001 0.0015]*max(max(g1));

[const1,h]=contour(x1,x2,g1,cv1,'k','LineWidth',2);

cv1=[0.001 0.0015]*max(max(g2));
[const,h]=contour(x1,x2,g2,cv1,'k','LineWidth',2);

cv1=[0.001 0.0015]*max(max(g3));
[const,h]=contour(x1,x2,g3,cv1,'k','LineWidth',2);

cv1=[0.001 0.0015]*max(max(g4));
[const,h]=contour(x1,x2,g4,cv1,'k','LineWidth',2);

DVn=2;
xk=[-1 3.5]
Rk=2;
delta=0.01; 
tol=0.001;

iter_x_history=[];
iter_y_history=[];
iter_g_history=[];
iterloop=0;
while 1 
    
    [fk,dfk]=feval(fx_fun,xk);
    [gk,dgk,hk,dhk]=feval(gx_fun,xk);
    
    iter_x_history=[iter_x_history;xk];
    iter_y_history=[iter_y_history;fk];
    iter_g_history=[iter_g_history;gk];
    scatter3(xk(1),xk(2),fk,50,'MarkerEdgeColor','b','MarkerFaceColor','b')
    
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
        Rk=2.0;
    end
    
    %Rk=max(Rk,rk);
    phi0=fk+Rk*Vk;
    
    [ak,fak,x_history,y_history]=SQP_ch12_golden_section_with_x_mod01(fx_fun,gx_fun,xk,dk,Rk,delta);
    
    xlp=xk+ak*dk;
    scatter3(xlp(1),xlp(2),fval,50,'MarkerEdgeColor','r','MarkerFaceColor','r')
    xk=xlp;
    iterloop=iterloop+1
    if norm(dk)<tol 
        break;
    end
    xk;
    
end



xk
    