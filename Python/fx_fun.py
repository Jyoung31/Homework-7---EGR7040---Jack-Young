def fx_fun(xi):
    from numpy import asarray
    x1=xi[0]
    x2=xi[1]
    #y= -(2*xi[0]+3*xi[1]-xi[0]**3-2*xi[1]**2);
    #y=0.5*xi[0]**2+xi[1]**2-xi[0]*xi[1]-7*xi[0]-7*xi[1];
    y=10*x1**4-(20*x1**2)*x2+10*x2**2+(x1**2)-2*x1+5
    dydx1=40*x1**3 - 40*x1*x2 + 2*x1 - 22;
    dydx2=-20*x1**2 + 20*x2;
    df= asarray([dydx1,dydx2])
    return y,df