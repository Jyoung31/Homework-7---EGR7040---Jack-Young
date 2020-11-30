def gx_fun(x):
    from numpy import transpose
    x1=x[0]
    x2=x[1]
    g1=-2-x1
    g2=x1-0.5
    g3=-0.5-x2
    g4=x2-4.5
    g=[g1,g2,g3,g4]
    dg1dx=[-1., 1., 0., 0.]
    dg2dx=[0., 0., -1., 1.]
    dg=transpose([dg1dx,dg2dx])
    h=[]
    dh=[]
    return g,dg,h,dh