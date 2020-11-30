def golden_section_fx_gx(fx_fun, gx_fun, xk, dk, Rk, delta):
    # ####comment after initial debugging
    # from numpy import transpose
    # xk=[-1,3.5]
    # #print("xk:"+str(xk))
    # delta=.01
    # dk=[-1,-1]
    # Rk=2
    # #####

    from numpy import asarray, transpose, max, amax, matrix, append, sum, arange
    from math import sqrt
    x_history = []
    y_history = []
    tol = .0001
    r = (1 + sqrt(5)) / 2
    ir = 1 / r
    a0 = delta
    # a0=.000000001
    x0 = xk + a0 * asarray(dk)
    x0 = x0  # [:,1]
    #print("x0:" + str(x0))
    f0, df0 = fx_fun(x0)
    g0, gf0,gx3,gx4 = gx_fun(x0)
    #print("f0:" + str(f0))
    #print("g0:" + str(g0))
    Vt = max([0, amax(g0)])
    f0 = f0 + Rk * Vt
    a1 = delta + delta * r
    x1 = xk + a1 * asarray(dk)
    f1, df1 = fx_fun(x1)
    g1, dg1,gx3,gx4 = gx_fun(x1)
    #print("f1:" + str(f1))
    #print("g1:" + str(g1))
    #print("Vt:" + str(Vt))
    #print("g0:" + str(g0))
    #print("f0:" + str(f0))
    Vt = max([0, amax(g1)])
    f1 = f1 + Rk * Vt
    #print("f1:" + str(f1))
    id = 2

    if f1 > f0:
        delta = -delta

    x_history = append(append(x_history, x0), x1)
    y_history = append(append(y_history, f0), f1)
    # print("x_history"+str(x_history))
    # print("y_history"+str(y_history))

    whileN = 0

    while 1:
        a2 = delta * sum(r ** arange(0, id, 1))
        # print("a2:"+str(a2))
        x2 = xk + a2 * asarray(dk)
        # print("x2:"+str(x2))
        f2, df2 = fx_fun(x2)
        g2, dg2,gx3,gx4 = gx_fun(x2)
        Vt = max([0, amax(g2)])
        f2 = f2 + Rk * Vt
        # print("f2:"+str(f2))
        # print('g2:'+str(g2))
        # print("vt:"+str(Vt))
        # print("f2:"+str(f2))
        x_history = append(x_history, x2)
        y_history = append(y_history, f2)
        # print("x_history"+str(x_history))
        # print("y_history"+str(y_history))

        if f0 > f1 and f1 < f2 or whileN > 100:
            break
        else:
            id = id + 1
            a0 = a1
            a1 = a2
            f0 = f1
            f1 = f2
        #             print("a0:"+str(a0))
        #             print("a1:"+str(a1))
        #             print("f0:"+str(f0))
        #             print("f1:"+str(f1))
        whileN = whileN + 1

    # print("whileN:"+str(whileN))
    aL = a0
    aA = a1
    aU = a2
    fL = f0
    fA = f1
    fU = f2
    # print("aL:"+str(aL))
    # print("aA:"+str(aA))
    # print("aU:"+str(aU))
    # print("fL:"+str(fL))
    # print("fA:"+str(fA))
    # print("fU:"+str(fU))
    Intv0 = aU - aL
    aB = aL + ir * Intv0
    xB = xk + aB * asarray(dk)
    # print("xB:"+str(xB))

    fB, dfB = fx_fun(xB)
    gB, dgB,gx3,gx4 = gx_fun(xB)
    Vt = max([0, amax(gB)])
    fB = fB + Rk * Vt

    x_history = append(x_history, xB)
    y_history = append(y_history, fB)

    while 1:

        if fA < fB:
            aL = aL
            aU = aB
            aB = aA
            fL = fL
            fU = fB
            fB = fA
            Intv1 = aU - aL
            aA = aL + (1 - ir) * Intv1
            xA = xk + aA * asarray(dk)
            fA, dfA = fx_fun(xA)
            gA, dgA,dg3,dg4 = gx_fun(xA)
            Vt = max([0, amax(gA)])
            fA = fA + Rk * Vt
            #print("fA:" + str(fA))
            x_history = append(x_history, xA)
            y_history = append(y_history, fA)
        else:
            aL = aA
            aU = aU
            aA = aB
            fL = fA
            fU = fU
            fA = fB
            Intv1 = aU - aL
            aB = aL + ir * Intv1
            xB = xk + aB * asarray(dk)
            fB, dfB = fx_fun(xB)
            gB, dgB,dg3,dg4 = gx_fun(xB)
            Vt = max([0, amax(gB)])
            fB = fB + Rk * Vt
            #print("fB:" + str(fB))
            x_history = append(x_history, xB)
            y_history = append(y_history, fB)
            if abs(Intv1 - Intv0) < tol:
                break
            else:
                Intv0 = Intv1
            # print("x_history:"+str(x_history))
    a_opt = (aU + aL) / 2
    f_opt = (fU + fL) / 2
    # print(y_history)
    # print("fU:"+str(fU))
    # print("a_opt:"+str(a_opt))
    # print("f_opt:"+str(f_opt))
    return a_opt, f_opt, x_history, y_history  # y history is the objective function value history


from gx_fun import *
from fx_fun import *

# def fx_fun(xi):
#     from numpy import asarray
#     x1=xi[0]
#     x2=xi[1]
#     #y= -(2*xi[0]+3*xi[1]-xi[0]**3-2*xi[1]**2);
#     #y=0.5*xi[0]**2+xi[1]**2-xi[0]*xi[1]-7*xi[0]-7*xi[1];
#     y=10*x1**4-(20*x1**2)*x2+10*x2**2+(x1**2)-2*x1+5
#     dydx1=40*x1**3 - 40*x1*x2 + 2*x1 - 22;
#     dydx2=-20*x1**2 + 20*x2;
#     df= asarray([dydx1,dydx2])
#     return y,df
# def gx_fun(x):
#     from numpy import transpose
#     x1=x[0]
#     x2=x[1]
#     g1=-2-x1
#     g2=x1-0.5
#     g3=-0.5-x2
#     g4=x2-4.5
#     g=[g1,g2,g3,g4]
#     dg1dx=[-1, 1, 0, 0]
#     dg2dx=[0, 0, -1, 1]
#     dg=transpose([dg1dx,dg2dx])
#     h=[]
#     dh=[]
#     return g,dg#,h,dh

#
# if __name__ == '__main__':
#     xk = [-1, 3.5]
#     dk = [-1, -1]
#     Rk = 2
#     delta = .01
#     a_opt,f_opt,x_history,y_history=golden_section_fx_gx(fx_fun, gx_fun, xk, dk, Rk, delta)
#     print("a_opt:"+str(a_opt))
#     print("f_opt:"+str(f_opt))
