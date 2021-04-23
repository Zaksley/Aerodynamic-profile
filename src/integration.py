import numpy as np

#integrations of f between start and end, the functions stop if the difference
#between two iteration if lower than epsilon or if the number of iteration
#reach iMax


def leftRectangleIntegration(f,start,end,iMax,epsilon):
    n = 1
    i = 1
    r1 = f(start)*(end-start)
    r0 = r1+epsilon+1 #an arbitrary value to enter in the first iteration
    while (i<iMax and abs(r1-r0) > epsilon):
        r0 = r1
        r1 = 0.5*r0
        for k in range(n):
            r1 += ((end-start)/(2*n))*f( start+(k+0.5)*(end-start)/n )
            n = n*2
            i += 1
    return r1

def midRectangleIntegration(f,start,end,iMax,epsilon):
    n = 1
    i = 1
    r1 = f((start+end)/2)*(end-start)
    r0 = r1+epsilon+1
    while (i<iMax and abs(r1-r0) > epsilon):
        r0 = r1
        r1 = (1/3)*r0
        for k in range(n):
            r1 += ((end-start)/(3*n))*f( start+(1/6+k)*(end-start)/n )
            r1 += ((end-start)/(3*n))*f( start+(5/6+k)*(end-start)/n )
            n = n*3
            i += 1
    return r1

def trapezeIntegration(f,start,end,iMax,epsilon):
    n = 1
    i = 1
    r1 = f((start+end)/2)*(end-start)
    r0 = r1+epsilon+1
    while (i<iMax and abs(r1-r0) > epsilon):
        r0 = r1
        r1 = 0.5*r0
        for k in range(n):
            r1 += ((end-start)/(2*n))*f( start+(k+0.5)*(end-start)/n )
            n = n*2
            i += 1
    return r1 

def simpsonIntegration(f,start,end,iMax,epsilon):
    I0 = (f(start)+f(end))/6
    I1 = 0
    I2 = 0
    n = 1
    i = 1
    r1 = f(start)*(end-start)
    r0 = r1+epsilon+1
    while (i<iMax and abs(r1-r0) > epsilon):
        r0 = r1
        I1 = I1+I2
        I2 = 0
        for k in range(2*n):
            I2 += f( start+(k/2+0.25)*(end-start)/n)
        r1 = (I0+(1/3)*I1+(2/3)*I2)*(end-start)/(2*n)
        n = n*2
        i += 1
    return r1

#calculation of the curve f between start and end
#iMax and epsilon are the variables describe in the integration methods
#integrate is the integration function that will be used


def curveLength(f,h,start,end,iMax,epsilon,integrate):
    derive = lambda x : (f(x+h)-f(x-h))/(2*h)
    fline = lambda x : np.sqrt(1+derive(x)**2)
    return integrate(fline,start,end,iMax,epsilon)