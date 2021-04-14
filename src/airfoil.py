## Modules

import numpy as np
import matplotlib.pyplot as plt
import math

import load_foil

### Airfoil refinement

## Algorithms

# Cubic spline algorithm
#
# @param x array of x-coordinates
# @param y array of y-coordinates
# @param yp1 first derivative at point 1
# @param ypn first derivative at point n
#
def spline(x, y, yp1, ypn):
    n = len(x)
    y2 = [0 for k in range (n)]
    u = [k for k in range(n-1)]

    # Lower boundary initialization
    if(yp1 > 0.99e30):
        y2[0] = 0
        u[0] = 0
    else:
        y2[0] = -0.5
        u[0] = (3./x[1] - x[0]) * ((y[1] - y[0])/(x[1] - x[0]) - yp1)

    # Main loop
    for i in range(1, n-1):
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        p = sig * y2[i-1] + 2.
        y2[i] = (sig-1.) / p
        u[i] = (y[i+1] - y [i])/(x[i+1] - x[i]) - (y[i] - y [i-1])/(x[i] - x[i-1])
        u[i] = (6. * u[i]/(x[i+1] - x[i-1]) - sig*u[i-1]) / p

    # Upper boundary set
    if(ypn > 0.99e30):
        qn = 0.
        un = 0.
    else:
        qn = 0.5
        un = (3./(x[n-1] - x[n-2])) * (ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]))

    y2[n-1] = (un - qn*u[n-2]) / (qn*y2[n-2] + 1.)

    # Backsubstitution loop
    for k in range(n-2, -1, -1):
        y2[k] = y2[k]*y2[k+1] + u[k]

    return y2


# Splint Interpolation algorithm
#
# @param xa array of x-coordinates
# @param ya array of y-coordinates
# @param y2a result of spline algorithm on points xa;ya
# @param x value at which to compute f(x)
#
def splint(xa, ya, y2a, x):
    n = len(xa)
    klo = 0
    khi = n-1

    #Main dichotomy loop to find in which section x is
    while(khi - klo > 1):
        k = (khi + klo)//2
        print(k)
        if(xa[k] > x):
            khi = k
        else:
            klo = k

    h = xa[khi] - xa[klo]
    # Bad input (Duplicate xa values)
    assert(h != 0)

    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    den = ((a**3 - a)*y2a[klo] + (b**3 - b)*y2a[khi])*h*h / 6.

    return a*ya[klo] + b*ya[khi] + den



## Tests

# Montrer des courbes : ce qui marche et ce qui marche pas (histoire de voir ce qui fonctionne)

def test__spline():
    #Todo
    return


## Main

if __name__ == '__main__':
    print(test__spline())





















