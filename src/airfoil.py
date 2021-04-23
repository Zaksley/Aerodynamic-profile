## Modules

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

from load_foil import load_foil

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

# Foil interploation getter
#
# @param sys.argv[1] a .dat file with points coordinates
#
def get_foil_spline():
    if len(sys.argv) == 1:
        print( "Please provide an input file")
        exit(1)

    (dim,ex,ey,ix,iy) = load_foil( sys.argv[1] )

    # First derivative computation
    de1 = (ey[1]-ey[0])/(ex[1]-ex[0])
    de2 = (ey[-1]-ey[-2])/(ex[-1]-ex[-2])
    di1 = (iy[1]-iy[0])/(ix[1]-ix[0])
    di2 = (iy[-1]-iy[-2])/(ix[-1]-ix[-2])

    # Cubic spline on extrados & intrados
    extrados = spline(ex,ey,de1,de2)
    intrados = spline(ix,iy,di1,di2)

    return (extrados, intrados, ex, ey, ix, iy)


## Tests

# Test on an unkown function (airfoil)
#
def test__spline_foil():
    
    (extrados, intrados, ex, ey, ix, iy) = get_foil_spline()

    de1 = (ey[1]-ey[0])/(ex[1]-ex[0])
    de2 = (ey[-1]-ey[-2])/(ex[-1]-ex[-2])

    # Nb of points multiplicator
    N = 100

    # New points computation
    ex_new = np.linspace(ex[0],ex[-1],len(ex)*N)
    ix_new = np.linspace(ix[0],ix[-1],len(ix)*N)
    ey_new = [splint(ex,ey,extrados,k) for k in ex_new]
    iy_new = [splint(ix,iy,intrados,k) for k in ix_new]

    # Points given
    plt.scatter(ex,ey, c='b')
    plt.scatter(ix,iy, c='b')

    # First interpolation
    plt.plot(ex, ey, c='g')
    plt.plot(ix, iy, c='g')

    # Beter interpolations
    plt.plot(ex_new, ey_new, c='r')
    plt.plot(ix_new, iy_new, c='r')

    plt.show()

    return

# Test on a known function (x^3)
#
def test__spline_knownfunc():

    #Tests on x^3 :
    ex = [k for k in np.linspace(0,1,10)]
    ey = [(k**3) for k in ex]
    de1 = 0
    de2 = 3
    extrados = spline(ex,ey,de1,de2)
    better = [k for k in np.linspace(0,1,1000)]
    bettery = [k**3 for k in better]

    # Nb of points multiplicator
    N = 100

    # New points computation
    ex_new = np.linspace(ex[0],ex[-1],len(ex)*N)
    ey_new = [splint(ex,ey,extrados,k) for k in ex_new]

    plt.plot(better, [abs(bettery[k] - ey_new[k]) for k in range(len(better))])
    plt.yscale('log')

    plt.show()

    return


## Main

if __name__ == '__main__':
    test__spline_foil()
