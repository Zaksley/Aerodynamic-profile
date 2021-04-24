import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from math import floor
from matplotlib.colors import Normalize

from load_foil import *
from airfoil import *
from integration import *


# Calcul f lambda
# ==> f = f_down & h = h_min for lower airflow
# ==> f = f_up & h = h_max for upper airflow 
def f_lambda(f, h, l, x):
    return ((1-l) * f(x)) + l*3*h

def derivative(f, X):
         # Initial values
    x1 = X[1]
    x2 = X[10]
    x3 = X[20]
    x4 = X[30]

    y1 = f(x1)
    y2 = f(x2)
    y3 = f(x3)
    y4 = f(x4)

        #Solving AX^3
    A = np.array([[x1**3, x1**2, x1, 1],
                  [x2**3, x2**2, x2, 1],
                  [x3**3, x3**2, x3, 1],
                  [x4**3, x4**2, x4, 1]])
    B = np.array([y1, y2, y3, y4])
    C = np.linalg.solve(A, B)

    def func(x):
        return (3*C[0]*x*x + 2*C[1]*x + C[2])
    return func

def length_graph(f, X_down, integrate_func):
    fp = derivative(f, X_down)
    def func(x):
        f = sqrt(1+fp(x)*fp(x))
        return f 
    return integrate_func(func, 0, 1, 5, 10 ** -10)

def zero(x):
    return 0

def create_map(file): 

    # Precision in the calcul 
    size = 1000
    A = np.zeros( (size, size) )
    Xp = np.arange(0.0, 1.0, 0.001)

        # Initialization
    # ======
    (dim, X_up, Y_up, X_down, Y_down) = load_foil(file)
    Y_up = 1.5*Y_up
    Y_down = 1.5*Y_down

    f_up = splint_f(X_up, Y_up, 1.0, 1.0)
    f_down = splint_f(X_down, Y_down, 1.0, 1.0)

    h_max = np.amax(Y_up)
    h_min = np.amin(Y_down)
    l_max = length_graph(f_up, X_up, simpsonIntegration)
    l_min = length_graph(zero, X_down, simpsonIntegration)
    # ======


    # Upper surface
    for l in range(1, size+1):
        def get_func(x):
            return f_lambda(f_up, h_max, l/size, x)

        length = length_graph(get_func, X_up, simpsonIntegration)
        d = (length - l_min)/(l_max - l_min)

        for xi in Xp:
            A[size//2 - floor(get_func(xi)*size), floor(xi*size)] = d

    # Lower surface
    for l in range(1, size+1):
        def get_func(x):
            return f_lambda(f_down, h_min, l/size, x)

        length = length_graph(get_func, X_down, simpsonIntegration)
        d = (length - l_min)/(l_max - l_min)

        for xi in Xp:
            A[size//2 - floor(get_func(xi)*size), floor(xi*size)] = d

    return A


def plot_map(A):
    plt.imshow(A, interpolation='none', cmap=plt.hot(), norm=Normalize(0.0, A.max()) )
 
    rgn = np.arange(len(A))
    ticks = [i for i in range(0, len(A), 100)]
    plt.xticks( range(0, len(A), 100), ticks)
    plt.yticks( range(0, len(A), 100), ticks)
 
    fig = plt.gcf()
    ax = fig.add_axes([0.85, 0.5, 0.04, 0.4])
    plt.colorbar(cax=ax)
 
    plt.show()






if __name__ == "__main__":
    if len(sys.argv) == 1:
        print( "Please provide a input file")
        exit(1)

    A = create_map(sys.argv[1])
    plot_map(A)