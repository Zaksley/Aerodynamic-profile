import numpy as np
from load_foil import *
from airfoil import *
from matplotlib.pyplot import figure, show

def hMaxMin(ey,iy,dim):
    eyMax = 0
    iyMin = 0
    for i in range(int(dim[0])):
        if ey[i] > eyMax:
            eyMax = ey[i]
        if iy[i] < iyMin:
            iyMin = iy[i]
    return (eyMax,iyMin)



def nextCurve(airflow,F,dim,hm,lamb):
    for i in range (int(dim[0])):
        airflow[i] = ((1-lamb)*F[i] + lamb*3*hm)

    

def airflowModelling(file,nbUpperCurves,nbLowerCurves):
    (dim,ex,ey,ix,iy) = load_foil(file)
    upperOrdinate = get_foil_spline()[0]
    upperOrdinate = get_foil_spline()[1]
    upperF = []
    lowerF = []
    for i in range (int(dim[0])):
        upperF.append(splint(ex,ey,upperOrdinate,float(i)/dim[0]))
        lowerF.append(splint(ix,iy,upperOrdinate,float(i)/dim[0]))
    (hmax,hmin) = hMaxMin(upperF,lowerF,dim)
    upperAirflow = np.zeros((nbUpperCurves , int(dim[0])) )
    lowerAirflow = np.zeros((nbLowerCurves , int(dim[0])) )
    for i in range(nbUpperCurves):
        nextCurve(upperAirflow[i],upperF,dim,hmax,float(i)/nbUpperCurves)
    for i in range(nbLowerCurves):
        nextCurve(lowerAirflow[i],lowerF,dim,hmin,float(i)/nbLowerCurves)
        
    return(dim,ex,ix,lowerAirflow,upperAirflow,upperF,lowerF)
 

def plot_curves(file,nbUpperCurves,nbLowerCurves):
    plt.axis('equal' )
    (dim,ex,ix,lowerAirflow,upperAirflow,F,res_down) = airflowModelling(file,nbUpperCurves,nbLowerCurves)
    for i in range(nbUpperCurves):
        plt.plot(ex,upperAirflow[i], '-k')
    for i in range(nbLowerCurves):
        plt.plot(ix,lowerAirflow[i], '-k')
    
    plt.plot(ex,F,'-k')
    plt.plot(ix,res_down,'-k')
    plt.show()
    



if __name__ == "__main__":
    if len(sys.argv) == 1:
        print( "Please provide a input file")
        exit(1)
    plot_curves( sys.argv[1], 20, 7 )
