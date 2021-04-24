import numpy as np
from load_foil import *
from airfoil import *
from matplotlib.pyplot import figure, show

def hMaxMin(ey,iy,dim):
    """
    Calculer h_max et h_min.

    paramètres:
    ey: liste des ordonnées des points de l'extrados
    iy: liste des ordonnées des points de l'intrados
    dim: nombre de points.

    retourne: 
    eyMax: h_max
    iyMin: h_min

    """
    eyMax = 0
    iyMin = 0
    for i in range(int(dim[0])):
        if ey[i] > eyMax:
            eyMax = ey[i]
        if iy[i] < iyMin:
            iyMin = iy[i]
    return (eyMax,iyMin)



def nextCurve(airflow,F,dim,hm,lamb):
    """
    Calculer toutes les lignes du flux d'air.

    paramètres:
    airflow: tableau vide qui va contenir les valeurs de la fonction f_lambda pour un lambda donné.
    F: tableau des valeurs de f.
    dim: nombre de points
    hm: h_max ou h_min.
    lambda: un réel € [0,1]

    """
    for i in range (int(dim[0])):
        airflow[i] = ((1-lamb)*F[i] + lamb*3*hm)

    

def airflowModelling(file,nbUpperCurves,nbLowerCurves):

    """
    Calculer les ordonnées de chaque ligne au dessus et en dessous de l'aile.

    paramètres:
    file: le fichier .dat qui contient les données.
    nbUpperCurves: nombre de lignes au dessus de l'aile.
    nbLowerCurves: nombre de lignes en dessous de l'aile.

    retourne:
    ex: tableau des valeurs des abscisses de l'extrados.
    ix: tableau des valeurs des abscisses de l'intrados.
    lowerAirflow: tableu des tableaus des valeurs des ordonnées de chaque ligne en dessous de l'aile.
    upperAirflow: tableu des tableaus des valeurs des ordonnées de chaque ligne au dessus de l'aile.
    """

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
        
    return(ex,ix,lowerAirflow,upperAirflow)
 

def plot_curves(file,nbUpperCurves,nbLowerCurves):
    plt.axis('equal' )
    (ex,ix,lowerAirflow,upperAirflow) = airflowModelling(file,nbUpperCurves,nbLowerCurves)
    for i in range(nbUpperCurves):
        plt.plot(ex,upperAirflow[i], '-k')
    for i in range(nbLowerCurves):
        plt.plot(ix,lowerAirflow[i], '-k')
    plt.show()
    



if __name__ == "__main__":
    if len(sys.argv) == 1:
        print( "Please provide a input file")
        exit(1)
    plot_curves( sys.argv[1], 20, 7 )
