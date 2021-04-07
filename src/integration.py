import numpy as np

#integrate f between start and end by cuting the interval in n pieces
#and calling the aprox function on each piece 
#biger is n, more accurate the result will be
#Complexicity: O(n)
def integration(f,start,end,aprox,n):
    s = 0
    for k in range(n):
        s+= aprox(f(start+k*(end-start)/n), f(start+(k+1)*(end-start)/n), (end-start)/n)
    return s

def leftRectangle(y1,y2,step):
    return step*y1

def rightRectangle(y1,y2,step):
    return step*y2

def trapeze(y1,y2,step):
    return step*(y1+y2)/2
