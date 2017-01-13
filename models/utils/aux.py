import numpy as np

def getAngleFromVector(dx, dy) :
    if abs(dx) < 1.0e-8 :
        if dy > 0.0 :
            phi = np.pi/2
        else :
            phi = 3.0*np.pi/2.0
    else :
        phi = (2.0*np.pi + np.arctan2(dy,dx)) % (2.0*np.pi)
    return phi

def isColinear(p0,p1,q) :
    if abs(p0[0]*(p1[1]-q[1]) + p1[0]*(q[1]-p0[1]) + q[0]*(p0[1]-p1[1])) < 1.0e-10 :
        return True
