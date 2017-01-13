import numpy as np

def recursiveIntegral(v0,v1,v2,f,tol) :
    int0 = triMidPntInt(v0,v1,v2,f)

    w0 = ( (v0[0]+v1[0])/2.0, (v0[1]+v1[1])/2.0 )
    w1 = ( (v1[0]+v2[0])/2.0, (v1[1]+v2[1])/2.0 )
    w2 = ( (v2[0]+v0[0])/2.0, (v2[1]+v0[1])/2.0 )

    int1 = triMidPntInt(v0,w0,w2,f)
    int2 = triMidPntInt(w0,v1,w1,f)
    int3 = triMidPntInt(w0,w1,w2,f)
    int4 = triMidPntInt(w2,w1,v2,f)

#    int1 = triVertAvg(v0,w0,w2,f)
#    int2 = triVertAvg(w0,v1,w1,f)
#    int3 = triVertAvg(w0,w1,w2,f)
#    int4 = triVertAvg(w2,w1,v2,f)

    intRef = int1+int2+int3+int4
    if intRef < tol :
        return (intRef, 1)
    
    if abs((int0 - intRef)/intRef) > tol :
        (int1,d1) = recursiveIntegral(v0,w0,w2,f,tol)
        (int2,d2) = recursiveIntegral(w0,v1,w1,f,tol)
        (int3,d3) = recursiveIntegral(w0,w1,w2,f,tol)
        (int4,d4) = recursiveIntegral(w2,w1,v2,f,tol)
        intRef = int1+int2+int3+int4
        d = max(d1,d2,d3,d4)+1
    else :
        d = 1

    return (intRef, d)


def triMidPntInt(v0,v1,v2,f) :
    A = triArea(v0,v1,v2)
    g = ( (v0[0]+v1[0]+v2[0])/3, (v0[1]+v1[1]+v2[1])/3 )
    return A*f(g[0],g[1])

def triVertAvg(v0,v1,v2,f) :
    A = triArea(v0,v1,v2)
    return A*(f(v0[0],v0[1])+f(v1[0],v1[1])+f(v2[0],v2[1]))/3.0

def triArea(v0,v1,v2) :
    return abs(v0[0]*(v1[1]-v2[1]) - v0[1]*(v1[0]-v2[0]) + v1[0]*v2[1] - v2[0]*v1[1])/2.0



def test() :
    v0 = [0.0, 0.0]
    v1 = [1.0, 0.0]
    v2 = [0.2,-1.0]
    def f(x,y) : return np.sin(x*np.pi/10)
    
    (intRef,d) = recursiveIntegral(v0,v1,v2,f,1.0e-6)
    A = triArea(v0,v1,v2)
    print intRef/A,d
    
