import numpy as np
import numpy as np
from scipy.spatial import Delaunay
import h5py
import matplotlib.pyplot as plt

from geometry import simpleMesh

class triangleGeometry :
    def __init__(self, angles, side0) :
        if len(angles) >= 2 and type(side0) == float :
            self.angles = (angles[0], angles[1], np.pi-angles[0]-angles[1])
            self.sides = (side0,
                          side0*np.sin(self.angles[0])/np.sin(self.angles[2]),
                          side0*np.sin(self.angles[1])/np.sin(self.angles[2]))
        else :
            raise RuntimeError("Bad triangle specification")
                
        self.mesh = None
        self.mat = None


def getTriangleTriangulation(tri,n,block=1) :
    nVerts = n*(n+1)/2
    verts = np.zeros((nVerts,2))
    hx = float(tri.sides[0])/(n-1)
    hy = float(tri.sides[2]*np.sin(tri.angles[0]))/(n-1)
    xStep = float(tri.sides[2]*np.cos(tri.angles[0]))/(n-1)
    
    vc = 0
    for j in range(n) :
        for i in range(n-j) :
            verts[vc,0] = j*xStep + i*hx
            verts[vc,1] = j*hy
            vc += 1
            
    delMesh = Delaunay(verts)
    materials = np.ones(delMesh.simplices.shape[0], dtype=np.int32)*block
    #convex_hull = [[delMesh.simplices[s,i] for i in range(3) if delMesh.neighbors[s,i] > -1]
    #               for s in range(len(delMesh.simplices)) if -1 in delMesh.neighbors[s,:]]
    #convex_hull = np.array(convex_hull)
    convex_hull = delMesh.convex_hull
    bdries = [s for s in range(len(delMesh.simplices)) if -1 in delMesh.neighbors[s,:]]
    return simpleMesh(delMesh.points,
                      delMesh.simplices,
                      materials,
                      delMesh.neighbors,
                      convex_hull)
            

