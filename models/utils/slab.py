import numpy as np
from scipy.spatial import Delaunay
import h5py
import matplotlib.pyplot as plt

from geometry import simpleMesh

class slabGeometry :
    def __init__(self, slabWidth) :
        if len(slabWidth) == 1 :
            self.widthX = float(slabWidth[0])
            self.widthY = float(slabWidth[0])
        else :
            self.widthX = float(slabWidth[0])
            self.widthY = float(slabWidth[1])

        # Boundaries are given CCW starting north: N, W, S, E
        self.boundaries = ((self.widthX, self.widthY),
                           (        0.0, self.widthY),
                           (        0.0,       0.0),
                           (self.widthX,       0.0),
                           (self.widthX, self.widthY))
        self.boundaryPoints = [2, 2, 2, 2];
        self.mesh = None
        self.mat = None

def getDelaunayTriangulation(dims, disc) :
    if len(dims)==1 :
        widthX = float(dims[0])
        widthY = float(dims[0])
    else :
        widthX = float(dims[0])
        widthY = float(dims[1])

    if len(disc)==1 :
        nx = int(disc[0])
        ny = int(disc[0])
    else :
        nx = int(disc[0])
        ny = int(disc[1])
        
    dx = float(widthX)/(nx-1)
    dy = float(widthY)/(ny-1)

    nVerts = nx*ny
    verts = np.empty((nVerts,2))

    for j in range(0, ny) :
        for i in range(0, nx) :
            verts[ny*j + i, 0] = i*dx 
            verts[ny*j + i, 1] = j*dy

    #delMesh = Delaunay(verts)
    simplices = np.zeros((2*(nx-1)*(ny-1),3), dtype=np.int32)
    neighbors = np.zeros((2*(nx-1)*(ny-1),3), dtype=np.int32)
    for j in range(0, ny-1) :
        for i in range(0, nx-1) :
            simplices[j*2*(nx-1)+2*i,0] = j*nx+i
            simplices[j*2*(nx-1)+2*i,1] = j*nx+i+1
            simplices[j*2*(nx-1)+2*i,2] = (j+1)*nx+i
            neighbors[j*2*(nx-1)+2*i,0] = j*2*(nx-1)+2*i+1
            if i > 0 :
                neighbors[j*2*(nx-1)+2*i,1] = j*2*(nx-1)+2*i-1
            else :
                neighbors[j*2*(nx-1)+2*i,1] = -1
            if j > 0 :
                neighbors[j*2*(nx-1)+2*i,2] = (j-1)*2*(nx-1)+2*i+1
            else :
                neighbors[j*2*(nx-1)+2*i,2] = -1
                
            simplices[j*2*(nx-1)+2*i+1,0] = j*nx+i+1
            simplices[j*2*(nx-1)+2*i+1,1] = (j+1)*nx+i+1
            simplices[j*2*(nx-1)+2*i+1,2] = (j+1)*nx+i
            if j < ny-2 :
                neighbors[j*2*(nx-1)+2*i+1,0] = (j+1)*2*(nx-1)+2*i
            else :
                neighbors[j*2*(nx-1)+2*i+1,0] = -1
            neighbors[j*2*(nx-1)+2*i+1,1] = j*2*(nx-1)+2*i
            if i < nx-2 :
                neighbors[j*2*(nx-1)+2*i+1,2] = j*2*(nx-1)+2*i+2
            else :
                neighbors[j*2*(nx-1)+2*i+1,2] = -1
                
            materials = np.ones(simplices.shape[0], dtype=np.int32)
    return simpleMesh(verts, simplices, materials, neighbors, None)




def getDelaunayTriangulation2(dims, disc, block=1, origin=(0.0,0.0), pointsUp=True) :
    if len(dims)==1 :
        widthX = float(dims[0])
        widthY = float(dims[0])
    else :
        widthX = float(dims[0])
        widthY = float(dims[1])

    if len(disc)==1 :
        nx = int(disc[0])
        ny = int(disc[0])
    else :
        nx = int(disc[0])
        ny = int(disc[1])
        
    dx = float(widthX)/(nx-1)
    dy = float(widthY)/(ny-1)

    if ny%2 == 0 :
        nVerts = ny/2*(2*nx+1)
    else :
        nVerts = (ny-1)/2*(2*nx+1) + nx
        
    verts = np.empty((nVerts,2))

    vc = 0
    if pointsUp :
        for j in range(0, ny) :
            if j%2 == 0 :
                for i in range(0, nx) :
                    verts[vc,0] = i*dx
                    verts[vc,1] = j*dy
                    vc = vc + 1
            else :
                verts[vc,0] = 0.0
                verts[vc,1] = j*dy
                vc = vc + 1
                for i in range(0, nx-1) :                            
                    verts[vc,0] = (i+0.5)*dx
                    verts[vc,1] = j*dy
                    vc = vc + 1
                verts[vc,0] = widthX
                verts[vc,1] = j*dy
                vc = vc + 1
    else :
        for j in range(0, ny) :
            if j%2 == 1 :
                for i in range(0, nx) :
                    verts[vc,0] = i*dx
                    verts[vc,1] = j*dy
                    vc = vc + 1
            else :
                verts[vc,0] = 0.0
                verts[vc,1] = j*dy
                vc = vc + 1
                for i in range(0, nx-1) :                            
                    verts[vc,0] = (i+0.5)*dx
                    verts[vc,1] = j*dy
                    vc = vc + 1
                verts[vc,0] = widthX
                verts[vc,1] = j*dy
                vc = vc + 1

    delMesh = Delaunay(verts)
    materials = np.ones(delMesh.simplices.shape[0], dtype=np.int32)*block
    delMesh.points[:,0] += origin[0]
    delMesh.points[:,1] += origin[1]
    convex_hull = [[delMesh.simplices[s,i] for i in range(3) if delMesh.neighbors[s,i] > -1]
                   for s in range(len(delMesh.simplices)) if -1 in delMesh.neighbors[s,:]]
    convex_hull = np.array(convex_hull)
    bdries = [s for s in range(len(delMesh.simplices)) if -1 in delMesh.neighbors[s,:]]
    return simpleMesh(delMesh.points,
                      delMesh.simplices,
                      materials,
                      delMesh.neighbors,
                      convex_hull)
