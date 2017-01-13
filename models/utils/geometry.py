import numpy as np
from scipy.spatial import Delaunay
import h5py
import matplotlib.pyplot as plt

import dag
import aux

import pdb

class simpleMesh :
    def __init__(self, vertices, connectivity, materials, simplexNeighbors, hull) :
        self.points = vertices
        self.simplices = connectivity
        self.materials = materials
        self.neighbors = simplexNeighbors
        self.convex_hull = hull

    def translate(self,dx,dy) :
        for i in range(len(self.points)) :
            self.points[i,0] += dx
            self.points[i,1] += dy

    def rotate(self,theta) :
        for i in range(len(self.points)) :
            x = self.points[i,0]
            y = self.points[i,1]
            self.points[i,0] = x*np.cos(theta) - y*np.sin(theta)
            self.points[i,1] = x*np.sin(theta) + y*np.cos(theta)

    def flip(self,m,b) :
        for i in range(len(self.points)) :
            x = self.points[i,0]
            y = self.points[i,1]
            d = (x + (y - b)*m)/(1.0 + m**2)
            self.points[i,0] = 2.0*d - x
            self.points[i,1] = 2.0*d*m - y + 2.0*b
        
    def plotMesh(self) :
        fig = plt.figure()
        plt.triplot(self.points[:,0],
                    self.points[:,1],
                    self.simplices.copy())
        plt.axis('equal')
        plt.show()

    def plotGeom(self) :
        fit = plt.figure()
        plt.tripcolor(self.points[:,0],
                      self.points[:,1],
                      self.simplices,
                      facecolors=self.materials,
                      edgecolors='k')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.axis('equal')
        plt.show()

    def write_h5m(self, fileName) :
   
        simplexShape  = self.simplices.shape
        nTris = simplexShape[0]
        nVerts = self.points.shape[0]
    
        outputFile = h5py.File(fileName,'w')
        tsttGroup = outputFile.create_group("/tstt")
        tsttGroup.attrs["max_id"] = nVerts + nTris
    
        elemEnum = {"Edge":1,
                    "Tri":2,
                    "Quad":3,
                    "Polygon":4,
                    "Tet":5,
                    "Pyramid":6,
                    "Prism":7,
                    "Knife":8,
                    "Hex":9,
                    "Polyhedron":10}


        nodeList = outputFile.create_dataset("/tstt/nodes/coordinates",
                                             (nVerts,2), dtype='f8',
                                             data=self.points)
        nodeList.attrs["start_id"] = 1
    
        elementGroup = outputFile.create_group("/tstt/elements/Tri3")
        elementGroup.attrs.create("element_type",
                                  elemEnum["Tri"],
                                  dtype=h5py.special_dtype(enum=("i1",elemEnum)))
        elementList = outputFile.create_dataset("tstt/elements/Tri3/connectivity",
                                                simplexShape, dtype='i',
                                                data=self.simplices+1)
        elementList.attrs["start_id"] = nVerts + 1
    
        tagList = outputFile.create_group("/tstt/tags")
        setList = outputFile.create_group("/tstt/sets")
        outputFile.close()

def mergeMeshList(meshList) :
    bigMesh = meshList.pop(0)
    for mesh in meshList :
        bigMesh = mergeMesh(bigMesh, mesh)
    return bigMesh

def mergeMesh(mesh1, mesh2) :
    nVerts1 = mesh1.points.shape[0]
    nTris1 = mesh1.simplices.shape[0]
    nVerts2 = mesh2.points.shape[0]
    nTris2 = mesh2.simplices.shape[0]

    # Identify duplicate points between slab1 and slab2
    eps = 1.0e-14
    dupPoints = []
    for p2 in set(mesh2.convex_hull.flatten()) :
        xy2 = mesh2.points[p2]
        matchX = np.where(abs(mesh1.points[:,0]-xy2[0]) < eps)
        matchY = np.where(abs(mesh1.points[:,1]-xy2[1]) < eps)
        p1 = list(set(matchX[0]) & set(matchY[0]))
        if len(p1)==0 : continue
        assert len(p1) == 1,p1
        dupPoints.append([p1[0],p2])

    # Create mask for excluding duplicate points
    dupPoints = np.array(dupPoints)
    mask = np.ones(nVerts2, dtype=bool)
    mask[dupPoints[:,1]] = False
    p2Lookup = np.zeros(nVerts2, dtype=np.int32)
    p2Lookup[mask] = np.array(range(nVerts2-len(dupPoints))) + nVerts1
    for p in dupPoints :
        p2Lookup[p[1]] = p[0]

    # Form the merged arrays and offset where necessary
    points = np.concatenate((mesh1.points, mesh2.points[mask]))
    simplices = np.concatenate((mesh1.simplices, mesh2.simplices))
    materials = np.concatenate((mesh1.materials, mesh2.materials))
    neighbors = np.concatenate((mesh1.neighbors, mesh2.neighbors))
    for s in range(nTris2) :
        for v in range(3) :
            simplices[nTris1+s,v] = p2Lookup[mesh2.simplices[s,v]]
            if neighbors[nTris1+s,v]>=0 : neighbors[nTris1+s,v] += nTris1

    # Identify the new convex_hull and the new interior surface
    interior = []
    convex_hull = []
    for face in mesh1.convex_hull :
        if face[0] in dupPoints[:,0] and face[1] in dupPoints[:,0] :
            interior.append(face)
        else :
            convex_hull.append(face)
    for face in mesh2.convex_hull :
        if face[0] not in dupPoints[:,1] or face[1] not in dupPoints[:,1] :
            convex_hull.append([p2Lookup[face[0]], p2Lookup[face[1]]])
    convex_hull = np.array(convex_hull)
    

    # Fix-up new neighboring simplices
    bdrySimps = [s for s in range(nTris1+nTris2) if -1 in neighbors[s]]
    nghborMatch = {tuple(face):[] for face in interior}
    for face in nghborMatch :
        for s in bdrySimps :
            if set(face) < set(simplices[s]) : nghborMatch[face].append(s)
    
    for face,pair in nghborMatch.iteritems() :
        for v in range(3) :
            if set([simplices[pair[0],v]])|set(face) == set(simplices[pair[0],:]) : s0n = v
            if set([simplices[pair[1],v]])|set(face) == set(simplices[pair[1],:]) : s1n = v
        neighbors[pair[0],s0n] = pair[1]
        neighbors[pair[1],s1n] = pair[0]
        for n in neighbors[pair[0],:] :
            if n < 0 : continue
        for n in neighbors[pair[1],:] :
            if n < 0 : continue

    return simpleMesh(points, simplices, materials, neighbors, convex_hull)
    
def getSweepOrder(mesh, omega_x, omega_y, omega_z, numProcs = 1) :
    theta = np.arccos(omega_x/np.sqrt(1.0-omega_z**2))
    sweepOrder = np.zeros((mesh.simplices.shape[0],len(theta)),dtype=np.int32)
    for (i,thetai) in enumerate(theta) :
        if omega_y[i] < 0.0 : thetai = 2.0*np.pi - thetai
        order = getSweepOrderTheta(mesh, thetai)
        if not len(order)==mesh.simplices.shape[0] :
            print thetai, len(order), mesh.simplices.shape[0]
        sweepOrder[:,i] = np.array(order)
    return sweepOrder

        
def getSweepOrderTheta(mesh, thetaInc) :
    g = dag.Graph()

    for currentCell in range(mesh.simplices.shape[0]) :
        v = mesh.simplices[currentCell]
        n = mesh.neighbors[currentCell]

        phi1 = aux.getAngleFromVector(mesh.points[v[1],0] - mesh.points[v[0],0],
                                      mesh.points[v[1],1] - mesh.points[v[0],1])
        phi2 = aux.getAngleFromVector(mesh.points[v[2],0] - mesh.points[v[0],0],
                                      mesh.points[v[2],1] - mesh.points[v[0],1])
        if phi1 < 1.0e-10 and phi2 > np.pi : phi1 = 2*np.pi
        if phi2 < 1.0e-10 and phi1 > np.pi : phi2 = 2*np.pi

        theta = np.zeros(3)
        d = np.zeros(3)
        for i in range(3) :
            dx0 = mesh.points[v[(i+1)%3],0] - mesh.points[v[i%3],0]
            dy0 = mesh.points[v[(i+1)%3],1] - mesh.points[v[i%3],1]
            dx1 = mesh.points[v[(i+2)%3],0] - mesh.points[v[i%3],0]
            dy1 = mesh.points[v[(i+2)%3],1] - mesh.points[v[i%3],1]
            theta[i] = np.arccos((dx1*dx0 + dy1*dy0)/(np.sqrt( dx0**2 + dy0**2 )*np.sqrt( dx1**2 + dy1**2 )))
            d[i] = np.sqrt( dx0**2 + dy0**2 )

        # make sure vertices are oriented counter-clockwise
        det = ((mesh.points[v[1],0]-mesh.points[v[0],0])*(mesh.points[v[2],1]-mesh.points[v[0],1]) -
               (mesh.points[v[2],0]-mesh.points[v[0],0])*(mesh.points[v[1],1]-mesh.points[v[0],1]))
        if det < 0 :
            phitemp = phi1
            phi1 = phi2
            phi2 = phitemp
                
            thetatemp = theta[1]
            theta[1] = theta[2]
            theta[2] = thetatemp

            ntemp = n[1]
            n[1] = n[2]
            n[2] = ntemp

            vtemp = v[1]
            v[1] = v[2]
            v[2] = vtemp

        thetap = thetaInc - phi1;
        if thetap < 0.0: thetap += 2.0*np.pi
        if abs(thetap-2.0*np.pi) < 1.0e-9 : thetap = 0.0

        if thetap % np.pi <= theta[0] - 1.0e-9:
            if thetap < np.pi :
                # v0 to e12
                nextCell = (n[0],)
            else :
                # e12 to v0
                nextCell = (n[1],n[2])
        elif thetap % np.pi <= np.pi - theta[1] - 1.0e-9 :
            if thetap < np.pi :
                # e01 to v2
                nextCell = (n[0],n[1])
            else :
                # v2 to e01
                nextCell = (n[2],)
        else :
            if thetap < np.pi :
                # v1 to e20
                nextCell = (n[1],)
            else :
                # e20 to v1
                nextCell = (n[0],n[2])

        for e in nextCell :
            if e >=0 : g.addEdge(currentCell, e)

    #for v in g :
    #    print v
    #print "The roots are ", g.getRoots()
    #print "Order = ", g.orderVertices()
    #order = g.getBreadthFirstOrder()
    order = g.orderVertices()
    return order
