# Handmade tools
import utils.slab as slab
import utils.geometry as geometry
import utils.quadrature as quadrature
from utils.inputWriter import *
import utils.vtkReader as vtkReader

# Standard tools
import numpy as np
import h5py
import os
import pickle

import matplotlib.pyplot as plt

def makeProblem(slabWidth, sigma_t, nx, ny, extSource, solverMethod) :

    meshFileName = "hs_mesh_" + str(nx) + ".h5m"

    # Get quadrature
    (omega_x,omega_y,omega_z,wgt) = quadrature.getQuadrature3DHalf(2)

    # Make geometry
    print "  Making geometry"
    geom = slab.slabGeometry([slabWidth])
    geom.mesh = slab.getDelaunayTriangulation2([slabWidth], (nx,ny))
    geom.mesh.write_h5m(meshFileName)

    # Get sweep order
    print "  Getting sweep order"
    sweepFileName = "hs_sweeporder"+str(nx)
    if os.path.isfile(sweepFileName+".h5") :
        print "    Sweep file already exists... not recreating."
    else :
        sweepOrder = geometry.getSweepOrder(geom.mesh, omega_x, omega_y, omega_z)
        outputFile = h5py.File(sweepFileName+".h5",'w')
        outputFile.create_dataset("sweeporder",
                                  sweepOrder.shape,
                                  dtype='i8',
                                  data=sweepOrder)
        outputFile.close()
        print sweepOrder.shape

    # Make material file
    print "  Writing materials"
    materialFile = h5py.File("material.h5",'w')
    materialFile.create_dataset("/materials",
                                geom.mesh.materials.shape,
                                dtype='i8',
                                data=geom.mesh.materials)
    materialFile.close()

    # Make source
    source = np.zeros((geom.mesh.simplices.shape[0],len(omega_x)))
    for n in range(len(omega_x)) :
        for (e,conn) in enumerate(geom.mesh.simplices) :
            source[e,n] =  extSource

    sourceFile = h5py.File("source.h5",'w')
    sourceFile.create_dataset("/source",
                              source.shape,
                              dtype='f8',
                              data=source)
    sourceFile.close()



    # Make input
    # Mesh Input###############################################
    mesh = InputSet("Mesh")
    mesh.addData("name", "mesh1")
    mesh.addData("type", "moab")
    mesh.addData("sweep", sweepFileName)
    mesh.addData("file", meshFileName)
    ###########################################################

    # Material Input ##########################################
    material = InputSet("Material")
    material.addData("name",            "mat1")
    material.addData("numGroups",       1)
    material.addData("block",           1)
    material.addData("sigmaT",          sigma_t)
    material.addData("sigmaS",          sigma_s)
    ###########################################################

    # Transport Problem Input #################################
    transport = InputSet("Transport")
    transport.addData("mesh",      "mesh1")
    transport.addData("quadOrder", len(wgt))
    transport.addData("omega_x",   omega_x.tolist())
    transport.addData("omega_y",   omega_y.tolist())
    transport.addData("omega_z",   omega_z.tolist())
    transport.addData("weights",   wgt.tolist())

    externalSource = InputSet("ExternalSource")
    externalSource.addData("type",      "uniform")
    #externalSource.addData("magnitude", extSource)
    externalSource.addData("type", "file")
    externalSource.addData("name", "source")
    transport.add(externalSource)
    ###########################################################

    # Solver Input ############################################
    solver = InputSet("Solver")
    solver.addData("type", solverMethod)
    ###########################################################

    # Problem Manager Input ###################################
    problemManager = InputSet("FixedSource")
    ###########################################################

    inputRoot = InputSet("Input")
    inputRoot.add(mesh)
    inputRoot.add(material)
    inputRoot.add(transport)
    inputRoot.add(solver)
    inputRoot.add(problemManager)

    inputFile = open("input", "w")
    inputRoot.write(inputFile)
    inputFile.close()

    # Write session file
    with open("session.pickle", 'wb') as session :
        pickle.dump(geom, session, pickle.HIGHEST_PROTOCOL)

        
def postProcess() :
    results = {}

    # Read session file
    with open("session.pickle", 'rb') as session :
        geom = pickle.load(session)

    # Read output
    output = vtkReader.vtkFile("output.vtk")

    # Average fluxes by region/block
    regFluxes = [0.0]
    regVolumes = [0.0]
    for (si,s) in enumerate(geom.mesh.simplices) :
        v0 = geom.mesh.points[s[0]]
        v1 = geom.mesh.points[s[1]]
        v2 = geom.mesh.points[s[2]]
        area = abs((v0[0]-v2[0])*(v1[1]-v0[1])-(v0[0]-v1[0])*(v2[1]-v0[1]))/2.0
        regVolumes[geom.mesh.materials[si]-1] += area
        regFluxes[geom.mesh.materials[si]-1] += output.variables['scalarFlux_1'][si]*area

    for ri in range(len(regFluxes)) :
        results["reg"+str(ri+1)+"Flux"] = regFluxes[ri]/regVolumes[ri]

    plt.figure()
    plt.tripcolor(geom.mesh.points[:,0],geom.mesh.points[:,1], geom.mesh.simplices, output.variables['scalarFlux_1'])
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Numerical Solution')
    plt.savefig('flux.png', bbox_inches='tight')
    plt.show()

    return results

# Problem definition
if __name__ == "__main__" :
    import sys

    if "make" in sys.argv :
        print "Making problem input"
        slabWidth = 1.0
        sigma_t = 0.1
        c = 0.0;  sigma_s = sigma_t*c;
        nPoints = (17,17)
        extSource = 1.0/(4.0*np.pi)
        solver = "localMOC"
        solver = "regMOC"
        makeProblem(slabWidth, sigma_t, nPoints[0], nPoints[1], extSource, solver)

    if "run" in sys.argv :
        print "Running problem"
        os.system("time ../build/Transport input | tee output.screen")

    if "process" in sys.argv :
        print "Post processing results"
        results = postProcess()
        print results

    if "parametric" in sys.argv :
        import shutil
        import json
        print "Running parametric study"
        slabWidth = 10.0
        extSource = 1.0/(4.0*np.pi)
        sigma_ts = (0.1,0.01)
        cs = (0.0, 0.5, 0.9)
        nPoints = ((5,5),(9,9),(17,17),(33,33),(65,65),(129,129),(1025,1025))
        solvers = ("localMOC","regMOC")
        outputFiles = ("output.screen",
                       "output.log",
                       "output.vtk",
                       "flux.png")
        results = []
        name = "homogSlabParam"
        for disc in nPoints :
            for solver in solvers :
                for c in cs :
                    for sigma_t in sigma_ts :
                        sigma_s = sigma_t*c;
                        # Setup this case
                        pid = name+"/_"+solver+"_"+str(disc[0])+"_"+str(c)+"_"+str(sigma_t)
                        os.makedirs(pid)
                        makeProblem(slabWidth, sigma_t, disc[0], disc[1], extSource, solver)

                        # Run the problems (this will only work *NIX-based systems)
                        print "  Case "+pid
                        os.system("time ../build/Transport input | tee output.screen")

                        # Get the results and copy output files
                        caseResults = {"solver":solver,
                                       "nx":disc[0],
                                       "ny":disc[1],
                                       "c":c,
                                       "sigma_t":sigma_t}
                        caseResults.update(postProcess())
                        results.append(caseResults)
                        for outFile in outputFiles :
                            shutil.copy(outFile, os.path.join(pid, outFile))

                        # Save a database of the results
                        # this JSOB DB can be read several ways
                        # In R,
                        # > library(jsonlite)
                        # > results <- fromJSON(readLines("results.json"))
                        f = open("results.json",'w')
                        json.dump(results, f, indent=4, sort_keys=True)
                        f.close()

        shutil.copy("results.json", os.path.join(name, "results.json"))            

