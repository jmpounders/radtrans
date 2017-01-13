#!/usr/bin/python

def getSolution(outputFileName) :
    import h5py
    import numpy

    outputFile = h5py.File(outputFileName, "r")

    x = outputFile["/Mesh/x"][:]
    flux = []

    group = 0
    while True :
        group = group + 1
        varName = "scalarFlux_" + str(group)
        if (varName not in outputFile) : break
        groupFlux = outputFile[varName][:]
        flux.append(groupFlux)

    if ("k" in outputFile) :
        k = outputFile["k"][:]
        k = k[0]
    else :
        k = float("nan")
        
    outputFile.close()

    fluxArray = numpy.empty([x.shape[0], len(flux)])
    for g in range(len(flux)) :
        fluxArray[:,g] = flux[g]
        
    return (x, fluxArray, k)


def getTransientSolution(outputFileName) :
    import h5py
    import numpy

    outputFile = h5py.File(outputFileName, "r")

    time = outputFile["time"][:]
    productionRate = outputFile["totalNeutronProduction"][:]

    outputFile.close()

    return (time, productionRate)
