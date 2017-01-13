import numpy as np

class vtkFile :
    def __init__(self, fileName) :
        f = open(fileName, 'r')
        self.version = f.readline().rstrip('\r\n')
        self.header = f.readline().rstrip('\r\n')
        frmat = f.readline().rstrip('\r\n')
        if frmat.upper() != "ASCII" :
            print "Cant read binary vtk files"
            return
        self.grid = f.readline().rstrip('\r\n').split()[1]
        if self.grid.upper() != "UNSTRUCTURED_GRID" :
            print "Cant handle anything except unstructured triangular grids"
            return
        self.variables = {}

        record = f.readline()
        while record :
            record = record.rstrip('\r\n').split()
            if record[0].upper() == "POINTS" :
                self.numPoints = int(record[1])
                pointType = record[2]
                self.points = np.zeros((self.numPoints,3))
                for i in range(self.numPoints) :
                    self.points[i,:] = f.readline().rstrip('\r\n').split()
                
            elif record[0].upper() == "CELLS" :
                self.numCells = int(record[1])
                connectivityLength = int(record[2])
                if float(connectivityLength)/self.numCells != 4 :
                    print "This does not appear to be a triangulation."
                    return
                self.connectivity = np.zeros((self.numCells,3),np.int64)
                for i in range(self.numCells) :
                    self.connectivity[i,:] = f.readline().rstrip('/r/n').split()[1:4]
                
            elif record[0].upper() == "CELL_TYPES" :
                numCellTypes = int(record[1])
                for i in range(numCellTypes) :
                    f.readline()
                
            elif record[0].upper() == "CELL_DATA" :
                numCellData = int(record[1])
                dataAttribute = f.readline().rstrip('\r\n').split()
                if dataAttribute[0].upper() == "SCALARS" :
                    varName = dataAttribute[1]
                    varType = dataAttribute[2]
                    lookupTable = f.readline().rstrip('\r\n').split()[0]
                    if lookupTable.upper() != "LOOKUP_TABLE" :
                        print "Expected to find LOOKUP_TABLE; not found."
                        return
                    print "Reading ", varName
                    var = np.zeros(numCellData)
                    for i in range(numCellData) :
                        var[i] = float(f.readline().rstrip('/r/n'))
                    self.variables[varName] = var
                
            record = f.readline()
        
        return

    def getSolutionAtPoint(self,varName,p) :
        print p
        for cell in range(self.numCells) :
            x = self.points[self.connectivity[cell,:],0]
            y = self.points[self.connectivity[cell,:],1]
            x[1] = x[1]-x[0]
            y[1] = y[1]-y[0]
            x[2] = x[2]-x[0]
            y[2] = y[2]-y[0]
            d = x[1]*y[2]-y[1]*x[2]
            a =  (p[0]*y[2] - p[1]*x[2] - x[0]*y[2] + y[0]*x[2])/d
            b = -(p[0]*y[1] - p[1]*x[1] - x[0]*y[1] + y[0]*x[1])/d

            if a>0.0 and b>0.0 and (a+b)<1.0 :
                return self.variables[varName][cell]
