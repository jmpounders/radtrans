
from collections import deque

class Vertex :
    def __init__(self, id) :
        self.id = id
        self.downwindVertices = []
        self.upwindVertices = []

    def __str__(self) :
        return ("Vertex " + str(self.id) + ":\n" +
                "   dw connected to: " + str([x.id for x in self.downwindVertices]) + "\n" +
                "   uw connected to: " + str([x.id for x in self.upwindVertices]))
#        return ("Vertex " + str(self.id) + ":\n" +
#                "   dw connected to: " + str([x.id for x in self.downwindVertices]))
        
    def addDownwind(self, dwVert) :
        self.downwindVertices.append(dwVert)

    def removeDownwind(self, dwVert) :
        self.downwindVertices.remove(dwVert)

    def addUpwind(self, upVert) :
        self.upwindVertices.append(upVert)

    def removeUpwind(self, upVert) :
        self.upwindVertices.remove(upVert)

    def getDownwindVertices(self) :
        return self.downwindVertices.keys()

    def getID(self) :
        return self.id

class Graph :
    def __init__(self) :
        self.vertexList = {}
        self.numVertices = 0
        self.roots = []

    def __contains__(self, key) :
        return key in self.vertexList

    def __iter__(self) :
        return iter(self.vertexList.values())

    def addVertex(self, key) :
        self.numVertices += 1
        newVertex = Vertex(key)
        self.vertexList[key] = newVertex

    def getVertex(self, key) :
        if key in self.vertexList :
            return self.vertexList[key]
        else :
            return None

    def addEdge(self, parent, child) :
        if parent not in self.vertexList :
            self.addVertex(parent)
        if child not in self.vertexList :
            self.addVertex(child)
        self.vertexList[parent].addDownwind(self.vertexList[child])
        self.vertexList[child].addUpwind(self.vertexList[parent])

    def removeEdge(self, parent, child) :
        self.vertexList[parent].removeDownwind(self.vertexList[child])
        self.vertexList[child].removeUpwind(self.vertexList[parent])

    def getVertices(self) :
        return self.vertexList.keys()

    def getRoots(self) :
        self.roots = []
        for v in self.vertexList.values() :
            if len(v.upwindVertices) == 0 : self.roots.append(v.id)
        return self.roots

    def getBreadthFirstOrder(self) :
        self.getRoots()
        discovered = [False] * self.numVertices

        order = self.roots[:]
        for i in order : discovered[i] = True
        
        currentIndex = 0
        while currentIndex < self.numVertices :
            v = self.vertexList[order[currentIndex]]
            discovered[v.id] = True
            for w in v.downwindVertices :
                if not discovered[w.id] :
                    order.append(w.id)
                    discovered[w.id] = True
            currentIndex += 1
        return order

    def orderVertices(self) :
        L = []
        S = self.getRoots()
        while len(S) > 0 :
            n = S.pop()
            L.append(n)
            childList = self.vertexList[n].downwindVertices[:]
            for m in childList :
                self.removeEdge(n, m.id)
                if len(m.upwindVertices) == 0 :
                    S.append(m.id)
        return L
                

    
def test() :
    g = Graph()

    print g.vertexList
    g.addEdge(0,1)
    g.addEdge(0,2)
    g.addEdge(2,3)
    g.addEdge(4,0)
    g.addEdge(4,5)
    g.addEdge(5,6)
    g.addEdge(6,7)
    g.addEdge(6,8)
    g.addEdge(9,5)
    g.addEdge(9,10)

    for v in g :
        print v

    print "The roots are ", g.getRoots()
    print "BFO = ", g.getBreadthFirstOrder()
    print "Order = ", g.orderVertices()
