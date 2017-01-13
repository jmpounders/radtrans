#include "mesh.h"
#include "log.h"

Mesh::Mesh(long nNodesReserve, long nElementsReserve)
  : MeshInterface()
{
  _nodeList.reserve(nNodesReserve);
  _elementList.reserve(nElementsReserve);
  _dimension = 1;
}

void
Mesh::make1DUniform(double width, unsigned long numCells)
{
  double meshSpacing = width/numCells;

  for (long i=0; i<numCells+1; i++) {
    _nodeList.push_back( Node( i*meshSpacing ) );
    _numNodes++;
  }

  for (long i=0; i<numCells; i++) {
    _elementList.push_back( LineElement( &_nodeList[i], &_nodeList[i+1] ) );
    _numElements++;
  }

}

void
Mesh::writeMesh(Output* outputFile)
{
  // Write 1D mesh
  double* outputBuffer = new double [ _numElements ];

  // Write position vector
  for (long i=0; i<_numElements; i++) {
    outputBuffer[i] = (_elementList[i].getNode(0)->get_x() + _elementList[i].getNode(1)->get_x())/2.0;
  }
  std::string varName = "/Mesh/x";
  outputFile->writeData(varName, outputBuffer, _numElements);
  
  delete [] outputBuffer;
}

