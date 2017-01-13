#ifndef MESHINTERFACE_H
#define MESHINTERFACE_H

#include <string>

#include "element.h"
#include "material.h"
#include "output.h"

/// Virual mesh interface
class MeshInterface
{
 public:
  MeshInterface() : _numNodes(0), _numElements(0) {};
  virtual ~MeshInterface() {};

  virtual void loadMesh(std::string inputFileName) = 0;
  virtual void writeMesh(Output* outputFile) = 0;

  virtual void tagMesh(std::string tagName, double* tagDataBuffer, long tagSize) = 0;

  virtual Element* getElement(long elemID) = 0;
  virtual double getElementVolume( long elemID ) = 0;
  virtual Material* getElementMat ( long elemID ) = 0;

  long numNodes() { return _numNodes; };
  long numEdges() { return _numEdges; };
  long numElements() { return _numElements; };
  int dim() { return _dimension; };

  virtual void readMeshSweepOrder(const std::string fileName) = 0;
  virtual void createDefaultMeshSweepOrder() = 0;

 protected:
  int _dimension;
  long _numNodes;
  long _numEdges;
  long _numElements;
};

#endif
