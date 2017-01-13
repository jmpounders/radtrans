#ifndef MESH_H
#define MESH_H

#include<vector>

#include "meshinterface.h"

/// One-dimensional mesh
/**
 *  This "native" mesh could be generalized to higher dimension,
 *  but it has only been used for one-dimensional problems.  The
 *  structure allows for an arbitrary unstructured mesh.
 */
class Mesh : public MeshInterface
{
 public:
  Mesh(long nNodesReserve, long nElementsReserve);

  void make1DUniform(double width, unsigned long numCells);

  // Un-implemented base-class features
  void loadMesh(std::string inputFileName) {};
  void tagMesh(std::string tagName, double* tagDataBuffer, long tagSize) {};
  void readMeshSweepOrder(const std::string fileName) {};
  void createDefaultMeshSweepOrder() {};

  // Concrete implementation of abstract base
  void writeMesh(Output* outputFile);
  Element* getElement(long elemID) { return &_elementList[elemID]; };
  double getElementVolume( long elemID ) { return _elementList[elemID].getVolume(); };
  Material* getElementMat ( long elemID ) { return _elementList[elemID].getMat(); };

 private:
  std::vector<Node> _nodeList;
  std::vector<LineElement> _elementList;
};

#endif
