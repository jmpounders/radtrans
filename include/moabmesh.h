#ifndef MOABMESH
#define MOABMESH

#include "meshinterface.h"
#include "moab/Core.hpp"

/// Unstructured mesh using the MOAB library
class MoabMesh : public MeshInterface
{
  friend class Sweeper;
 public:
  MoabMesh();
  ~MoabMesh();

  long getEdgeID(long node1, long node2);

  void setElementMat(long elementID, Material* mat) { _materialList[elementID] = mat; };
  Material* getElementMat ( long elementID ) { return _materialList[elementID]; };

  long queryPointLocation(  double* point );
  long getAdjacentTriangle( long triangle, long vertex1, long vertex2 );

  void getCurrentElementFromID(long elementID, UltraLightElement& elementUL);


  // Un-implemented base-class features
  Element* getElement(long elemID) { return NULL; };
  
  // Concrete implementation of abstract base class
  void loadMesh(std::string inputFileName);
  void writeMesh(Output* outputFile);
  void tagMesh(std::string tagName, double* tagDataBuffer, long tagSize);
  void readMeshSweepOrder(const std::string fileName);
  void createDefaultMeshSweepOrder() {} ;
  double getElementVolume( long elemID );

  // MOAB-specific implementations
  void logMemoryUse();

  std::vector<long> boundaryElements;

 private:
  moab::Interface* _mb;
  
  std::vector<moab::EntityHandle> _orderedMeshSet;
  std::vector<moab::EntityHandle> _tris;
  moab::Range _verts;

  void _identifyEdges();
  void _identifyBoundaryElements();
  std::map<std::pair<long, long>, long> _edgeID;

  std::vector<Material*> _materialList;

  moab::EntityHandle* _connectivity;
  double* _x;
  double* _y;
  double* _z;
  
  moab::ErrorCode _rval;
};

#endif
