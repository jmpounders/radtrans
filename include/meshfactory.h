#ifndef MESHFACTORY_H
#define MESHFACTORY_H

#include "meshinterface.h"

enum mesh_t {NATIVE, MOAB};

/// Factory for creating mesh and keeping a list of available mesh
class MeshFactory
{
 public:
  void createMeshFromInput(InputParser& input, MaterialFactory& materialFactory);
  void deleteAllMesh();
  MeshInterface* createNewMesh(std::string& meshName, mesh_t meshType, int nNodesReserve=0, int nElementsReserve=0);
  MeshInterface* getMesh(std::string& materialName);
  void list();
  
 private:
  std::map<std::string, MeshInterface*> _meshMap;

};

#endif
