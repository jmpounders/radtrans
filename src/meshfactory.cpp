#include "meshfactory.h"
#include "mesh.h"
#include "moabmesh.h"
#include "log.h"

void
MeshFactory::createMeshFromInput(InputParser& input, MaterialFactory& materialFactory)
{
  std::cout << "Making mesh... ";

  std::vector<std::string> path(1,"Mesh");
  std::vector<double> v;

  // Get the name of the mesh
  std::string meshName = input.getString(path, "name");
  if (meshName == "empty") {
    LOG_ERR("Must give the mesh a name; use the name parameter.");
  }

  // Load or create a mesh based on the mesh type
  // The mesh type determines which derivative of the MeshInterface class gets instantiated
  mesh_t meshType;
  std::string meshTypeStr = input.getString(path, "type");
  if (meshTypeStr == "moab") {
    // Instantiate the class based on the MOAB mesh library
    LOG("Creating MOAB mesh.");
    meshType = MOAB;
    MeshInterface* moabMesh = createNewMesh(meshName, meshType);

    std::string meshFileName = input.getString(path, "file");
    moabMesh->loadMesh(meshFileName);

    std::string sweepFileName = input.getString(path, "sweep");
    if (sweepFileName == "empty")
      moabMesh->createDefaultMeshSweepOrder();
    else
      moabMesh->readMeshSweepOrder(sweepFileName);

    std::string fileName("material");
    std::string varName("/materials");
    unsigned int* matArray;
    HDF5Interface h5;
    h5.open(fileName,'R');
    HDFData* materialData = h5.readData(fileName, varName);
    matArray = (unsigned int*)materialData->data;
    h5.close(fileName);

    MoabMesh* meshp = static_cast<MoabMesh*>(moabMesh);
    for (int elementID=0; elementID<meshp->numElements(); elementID++)
      meshp->setElementMat( elementID, materialFactory.getMaterial(matArray[elementID]) );

    delete [] (unsigned int*)materialData->data;
    delete materialData;

  }
  else if (meshTypeStr == "native" || meshTypeStr == "empty") {
    // Instantiate the native mesh class (currently 1D only)
    LOG("Creating native mesh.");
    meshType = NATIVE;
    
    v = input.getVector(path, "slabWidth");
    if (v.size() == 0) {
      LOG_ERR("Must specify the slab width; use the slabWidth parameter.");
    }
    double slabWidth = v[0];

    v = input.getVector(path, "numCells");
    if (v.size() == 0) {
      LOG_ERR("Must specify number of cells in mesh; use the numCells parameter.");
    }
    int numCells = int(v[0]);

    MeshInterface* mesh1d = createNewMesh(meshName, meshType, numCells+1, numCells);
    static_cast<Mesh*>(mesh1d)->make1DUniform( slabWidth, numCells);

    // Assign materials to mesh
    std::string matName = "mat1";
    for (int i=0; i<numCells; i++) {
      mesh1d->getElement(i)->setMat( materialFactory.getMaterial(matName) );
    }
  }
  
  std::cout << "done" << std::endl;
}

MeshInterface*
MeshFactory::createNewMesh(std::string& meshName, mesh_t meshType, int nNodesReserve, int nElementsReserve)
{
  MeshInterface* mesh = getMesh(meshName);
  if ( mesh ) {
    std::cerr << "****** Mesh name collision" << std::endl;
    return mesh;
  }

  
  switch (meshType) {
    case NATIVE:
      mesh = new Mesh(nNodesReserve, nElementsReserve);
      break;
    case MOAB:
      mesh = new MoabMesh();
      break;
  }
  _meshMap[meshName] = mesh;
  return mesh;
}


MeshInterface*
MeshFactory::getMesh(std::string& meshName)
{
  std::map<std::string, MeshInterface*>::iterator _meshIt;
  _meshIt = _meshMap.find( meshName );
  if ( _meshIt == _meshMap.end() ) {
    return NULL;
  }
  else {
    return _meshIt->second;
  }

}

void
MeshFactory::deleteAllMesh()
{
  for (std::map<std::string, MeshInterface*>::iterator it=_meshMap.begin(); it!=_meshMap.end(); ++it)
    delete it->second;
}

void
MeshFactory::list()
{
  std::cout << std::endl << "Mesh generated:" << std::endl;
  for (std::map<std::string, MeshInterface*>::iterator it=_meshMap.begin(); it!=_meshMap.end(); ++it) {
    std::cout << "  " << it->first << std::endl;
  }
}
