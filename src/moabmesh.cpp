#include "moabmesh.h"
#include "moab/AdaptiveKDTree.hpp"
#include "global.h"

MoabMesh::MoabMesh()
  : MeshInterface()
{
  // Instantiate a new MOAB interface
  _mb = new moab::Core;
}

MoabMesh::~MoabMesh()
{
  // Delete the MOAB interface
  logMemoryUse();
  delete _mb;
  _mb = NULL;
  _materialList.clear();
}


void
MoabMesh::loadMesh(std::string inputFileName)
{
  _rval = _mb->load_file(inputFileName.c_str());


  moab::Range tris, edges, verts;
  int cnt, vertsPerEnt;
  _rval = _mb->get_entities_by_dimension(0, 0, verts);
  _rval = _mb->coords_iterate(verts.begin(), verts.end(), _x, _y, _z, cnt);
  
  _rval = _mb->get_entities_by_dimension(0, 2, tris);
  _rval = _mb->connect_iterate(tris.begin(), tris.end(), _connectivity, vertsPerEnt, cnt);

  _dimension = 2;
  _numNodes = verts.size();
  _numElements = tris.size();
  _rval = _mb->get_entities_by_dimension(0, 2, _tris);

  _rval = _mb->get_entities_by_dimension(0, 0, _verts);

  //moab::Range vertAdj;
  //_rval = _mb->get_adjacencies(_verts, 1, true, vertAdj, moab::Interface::UNION);
  // This should be " - 2" on RHS.  Overcounting nodes or elements?
  _numEdges = _numNodes + _numElements - 1;
  _identifyEdges();
  _identifyBoundaryElements();

  LOG("MOAB mesh read from file.");
  LOG("  Nodes:   ", _numNodes);
  LOG("  Edges:   ", _numEdges);
  LOG("  Elements:", _numElements);
  logMemoryUse();
  
  _materialList.reserve(_numElements);

  //double point[3] = {5.0, 5.9, 0.0};
  //queryPointLocation( point );
}

void
MoabMesh::getCurrentElementFromID(long elementID, UltraLightElement& elementUL)
{
  elementUL.elementID = elementID;

  // Setup the current element
  int nxtNghbr[3] = {1, 2, 0};
  long nghbrID;
  
  for (int v=0; v<3; v++) {
    long vertexID = _connectivity[elementID*3 + v]-1;
    long nextVertexID = _connectivity[elementID*3 + nxtNghbr[v]]-1;
    elementUL.x[v] = _x[vertexID];
    elementUL.y[v] = _y[vertexID];
    elementUL.z[v] = _z[vertexID];
    elementUL.vertexID[v] =  vertexID;
    nghbrID = getAdjacentTriangle(elementID,vertexID,nextVertexID);
    if (nghbrID < 0.0) nghbrID = -(v+1);
    elementUL.neighborID[v] = nghbrID;
  }
}


void
MoabMesh::_identifyEdges()
{
  for (std::vector<moab::EntityHandle>::iterator it = _tris.begin(); it != _tris.end(); ++it) {
    long elementID = _mb->id_from_handle(*it)-1;
    int nxtNghbr[3] = {1, 2, 0};
    for (int v=0; v<3; v++) {
      // Get the two verts bounding an edge
      moab::Range _boundingVerts;
      _boundingVerts.insert( _connectivity[elementID*3 + v] );
      _boundingVerts.insert( _connectivity[elementID*3 + nxtNghbr[v]] );

      // Get the tri on the other side of the edge
      moab::Range _neighbor;
      _rval = _mb->get_adjacencies(_boundingVerts, 2, false, _neighbor);
      _neighbor.erase( *it );
      long nghbrID =  _mb->id_from_handle(_neighbor[0])-1;
      if (nghbrID<0) nghbrID = -(v+1);

      // Add this unique pair to the map
      std::pair<long,long> nghbrPair =
        elementID>nghbrID ? std::make_pair(elementID,nghbrID) : std::make_pair(nghbrID,elementID);
      _edgeID[nghbrPair] = 0;

      _boundingVerts.clear();
      _neighbor.clear();
    }
  }

  // Number the surfaces in the order in which they were added
  long edgeCntrID = 0;
  for (std::map<std::pair<long, long>, long>::iterator it = _edgeID.begin(); it != _edgeID.end(); ++it) {
    it->second = edgeCntrID;
    edgeCntrID++;
  }
  LOG_DBG("identified ", _edgeID.size(), " edges.");
}

void
MoabMesh::_identifyBoundaryElements()
{
  if (_edgeID.size() == 0) _identifyEdges();

  for (std::vector<moab::EntityHandle>::iterator it = _tris.begin(); it != _tris.end(); ++it) {
    long elementID = _mb->id_from_handle(*it)-1; 
    if (_edgeID.count(std::make_pair(elementID,-1)) > 0
        || _edgeID.count(std::make_pair(elementID,-2)) > 0
        || _edgeID.count(std::make_pair(elementID,-3)) > 0)
      boundaryElements.push_back(elementID);
  }

  boundaryElements.resize(boundaryElements.size());
  LOG_DBG("identified ", boundaryElements.size(), " boundary elements.");
}

long
MoabMesh::getEdgeID(long node1, long node2)
{
  // Make the (node1, node2) pair
  std::pair<long,long> nghbrPair =
    node1>node2 ? std::make_pair(node1,node2) : std::make_pair(node2,node1);
  // Find the edge by inserting it.  This is a built in check to verify that it is actually there.
  // Otherwise, you would have to do a "find" to verify existence followed by a look up.
  std::pair<std::map<std::pair<long, long>, long>::iterator, bool> it = _edgeID.insert( std::pair<std::pair<long,long>, long>(nghbrPair, 0));
  if (it.second) // It wasn't there
    return -1;
  else           // It was there, and here it is
    return it.first->second;
}


void
MoabMesh::readMeshSweepOrder(const std::string fileName)
{
  LOG("Reading sweep order from file.");
  HDF5Interface hdf;
  hdf.open(fileName, 'R');
  HDFDataStruct<int> data;
  data.data = hdf.readData(fileName, std::string("sweeporder"));
  hdf.close(fileName);

  LOG_DBG("Sweep data is ", data.data->dims[0], " by ", data.data->dims[1]);

  long nCells = data.data->dims[0];
  int nAngles = data.data->dims[1];

  _orderedMeshSet.reserve(nAngles);
  _orderedMeshSet.resize(nAngles);

  // Create sweep ordering array (need to make sure this gets deallocated somehwere)
  for (int j=0; j<nAngles; j++) {
    _mb->create_meshset(moab::MESHSET_ORDERED, _orderedMeshSet[j]);
    for (long i=0; i<nCells; i++) {
      moab::EntityHandle entity;
      _mb->handle_from_id(moab::MBTRI, data(i,j)+1, entity);
      _mb->add_entities(_orderedMeshSet[j], &entity, 1);
    }
  }
  data.clear();
}


long
MoabMesh::queryPointLocation( double* point)
{
  // Get all 2-dimensional entities in the databse
  moab::Range tris;
  _rval = _mb->get_entities_by_dimension(0, 2, tris);
  
  // Build a KD tree to do the search
  moab::EntityHandle treeRoot = 0;
  moab::AdaptiveKDTree kdTree(_mb);
  _rval = kdTree.build_tree(tris, &treeRoot);

  // Find the point in a tree leaf (each leaf contains multiple elements)
  moab::AdaptiveKDTreeIter treeIter;
  LOG_DBG("Point = ",point[0], " ", point[1]);
  _rval = kdTree.point_search(point, treeIter);

  // Get the entities out of the leaf
  moab::Range leafTris;
  _rval = _mb->get_entities_by_dimension(treeIter.handle(), 2, leafTris, false);
  LOG_DBG("Found ", leafTris.size(), " triangles containing the point.");

  // Search through the entities to find the one containing the point
  moab::Range::iterator rit;
  for (rit = leafTris.begin(); rit != leafTris.end(); rit++) {
    moab::EntityHandle this_ent = *rit;

    const moab::EntityHandle* conn;
    int numConnected;
    _rval = _mb->get_connectivity(this_ent, conn, numConnected);

    double coords[numConnected*3];
    _rval = _mb->get_coords(conn, numConnected, coords);


    double v[2] = {point[0], point[1]};
    double v0[2] = {coords[0], coords[1]};
    double v1[2] = {coords[3]-v0[0], coords[4]-v0[1]};
    double v2[2] = {coords[6]-v0[0], coords[7]-v0[1]};
    double a = ((v[0]*v2[1]-v[1]*v2[0]) - (v0[0]*v2[1]-v0[1]*v2[0]))/(v1[0]*v2[1]-v1[1]*v2[0]);
    double b = -((v[0]*v1[1]-v[1]*v1[0]) - (v0[0]*v1[1]-v0[1]*v1[0]))/(v1[0]*v2[1]-v1[1]*v2[0]);
    if (a>0 && b>0 && a+b<1.0) {
      LOG_DBG(_mb->id_from_handle(this_ent));
      for (int v=0; v<numConnected; v++) {
        LOG_DBG(coords[3*v], " ",coords[3*v+1]," ", coords[3*v+2]);
      }
      return _mb->id_from_handle(this_ent);
    }
  }
  return 0;

}

void
MoabMesh::writeMesh(Output* outputFile)
{
  std::string fileName = outputFile->getName();
  fileName.append(".vtk");

  moab::Range fileEntities;
  fileEntities.insert( _verts.begin(), _verts.end());

  moab::Range tris;
  _rval = _mb->get_entities_by_dimension(0, 2, tris);
  fileEntities.insert( tris.begin(), tris.end());
  
  //_rval = _mb->write_file(fileName.c_str(), 0,0, &fileEntities.front(), fileEntities.size());
  _rval = _mb->write_file(fileName.c_str());
}


void
MoabMesh::tagMesh(std::string tagName, double* tagDataBuffer, long tagSize)
{
  moab::Tag tag;
  double defaultValue = 0.0;
  _rval = _mb->tag_get_handle(tagName.c_str(), 1, moab::MB_TYPE_DOUBLE, tag,
                                           moab::MB_TAG_DENSE | moab::MB_TAG_CREAT,
                                           &defaultValue);

  if (tagSize == _tris.size()) {
    LOG_DBG("Tagging mesh, size = ", _tris.size());
    _rval = _mb->tag_set_data(tag, &_tris.front(), _tris.size(), tagDataBuffer);
  }
  else if (tagSize == _verts.size()) {
    _rval = _mb->tag_set_data(tag, _verts, tagDataBuffer);
  }
  else {
    LOG_ERR("Could not tag mesh");
    LOG_ERR("  Received input buffer with size ", tagSize);
    LOG_ERR("  There are ", _tris.size(), " triangles and ", _verts.size(), " vertices.");
  }
}


void
MoabMesh::logMemoryUse()
{
  unsigned long long totalStorage, totalAmortStorage;
  unsigned long long entityStorage, entityAmortStorage;
  unsigned long long adjacencyStorage, adjacencyAmortStorage;
  unsigned long long tagStorage, tagAmortStorage;

  _mb->estimated_memory_use(NULL, 0,
                            &totalStorage, &totalAmortStorage,
                            &entityStorage, &entityAmortStorage,
                            &adjacencyStorage, &adjacencyAmortStorage,
                            NULL, 0,
                            &tagStorage, &tagAmortStorage);
  LOG("MOAB Mesh Memory Usage (bytes)");
  LOG("  Entity:      ", entityStorage, "/", entityAmortStorage);
  LOG("  Adjacency:   ", adjacencyStorage, "/", adjacencyAmortStorage);
  LOG("  Tag:         ", tagStorage, "/", tagAmortStorage);
  LOG("  Total:       ", totalStorage, "/", totalAmortStorage);
}



long
MoabMesh::getAdjacentTriangle( long triangle, long vertex1, long vertex2 )
{
  moab::Range _boundingVerts;
  _boundingVerts.insert( vertex1+1 );
  _boundingVerts.insert( vertex2+1 );

  // Using the bounding verts to get neighbor data uses more memory than just getting all adjacent tris
  moab::Range _neighbor;
  moab::EntityHandle triangleHandle;
  _rval = _mb->handle_from_id(moab::MBTRI, triangle+1, triangleHandle);
  _rval = _mb->get_adjacencies(_boundingVerts, 2, true, _neighbor);
  _neighbor.erase( triangleHandle );
  return _mb->id_from_handle(_neighbor[0]) - 1;
}

double
MoabMesh::getElementVolume(long elementID)
{
  double x[3], y[3];
  for (int v=0; v<3; v++) {
    long vertexID = _connectivity[elementID*3 + v]-1;
    x[v] = _x[vertexID];
    y[v] = _y[vertexID];
  }
  return std::abs(x[0]*(y[1]-y[2])
                  + x[1]*(y[2]-y[0])
                  + x[2]*(y[0]-y[1]))/2.0;
}
