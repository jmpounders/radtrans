#include "transportproblem.h"
#include "log.h"

#include "hdf5interface.h"

TransportProblem::TransportProblem(InputParser& input, MeshFactory& meshFactory, int numGroups_in) :
  numGroups(numGroups_in), scatterAnisotropy(0)
{
  // Make problem
  std::cout << "Making problem... ";
  std::vector<std::string> path(1,"Transport");
  std::vector<double> v;
  std::string genString;

  std::string meshName = input.getString(path, "mesh");
  if (meshName == "empty") {
    LOG_ERR("ERROR: must input a mesh");
  }
  mesh = meshFactory.getMesh(meshName);
  numCells = mesh->numElements();
  numEdges = mesh->numEdges();
  numNodes = mesh->numNodes();

  std::string boundaryConditions = input.getString(path, "boundary");
  if (boundaryConditions == "empty")
    globalBC = vacuum;
  else if (boundaryConditions == "vacuum")
    globalBC = vacuum;
  else if (boundaryConditions == "reflecting")
    globalBC = reflecting;
  else if (boundaryConditions == "source")
    globalBC = source;
  else
    LOG_ERR("Unknown boundary condition");



  v = input.getVector(path, "quadOrder");
  quadOrder = int(v[0]);

  if (quadOrder > 0) {
    omega_x = input.getVector(path, "omega_x");
    if ( omega_x.size() != quadOrder)
      LOG_ERR("Number of angular quadrature points does not match the quadrature order.");

    omega_y = input.getVector(path, "omega_y");
    if ( omega_y.size() == 0 )
      omega_y.resize(quadOrder, 0.0);
    if ( omega_y.size() != quadOrder)
      LOG_ERR("Number of angular quadrature points does not match the quadrature order.");

    omega_z = input.getVector(path, "omega_z");
    if ( omega_z.size() == 0 )
      omega_z.resize(quadOrder, 0.0);
    if ( omega_z.size() != quadOrder)
      LOG_ERR("Number of angular quadrature points does not match the quadrature order.");
    
    weights = input.getVector(path, "weights");
    if ( weights.size() != quadOrder)
      LOG_ERR("Number of angular quadrature weights does not match the quadrature order.");
  }
  
  bcLeft = vacuum;
  bcRight = vacuum;

  // Always allocate space, because this is used by eigenvalue and transient
  extSource = new double [ numCells*quadOrder*numGroups ];
  for (long i=0; i<numCells*quadOrder*numGroups; i++)  {
    extSource[i] = 0.0;
  }

  // Check for fixed source in input
  path.push_back("ExternalSource");
  if (input.countSets(path) == 1) {
    hasFixedSource = true;

    // Assign source
    genString = input.getString(path, "type");
    // Uniform source
    if (genString == "uniform") {
      v = input.getVector(path, "magnitude");
      double srcMagnitude = v[0];
      for (int g=0; g<numGroups; g++) {
        for (int i=0; i<numCells; i++) {
          for (int n=0; n<quadOrder; n++) {
            putExtSource(i,n,g, srcMagnitude);
          }
        }
      }
    }
    else if (genString == "file") {
      std::string fileName = input.getString(path,"name");
      std::string varName("/source");
      double* srcArray;
      HDF5Interface h5;
      h5.open(fileName,'R');
      HDFData* sourceData = h5.readData(fileName, varName);
      srcArray = (double*)sourceData->data;
      for (int g=0; g<numGroups; g++) {
        for (int i=0; i<numCells; i++) {
          for (int n=0; n<quadOrder; n++) {
            putExtSource(i,n,g, srcArray[i*quadOrder + n]);
          }
        }
      }
      delete [] (double*)sourceData->data;
    }
    else if (genString == "boundary") {
      nxnyq = input.getVector(path, "nxnyq");
    }
    else {
      LOG_ERR("Invalid source type: ", genString);
    }
  }
  else {
    hasFixedSource = false;
  }
  
  // Allocate space for reference solution
  refSolution = new double [ numCells*quadOrder*numGroups ];
  for (long i=0; i<numCells*quadOrder*numGroups; i++)  {
    refSolution[i] = 1.0e30;
  }
  
  std::cout << "done" << std::endl;

}

TransportProblem::~TransportProblem()
{
  delete [] extSource;
  delete [] refSolution;
}
