#include "solutionmanager.h"
#include "solverlocalmoc.h"
#include "solverregmoc.h"
#include "timing.h"
#include "material.h"
#include "mesh.h"
#include "transportproblem.h"
#include "log.h"

#include <vector>
#include <string>

SolutionManager::SolutionManager(ProblemType problemType, InputParser& input)
  : problemType(problemType), _input(input)
{
  // Make materials
  int numGroups = materialFactory.createMaterialsFromInput(input);
  materialFactory.list();

  // Make mesh
  // This shouldn't need the materialFactory.  temp fix for mat assignments
  meshFactory.createMeshFromInput(input, materialFactory);
  meshFactory.list();

  // Form the transport problem
  _transportProblem = new TransportProblem(input, meshFactory, numGroups);

  // Setup the output
  _output = new Output("output", Output::ASCII | Output::HDF5 | Output::MESH);

  
  // Parse solver input
  std::vector<std::string> path(1,"Solver");
  std::string solverType = _input.getString(path, "type");
  LOG("Solver type = " + solverType);

  // Set solver pointer
  if (solverType == "localMOC")
    solver = new SolverLocalMOC(*_transportProblem);
  else if (solverType == "regMOC")
    solver = new SolverRegMOC(*_transportProblem);
  else
    LOG_ERR("Invalid solver type");
  
  std::vector<double> v;
  v = _input.getVector(path, "rInfTolerance");
  if (v.size() > 0)
    solver->setRInfTol( v[0] );

  v = _input.getVector(path, "maxIters");
  if (v.size() > 0)
    solver->setMaxIters( v[0] );
  
}

SolutionManager::~SolutionManager()
{
  delete solver;
  delete _output;
  delete _transportProblem;
  meshFactory.deleteAllMesh();
}

void
SolutionManager::_dumpBaseSolution(std::string fileName)
{
  Output* stateFile = new Output("state", Output::HDF5);
}
