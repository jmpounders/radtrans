#include "transient.h"
#include "material.h"
#include "perfstats.h"
#include "log.h"

#include <cmath>

#include <iostream>


Transient::Transient(InputParser& input)
  : SolutionManager(SolutionManager::Transient, input),
    prevTimeSolution(NULL), prevPrevTimeSolution(NULL)
{
  std::vector<std::string> path(1,"Transient");
  std::vector<double> v;
  std::string s;

  v = input.getVector(path, "maxTime");
  _tMax = v[0];

  v = input.getVector(path, "dt");
  _dt = v[0];

  // Read initial conditions
  s = input.getString(path, "initialCondition");
  if (s == "zero") {
    solver->zeroSolution();
  }
  else {
    // Assume the input string contains an HDF file containing the IC solution data
    HDF5Interface hdf;
    hdf.open(s, 'R');
    
    HDFDataStruct<double> data;
    std::string var;

    var = "k";
    data.data = hdf.readData(s, var);
    _criticalEigenvalue = data[0];
    data.clear();
    LOG("Critical eigenvalue = ", _criticalEigenvalue);

    var = "solution";
    data.data = hdf.readData(s, var);
    
    if (data.data->dims[0] == solver->getNumDOFs()) {
      solver->setSolution(data.getData());
      solver->normalizeSolution(_criticalEigenvalue);
    }
    else {
      LOG_ERR("IC solution has inconsistent number of DOFs.");
    }
 
    hdf.close(s);
  }

  solver->sourceScaling = 1.0;

  if (_transportProblem->hasFixedSource)
    LOG_WARN("Fixed source in input is being ignored!!!");
  solver->sourceConfig.hasExternalSource = true;
  solver->sourceConfig.hasScatterSource = true;
  solver->sourceConfig.hasFissionSource = true;
  solver->sourceConfig.hasTransientSource = true;
  
  // Allocate space for delayed neutron precursors
  numDelayedGroups = 6;
  delayedNeutronPrec = new double [ _transportProblem->numCells*numDelayedGroups ];
  _calculateEquilibriumPrecs();

  // Setup output
  _output->registerStreamingOutput("time");
  _output->registerStreamingOutput("totalNeutronProduction");

  _output->streamOutput("time", 0.0);
  _output->streamOutput("totalNeutronProduction",
                        solver->getTotalNeutronProduction() / _criticalEigenvalue);

}

Transient::~Transient()
{
}


void
Transient::_pushBackSolution()
{  
  prevTimeSolution = solver->copySolution(prevTimeSolution);
}


void
Transient::_calculateEquilibriumPrecs()
{
  Material* mat;
  unsigned int index = 0;
  for (int i=0; i<_transportProblem->numCells; i++) {
    double productionRate = solver->getNeutronProduction(i) / _criticalEigenvalue;
    mat = _transportProblem->mesh->getElementMat(i);
    for (int d=0; d<numDelayedGroups; d++) {
      double beta = mat->beta[d];
      double lambda = mat->lambda[d];
      if (beta < 1.0e-10) {
        delayedNeutronPrec[index] = 0.0;
      }
      else {
        delayedNeutronPrec[index] = beta/lambda*productionRate;
      }
      index++;
    }
  }
}

void
Transient::_updateDelayedNeutronPrecs()
{
  Material* mat;
  unsigned int index = 0;
  for (int i=0; i<_transportProblem->numCells; i++) {
    double productionRate = solver->getNeutronProduction(i) / _criticalEigenvalue;
    mat = _transportProblem->mesh->getElementMat(i);
    for (int d=0; d<numDelayedGroups; d++) {
      double beta = mat->beta[d];
      double lambda = mat->lambda[d];
      if (beta > 1.0e-10) {
        delayedNeutronPrec[index] = (delayedNeutronPrec[index] + _dt*beta*productionRate)/(1.0 + _dt*lambda);
      }
      index++;
    }
  }
}


void
Transient::_saveSolutionState()
{
  solver->writeSolution(_output);
  solver->writeScalarFlux(_output);
}


void
Transient::_calculateTransientSource()
{
  Material* mat;
  double src, dnSrc, chi;
  unsigned int index = 0;
  for (int i=0; i<_transportProblem->numCells; i++) {
    mat = _transportProblem->mesh->getElementMat(i);
    dnSrc = 0;
    for (int d=0; d<numDelayedGroups; d++) {
      if (mat->beta[d] > 1.0e-10) {
        dnSrc += mat->lambda[d] / (1.0 + _dt * mat->lambda[d]) * delayedNeutronPrec[index];
      }
      index++;
    }
    for (int g=0; g<_transportProblem->numGroups; g++) {
      chi = *_transportProblem->mesh->getElementMat(i)->getFissionSpectrum(g+1);
      for (int n=0; n<_transportProblem->quadOrder; n++) {
        src = chi*dnSrc + *solver->getSolutionValue(i, n, g)/(_dt * mat->_speed[g]);
        _transportProblem->putExtSource(i, n, g, src);
      }
    }
  }
}

