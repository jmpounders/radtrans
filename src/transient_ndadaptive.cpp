#include "transient_ndadaptive.h"
#include "material.h"
#include "perfstats.h"
#include "log.h"

#include <cmath>
#include <iostream>

Transient_NDAdaptive::Transient_NDAdaptive(InputParser& input)
  : Transient(input), solnSecDer(NULL)
{
  std::vector<std::string> path(1,"Transient");
  std::vector<double> v;
  std::string s;

  v = input.getVector(path, "dtMin");
  if (v.size() == 0)
    _dtMin = 1.0e-5;
  else
    _dtMin = v[0];

  v = input.getVector(path, "dtMax");
  if (v.size() == 0)
    _dtMax = 1.0e-1;
  else
    _dtMax = v[0];

  v = input.getVector(path, "lteTol");
  if (v.size() == 0)
    _delta = 1.0e-4;
  else
    _delta = v[0];

  v = input.getVector(path, "relChange");
  if (v.size() == 0)
    _dtRelChange = 0.05;
  else
    _dtRelChange = v[0];

  LOG("LTE Tolerance = ", _delta);
  LOG("dtMin = ", _dtMin);
  LOG("dtMax = ", _dtMax);

  _output->registerStreamingOutput("dtRec");

  _pushBackSolution();
}

Transient_NDAdaptive::~Transient_NDAdaptive()
{
}

void
Transient_NDAdaptive::execute()
{
  PerfStats X("Transient_NDAdaptive::execute");
  
  LOG("Executing transient.");
  // Neutron speed
  int numGroups = _transportProblem->numGroups;

  // Solve the transient
  double t = _dt;
  double totalNeutronProduction, totalNeutronProductionPrev;
  double _dtNext;
  unsigned int numTimeSteps = 0;
  while (t < _tMax) {
    // Setup cross section modifications
    for (std::map<std::string, Material*>::iterator it=MaterialFactory::_materialMap.begin();
         it!=MaterialFactory::_materialMap.end();
         ++it) {
      Material* mat = it->second;
      double dnFissMod = 0;
      for (int i=0; i<numDelayedGroups; i++) {
        if (mat->beta[i] < 1.0e-10) break;
        dnFissMod += _dt * mat->lambda[i] * mat->beta[i] / (1.0 + _dt * mat->lambda[i]);
      }
      for (int g=0; g<numGroups; g++) {
        mat->_sigma_t_adder[g] = 1.0/(_dt * mat->_speed[g]);
        mat->_nu_sigma_f_scaling[g] = (1.0 - mat->betaEff + dnFissMod) / _criticalEigenvalue;
      }
    }

    // Solve current time step
    _calculateTransientSource();
    if (totalNeutronProductionPrev > 0.0) {
      double extrapFactor = log(totalNeutronProduction/totalNeutronProductionPrev);
      //solver->extrapolateSolution(extrapFactor, prevTimeSolution);
      solver->extrapolateSolution(extrapFactor);
    }
    solver->solve();
    _updateDelayedNeutronPrecs();
    _dtNext = _getTimeStep();
    _pushBackSolution();
    numTimeSteps++;

    // Time step output
    _output->streamOutput("time", t);
    _output->streamOutput("dtRec", _dtNext);

    // Calculate the total neutron production
    totalNeutronProductionPrev = totalNeutronProduction;
    totalNeutronProduction = solver->getTotalNeutronProduction() / _criticalEigenvalue;
    _output->streamOutput("totalNeutronProduction", totalNeutronProduction);

    // Setup next time step
    double relChange = (_dtNext - _dt)/_dt;
    if (relChange > _dtRelChange) _dtNext = _dt*(1.0 + copysign(_dtRelChange, relChange));
    if (_dtNext > _dtMax) _dtNext = _dtMax;
    if (_dtNext < _dtMin) _dtNext = _dtMin;
    if (numTimeSteps > 10) _dt = _dtNext;
    t += t+_dt>_tMax ? _tMax-t : _dt;
  }

  // Free up allocated memory
  delete [] delayedNeutronPrec;
  if (!prevTimeSolution) delete [] prevTimeSolution;

  _saveSolutionState();
}

void
Transient_NDAdaptive::_pushBackSolution()
{
  if (!prevPrevTimeSolution)
    prevPrevTimeSolution = new double [solver->getNumDOFs()];
  if (prevTimeSolution)
    for (int i=0; i<solver->getNumDOFs(); i++)
      prevPrevTimeSolution[i] = prevTimeSolution[i];
  
  prevTimeSolution = solver->copySolution(prevTimeSolution);
  _dtPrev = _dt;
}


void
Transient_NDAdaptive::_estimateSecondDerivative()
{
  if (!solnSecDer)
    solnSecDer = new double [solver->getNumDOFs()];

  for (int i=0; i<solver->getNumDOFs(); i++)
    solnSecDer[i] = *solver->getSolutionValue(i)/(_dt*_dt)
      - (_dt + _dtPrev)*prevTimeSolution[i]/(_dt*_dt*_dtPrev)
      + prevPrevTimeSolution[i]/(_dt*_dtPrev);
}
  
double
Transient_NDAdaptive::_getTimeStep()
{
  _estimateSecondDerivative();
  double _dtRec = 1.0e30;
  for (int i=0; i<solver->getNumDOFs(); i++) {
    _dtRec = fmin(_dtRec, sqrt(2.0*_delta*fabs(*solver->getSolutionValue(i)/solnSecDer[i])));
  }

  return _dtRec;

}
