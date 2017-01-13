#include "transient_uts.h"
#include "material.h"
#include "perfstats.h"
#include "log.h"

Transient_UTS::Transient_UTS(InputParser& input)
  : Transient(input)
{
  _pushBackSolution();
}

Transient_UTS::~Transient_UTS()
{
}

void
Transient_UTS::execute()
{
  PerfStats X("Transient_UTS::execute");
  
  LOG("Executing transient.");
  // Neutron speed
  int numGroups = _transportProblem->numGroups;

  // Solve the transient
  double t = _dt;
  double totalNeutronProduction;
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
    solver->solve();
    _updateDelayedNeutronPrecs();
    _dtNext = _getTimeStep();
    _pushBackSolution();
    numTimeSteps++;

    // Time step output
    _output->streamOutput("time", t);

    // Calculate the total neutron production
    totalNeutronProduction = solver->getTotalNeutronProduction() / _criticalEigenvalue;
    _output->streamOutput("totalNeutronProduction", totalNeutronProduction);

    // Setup next time step
    _dt = _dtNext;
    t += t+_dt>_tMax ? _tMax-t : _dt;
  }

  // Free up allocated memory
  delete [] delayedNeutronPrec;
  if (!prevTimeSolution) delete [] prevTimeSolution;

  _saveSolutionState();
}

double
Transient_UTS::_getTimeStep()
{
  return _dt;
}
