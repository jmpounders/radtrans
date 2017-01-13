#ifndef TRANSIENT_H
#define TRANSIENT_H

#include "solutionmanager.h"

/// Transient solution manager
/**
 *  Solve a transient problem
 **/
class Transient : public SolutionManager
{
 public:
  Transient(InputParser& input);
  virtual ~Transient();

 protected:
  void _pushBackSolution();
  void _calculateEquilibriumPrecs();
  void _updateDelayedNeutronPrecs();
  void _saveSolutionState();
  void _calculateTransientSource();

  virtual double _getTimeStep() = 0;
  
  double* prevTimeSolution;
  double* prevPrevTimeSolution;

  double _tMax;
  double _dt;
  double _criticalEigenvalue;

  int numDelayedGroups;
  double* delayedNeutronPrec;

};

#endif
