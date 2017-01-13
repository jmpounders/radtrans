#ifndef TRANSIENT_NDADAPATIVE_H
#define TRANSIENT_NDADAPATIVE_H

#include "transient.h"

/// Transient solution manager with nested difference adaptive time stepping
class Transient_NDAdaptive : public Transient
{
 public:
  Transient_NDAdaptive(InputParser& input);
  ~Transient_NDAdaptive();

  void execute();

 private:
  void _pushBackSolution();
  void _estimateSecondDerivative();
  double _getTimeStep();

  double* solnSecDer;
  double _delta;
  double _dtMin;
  double _dtMax;
  double _dtPrev;
  double _dtRelChange;
};

#endif
