#ifndef TRANSIENT_UTS_H
#define TRANSIENT_UTS_H

#include "transient.h"

/// Transient solution manager with unform time steps
class Transient_UTS : public Transient
{
 public:
  Transient_UTS(InputParser& input);
  ~Transient_UTS();

  void execute();

 private:
  double _getTimeStep();
};

#endif
