#ifndef FIXEDSOURCE_H
#define FIXEDSOURCE_H

#include "solutionmanager.h"

/// Fixed source solution manager
/**
 *  Solve a fixed source problem.
 **/
class FixedSource : public SolutionManager
{
 public:
  FixedSource(InputParser& input);
  ~FixedSource();

  void execute();

 private:
  void _saveSolutionState();
};

#endif
