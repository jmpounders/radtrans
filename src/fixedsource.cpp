#include "fixedsource.h"
#include "log.h"


FixedSource::FixedSource(InputParser& input)
  : SolutionManager(SolutionManager::FixedSource, input)
{
  solver->sourceConfig.hasExternalSource = true;
  solver->sourceConfig.hasScatterSource = true;
  solver->sourceConfig.hasFissionSource = true;
  solver->sourceConfig.hasTransientSource = false;
}

FixedSource::~FixedSource()
{
}

void
FixedSource::execute()
{
  if (!_transportProblem->hasFixedSource) {
    LOG_ERR("No Fixed Source!!!!!");
    return;
  }
  
  solver->solve();

  _saveSolutionState();
}

void
FixedSource::_saveSolutionState()
{
  solver->writeSolution(_output);
  solver->writeScalarFlux(_output);
}
