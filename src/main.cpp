/*! \mainpage Title
 *  \section devnotes Development Notes:
 *   * Could add an extra layer of abstraction in the solver chain, i.e. SolverSIBase
 *
 *  \section todo To Do:
 *   * Improve display of output logging/output file
 *   * Add delayed neutron spectrum
 *   * Handle errors better (gracefully).
 *   * Add comments to input file.
 *   * Make a virtual factory to inherit from
 *   * Make an hdf writer of a dataset.
 *  
 */

#include <iostream>
#include <string>

#include "timing.h"
#include "inputparser.h"
#include "dataset.h"
#include "transient_ndadaptive.h"
#include "transient_uts.h"
#include "fixedsource.h"
#include "perfstats.h"
#include "log.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
  std::cout << "RADIATION TRANSPORT CODE" << std::endl;
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " input.file" << std::endl;
    return 0;
  }

#ifdef _OPENMP
  int maxNumThreads = omp_get_max_threads();
  LOG("OpenMP has been detected. The max number of threads is ", maxNumThreads, ".");
#endif
  
  ////////////////////////////////////////////////////////////////////
  // Read input
  ////////////////////////////////////////////////////////////////////
  LOG("Parsing input from file ", argv[1]);
  Timing::timer.addTimer("Input Parsing");
  Timing::timer.punch("Input Parsing");

  std::string inputFile(argv[1]);

  InputParser input;
  bool parseResult = input.parseInputFile(inputFile);
  if (!parseResult) {
    LOG_ERR("Parse failure.");
    input.echoInput();
    return 1;
  }
  input.echoInput();
  Timing::timer.punch("Input Parsing");


  ////////////////////////////////////////////////////////////////////
  // Set up problem and solver
  ////////////////////////////////////////////////////////////////////
  LOG("Setting up problem.");
  Timing::timer.addTimer("Problem Setup");
  Timing::timer.punch("Problem Setup");

  SolutionManager* solution;
  
  std::vector<std::string> path(1,"Transient");
  if ( input.countSets(path = std::vector<std::string>(1,"Transient")) > 0) {
    LOG("Transient problem detected.");
    std::string transientType = input.getString(path, "type");
    if (transientType == "uts") 
      solution = new Transient_UTS(input);
    else if (transientType == "ndadaptive")
      solution = new Transient_NDAdaptive(input);
  }
  else if (input.countSets(path = std::vector<std::string>(1,"FixedSource")) > 0) {
    LOG("Fixed source problem detected.");
    solution = new FixedSource(input);
  }
  else {
    LOG_ERR("No suitable solution manager found.");
    return 1;
  }

  Timing::timer.punch("Problem Setup");

  ////////////////////////////////////////////////////////////////////
  // Solve problem
  ////////////////////////////////////////////////////////////////////
  LOG("Starting problem solution.");
  Timing::timer.addTimer("Problem Solution");
  Timing::timer.punch("Problem Solution");
  solution->execute();
  Timing::timer.punch("Problem Solution");



  
  delete solution;
  Timing::timer.display();
  PerfStats::display();

  LOG("Complete.");
}



