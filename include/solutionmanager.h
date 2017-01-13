#ifndef SOLUTIONMANAGER_H
#define SOLUTIONMANAGER_H

#include "solverbase.h"
#include "inputparser.h"
#include "output.h"

/// Solution manager object
/**
 *  This is an abstract class that manages the setup and solution for particular applications
 *  (e.g., fixed source, eigenvalue, and transient calculations).  Each particular
 *  application must inherit from this class.  This class also creates all of the
 *  necessary sub-classes from input, i.e., materials, mesh, problem parameters, 
 *  and solver.
 **/
class SolutionManager
{
 public:
  enum ProblemType {FixedSource, Eigenvalue, Transient} problemType;
  
  virtual ~SolutionManager();
   
  virtual void execute() = 0;
    
 protected:
  SolutionManager(ProblemType problemType, InputParser& input);
  virtual void _saveSolutionState() = 0;

  MaterialFactory materialFactory;
  MeshFactory meshFactory;

  TransportProblem* _transportProblem;

  SolverBase* solver;
  InputParser& _input;
  Output* _output;

 private:
  void _dumpBaseSolution(std::string fileName);
};

#endif
