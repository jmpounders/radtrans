#ifndef SOLVERBASE_H
#define SOLVERBASE_H

#include<vector>
#include<map>

#include "transportproblem.h"
#include "dofobj.h"
#include "timing.h"
#include "output.h"

/// Source configuration description
/**
 *  This struct specifies which sources are "on" for a given calculation
 **/
struct SourceConfiguration
{
  bool hasExternalSource;
  bool hasScatterSource;
  bool hasFissionSource;
  bool hasTransientSource;
};


/// Abstract solver class
/**
 *  DOFlist provides an ordering of the DOFs
 *  DOFmap provides an easy way to look up the value associated with a DOF
 **/
class SolverBase
{
 public:
  virtual ~SolverBase();

  virtual void solve() = 0;

  // Error/convergence calculators
  double calculateL2SolutionNorm(double* solutionComp);
  double calculateLInfSolutionNorm(double* solutionComp);
  double calculateRInfSolutionNorm(double* solutionComp);
  double calculateL2SolutionError();
  double calculateRInfSolutionError();

  // Getters
  const double* getSolutionValue(long space_i, int quad_n, int group_g)
    { return DOFmap[ DOFObj(space_i, quad_n, group_g) ]; }
  const double* getSolutionValue(long i)
    { return &_solution[i]; };
  long getNumDOFs() { return _numDOF; };
  TransportProblem& getTransportProblem() { return _prob; };
  virtual double getScalarFlux(long space_i, int group_g) = 0;
  double getNeutronProduction(long space_i);
  double getTotalNeutronProduction();

  // Setters
  void setRInfTol(double tol) { _convRInfTol = tol; };
  void setMaxIters(int maxIters) { _maxIters = maxIters; };

  // Utilitiy functions
  void setSolution(double* solution);
  double* copySolution(double* solutionCopy);
  void normalizeSolution(double totalProductionRate);
  void zeroSolution();
  void extrapolateSolution(double expExtrapFactor);
  void extrapolateSolution(double expExtrapFactor, double* baseSolution);
  void writeSolution(Output* outputFile);
  void writeScalarFlux(Output* outputFile);

  void printIterStatus(std::string name, int iter, double err, double tol);

  SourceConfiguration sourceConfig;

  double sourceScaling;       //!< This scales the external source, but is used as a fission source in the power method
  double criticalEigenvalue;  //!< "Critical" eigenvalue for transient calculations based on ICs from an eigenvalue calculation

 protected:
  SolverBase(TransportProblem &tp);
  void _calculateResidual();
  virtual void _mapDOFs() = 0;
  virtual void _calculateMatrixAction(double* x, double* y) = 0;

  void _saveOldSolution();

  TransportProblem &_prob;                //!< Reference to the base transport problem

  long _numDOF;                            //!< Number of DOF
  long _numSpaceDOF;                       //!< Number of spatial DOF
  std::vector<DOFObj> DOFlist;            //!< List of DOF (currently unused)
  std::map<DOFObj, double*> DOFmap;       //!< DOF map

  double *_solution;                      //!< Solution vector
  double *_residual;                      //!< Residual vector
  double *_source;                        //!< Source vector
  double *_h;                             //!< Mesh spacing

  double *_solutionPrev;

  double _convRInfTol;
  int _maxIters;

  std::map<DOFObj, double> _bdryFlux;

  
};

#endif
