#ifndef SOLVERLOCALMOC_H
#define SOLVERLOCALMOC_H

#include "solverbase.h"
#include "moabmesh.h"


/// LocalMOC Solver
/**
 *  Still under active development
 */
class SolverDummy : public SolverBase
{
 public:
  SolverDummy( TransportProblem &tp );
  ~SolverDummy();
  void solve();
  double getScalarFlux(long space_i, int group_g);

 private:
  unsigned long _dofIndex(long i, int n, int g)
    { return _prob.quadOrder*_prob.numGroups*i + _prob.numGroups*n + g; };
  //{ return _prob.numCells*_prob.numGroups*n + _prob.numGroups*i + g; };
  void _mapDOFs();

  void _calculateSphericalQuadrature();
  std::vector<double> _mu;
  std::vector<double> _theta;
  std::vector<int> _negDir;

  MoabMesh* mesh;

  void _sweep(int n);
  double _getAngleFromVector(double dx, double dy);
  void _calculateMatrixAction(double* x, double* y);
};

#endif
