#ifndef SOLVERLOCALMOC_H
#define SOLVERLOCALMOC_H

#include "solverbase.h"
#include "moabmesh.h"


/// LocalMOC Solver
/**
 *  Still under active development
 */
class SolverLocalMOC : public SolverBase
{
 public:
  SolverLocalMOC( TransportProblem &tp );
  ~SolverLocalMOC();
  void solve();
  double* getResidual();
  double getScalarFlux(long space_i, int group_g);

 private:
  unsigned long _dofIndex(long i, int n, int g, int edgeLoc=0)
    { return 2*_prob.quadOrder*_prob.numGroups*i + 2*_prob.numGroups*n + 2*g + edgeLoc; };
  void _mapDOFs();

  void _calculateSphericalQuadrature();
  int _getReflectedDirection(int i, double nx, double ny, double& OmegaDotn);
  std::vector<double> _mu;
  std::vector<double> _theta;
  std::vector<int> _negDir;

  MoabMesh* mesh;

  double* _cellFlux;
  double* _surfacePosition;

  void _sweep(int n);
  void _getTriangleOrientation(UltraLightElement &element2, int n,
                               double &mu01, double &mu12, double &mu20,
                               double &pathDist, int &evDir, int &veDir,
                               long &edgeNeighbor, long &vertexNeighbor1, long &vertexNeighbor2,
                               int &xedge, int &v0, int &v1, int &v2, double &surfacePosition);
  double _getAngleFromVector(double dx, double dy);

  void _calculateSource(double* fissionFlux);
  void _getExternalSource();
  void _addScatterSource();
  void _calculateMatrixAction(double* x, double* y);

  void _applyBoundaryConditions();

  unsigned int _numPhaseSpaceDOF;
};

#endif
