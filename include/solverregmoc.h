#ifndef SOLVERREGMOC_H
#define SOLVERREGMOC_H

#include "solverbase.h"
#include "moabmesh.h"

// This needs to be updated to four subcells
struct TriangleDescriptorReg
{
  long edgeNeighbor, vertexNeighbor1, vertexNeighbor2;
  // subcell intercell boundary fluxes
  double psi01, psi13, psi21;
  // triangle boundary fluxes
  double psi0,psi1,psi2,psi3,psi4,psi5;
  // subcell fluxes
  double cell0,cell1,cell2,cell3;
  double sigma;
  double edge01[2], edge12[2], edge20[2];
  double theta0, theta1, theta2;
  double d01, d12, d20;
  int i0, i1, i2, i3;
};



/// LocalMOC Solver
/**
 *  Still under active development
 */
class SolverRegMOC : public SolverBase
{
 public:
  SolverRegMOC( TransportProblem &tp );
  ~SolverRegMOC();
  void solve();
  double* getResidual();
  double getScalarFlux(long space_i, int group_g);

 private:
  long _dofIndex(long i, int n, int g, int edgeLoc=0)
    { return 2*_prob.quadOrder*_prob.numGroups*i + 2*_prob.numGroups*n + 2*g + edgeLoc; };
  long _dofIndexPS(long i, int n, int g, int subCell=0)
    { return 4*_prob.quadOrder*_prob.numGroups*i + 4*_prob.numGroups*n + 4*g + subCell; };
  void _mapDOFs();

  void _calculateSphericalQuadrature();
  int _getReflectedDirection(int i, double nx, double ny, double& OmegaDotn);
  std::vector<double> _mu;
  std::vector<double> _theta;
  std::vector<int> _negDir;

  MoabMesh* mesh;

  double* _cellFlux;

  void _sweep(int n);
  void _getTriangleOrientation(UltraLightElement &element2, int& n,
                               TriangleDescriptorReg& tri, int& blockID,
                               double& theta, int &v0, int &v1, int &v2);
  void _applyBoundaryConditions();
  double _getAngleFromVector(double dx, double dy);

  void _calculateSource(double* fissionFlux);
  void _getExternalSource();
  void _addScatterSource();
  void _calculateMatrixAction(double* x, double* y);
  double getSubCellScalarFlux(long space_i, int group_g, int subCell);

  void triangleSolveA(TriangleDescriptorReg& tri, double& phi, double& theta, double& mu, long& elementID, int& n, int& g);
  void triangleSolveB(TriangleDescriptorReg& tri, double& phi, double& theta, double& mu, long& elementID, int& n, int& g);

  void triangleSolveEV(double& q, double& S, double& sigma, double& x,
                       double& psi1, double& psi2,
                       double& psi20, double& psi01, double& cell, double& psiv);
  void triangleSolveVE2(double& q, double& S, double& sigma,
                        double& w1, double& w2, double& w3,
                        double& psi0r, double& psi0l,
                        double& psi1, double& psi2,
                        double& psi12, double& cell);


  unsigned int _numPhaseSpaceDOF;
};


#endif
