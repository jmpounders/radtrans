#include "solverregmoc.h"
#include "global.h"
#include "sweeper.h"

SolverRegMOC::SolverRegMOC(TransportProblem &tp)
  : SolverBase(tp)
{
  // Define number of DOFs
  _numDOF = 2 * tp.numEdges * tp.quadOrder * tp.numGroups;
  _numPhaseSpaceDOF = 4 * tp.numCells * tp.quadOrder * tp.numGroups;
  _numSpaceDOF = tp.numCells;
  LOG_DBG("num dof = ",_numDOF);
  LOG_DBG("num space dof = ", tp.numCells);

  mesh = dynamic_cast<MoabMesh*>(_prob.mesh);
  if (!mesh) {
    LOG_ERR("This solver needs a MOAB mesh or a mesh interface with integrated sweep methods.");
    return;
  }

  // Allocate and initialize solution and source arrays
  _solution     = new double [_numDOF];
  _residual     = new double [_numDOF];
  _solutionPrev = new double [_numDOF];
  _source       = new double [_numPhaseSpaceDOF];
  _cellFlux     = new double [_numPhaseSpaceDOF];
  for (long i=0; i<_numDOF; i++) {
    _solution[i] = 0.0;
    _solutionPrev[i] = 0.0;
  }
  for (long i=0; i<_numPhaseSpaceDOF; i++) {
    _source[i] = 0.0;
    _cellFlux[i] = 1.0;
  }

  _mapDOFs();

  _calculateSphericalQuadrature();
}

SolverRegMOC::~SolverRegMOC()
{
  delete [] _solution;
  delete [] _residual;
  delete [] _solutionPrev;
  delete [] _source;
  delete [] _cellFlux;
}

/// Map degrees-of-freedom
void
SolverRegMOC::_mapDOFs()
{
  DOFlist.reserve(_numDOF);

  for (long i=0; i<_prob.numNodes; i++) {
    for (int n=0; n<_prob.quadOrder; n++) {
      for (int g=0; g<_prob.numGroups; g++) {
	for (int c=0; c<2; c++) {
	  DOFlist.push_back( DOFObj(i,n,g,c) );
	  DOFmap.insert( std::pair<DOFObj, double*>(DOFlist[_dofIndex(i,n,g,c)], &_solution[_dofIndex(i,n,g,c)]) );
	}
      }
    }
  }
}

/// Get spherical representation of ordinates
/**
 *  Calculate (theta, mu) coordinates from the input (omega_x, omega_y, omega_z).
 *  Also get the mapping between each ordinate, u, and its negative, -v.  Currently
 *  only two-dimensional problems are assumed, so the polar cosine, mu, is assumed
 *  symmetric about the xy-plane.
 */
void
SolverRegMOC::_calculateSphericalQuadrature()
{
  // Polar coordinate
  _mu = _prob.omega_z;

  // Azimuthal coordinate
  _theta.resize(_prob.quadOrder, 0.0);
  for (int i=0; i<_prob.quadOrder; i++)
    _theta[i] = _getAngleFromVector(_prob.omega_x[i], _prob.omega_y[i]);

  // Aggregate azimuthal negative directions
  _negDir.resize(_prob.quadOrder, -1);
  for (int i=0; i<_prob.quadOrder; i++) {
    bool found;
    found = false;
    int j;
    for (j=0; j<_prob.quadOrder; j++) {
      if (_mu[i] == _mu[j]) {
        if ( std::abs(_theta[i] - fmod(_theta[j] + pi, 2*pi)) < 1e-8) {
          found = true;
          break;
        }
      }
    }
    if (!found)
      LOG_ERR("Could not find negative direction for ordinate ", i);
    else
      _negDir[i] = j;
  }
}

int
SolverRegMOC::_getReflectedDirection(int i, double nx, double ny, double& OmegaDotn)
{
  double Omegapx = -2.0*OmegaDotn*nx + sqrt(1.0-pow(_mu[i],2))*cos(_theta[i]);
  double Omegapy = -2.0*OmegaDotn*ny + sqrt(1.0-pow(_mu[i],2))*sin(_theta[i]);
  double thetap = _getAngleFromVector(Omegapx, Omegapy);

  bool found;
  found = false;
  int j;
  for (j=0; j<_prob.quadOrder; j++) {
    if (_mu[j] == _mu[i]) {
      if ( std::abs(_theta[j] - thetap) < 1e-8) {
        found = true;
        break;
      }
    }
  }
  if (!found) {
    LOG_ERR("Could not find reflected direction for ordinate ", i);
    LOG_DBG(nx," ",ny);
    LOG_DBG(Omegapx," ",Omegapy);
    return -1;
  }
  else
    return j;
}


/// Implementation of the LocalMOC solver
void
SolverRegMOC::solve()
{
  PerfStats X("SolverRegMOC::solve");

  double convRInf = 1e10;
  int innerIter;
  int scatterIter;
  int fissionIter;
  double* fissionSrcFlux = NULL;
  sourceConfig.hasFissionSource = false;

  for (fissionIter=0; fissionIter<_maxIters; fissionIter++) {
    if (sourceConfig.hasFissionSource)
      fissionSrcFlux = copySolution(fissionSrcFlux);
    else
      fissionSrcFlux = _solution;

    // Perform scattering iterations
    for (scatterIter=0; scatterIter<_maxIters; scatterIter++) {
      _saveOldSolution();
      _calculateSource(fissionSrcFlux);
      _applyBoundaryConditions();
      //for (innerIter=0; innerIter<_maxIters; innerIter++) {
      #pragma omp parallel for
      for (int n=0; n<_prob.quadOrder; n++) {
        _sweep(n);
      }
      // Test for convergence of inner iterations
      convRInf = calculateRInfSolutionNorm(_solutionPrev);
      //LOG_DBG("  ",convRInf);
      if (convRInf < _convRInfTol) break;
    }
    // Test for convergence of fission iterations
    printIterStatus("scatter", scatterIter, convRInf, _convRInfTol);
    convRInf = calculateRInfSolutionNorm(fissionSrcFlux);
    if (convRInf < _convRInfTol) break;
  }
  if (sourceConfig.hasFissionSource)
    delete [] fissionSrcFlux;
}

/// Mesh sweep
/**
 *  This function sweeps the mesh in the order in which the elements are listed
 *  in memory--essentially arbitray.
 *
 *  A MoabMesh must be used.
 */
void
SolverRegMOC::_sweep(int n)
{
  UltraLightElement element2;
  long elementID;
  double x[3], y[3], z[3];
  long neighbor[3];
  Sweeper sweep(mesh, n);
  while ((elementID = sweep.getNextElementID()) >= 0) {
    mesh->getCurrentElementFromID(elementID, element2);

    int v0,v1,v2;
    int blockID;
    double theta;
    TriangleDescriptorReg tri;
    _getTriangleOrientation(element2, n, tri, blockID, theta, v0, v1, v2);

    double sigma;
    //double psi0b,psi0a,psi1b,psi1a,psi2a,psi2b;
    double psi0,psi1,psi2,psi3,psi4,psi5;
    long edgeIndex;
    for (int g = 0; g<_prob.numGroups; g++) {
      // EXTERNAL CALL 
      sigma = *mesh->getElementMat(elementID)->getSigma_t(g+1);
      if (blockID==1) {
        if (tri.vertexNeighbor1 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor1);
          psi4 = _solution[ _dofIndex(edgeIndex,n,g,0)];
          psi5 = _solution[ _dofIndex(edgeIndex,n,g,1)];
        }
        else {
          psi4 = 0.0;
          psi5 = 0.0;
          psi4 = _bdryFlux[DOFObj(6*elementID+2,n,g)];
          psi5 = _bdryFlux[DOFObj(6*elementID+3,n,g)];
        }
        if (tri.vertexNeighbor2 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor2);
          psi0 = _solution[ _dofIndex(edgeIndex,n,g,0)];
          psi1 = _solution[ _dofIndex(edgeIndex,n,g,1)];
        }
        else {
          psi0 = 0.0;
          psi1 = 0.0;
          psi0 = _bdryFlux[DOFObj(6*elementID+4,n,g)];
          psi1 = _bdryFlux[DOFObj(6*elementID+5,n,g)];
        }
        tri.sigma = sigma;
        tri.psi0 = psi0;
        tri.psi1 = psi1;
        tri.psi4 = psi4;
        tri.psi5 = psi5;
        triangleSolveA(tri, theta, _theta[n], _mu[n], elementID, n, g);
      }
      else if (blockID==2) {
        // e01 to v2
        if (tri.edgeNeighbor >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.edgeNeighbor);
          psi0 = _solution[ _dofIndex(edgeIndex,n,g,0) ];
          psi1 = _solution[ _dofIndex(edgeIndex,n,g,1) ];
        }
        else {
          psi0 = 0.0;
          psi1 = 0.0;
          psi0 = _bdryFlux[DOFObj(6*elementID,n,g)];
          psi1 = _bdryFlux[DOFObj(6*elementID+1,n,g)];
        }
        tri.sigma = sigma;
        tri.psi0 = psi0;
        tri.psi1 = psi1;
        triangleSolveB(tri, theta, _theta[n], _mu[n], elementID, n, g);
      }
      else if (blockID==3) {
        // v1 to e20
        if (tri.vertexNeighbor1 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor1);
          psi0 = _solution[ _dofIndex(edgeIndex,n,g,0)];
          psi1 = _solution[ _dofIndex(edgeIndex,n,g,1)];
        }
        else {
          psi0 = 0.0;
          psi1 = 0.0;
          psi0 = _bdryFlux[DOFObj(6*elementID+2,n,g)];
          psi1 = _bdryFlux[DOFObj(6*elementID+3,n,g)];
        }
        if (tri.vertexNeighbor2 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor2);
          psi2 = _solution[ _dofIndex(edgeIndex,n,g,0)];
          psi3 = _solution[ _dofIndex(edgeIndex,n,g,1)];
        }
        else {
          psi2 = 0.0;
          psi3 = 0.0;
          psi2 = _bdryFlux[DOFObj(6*elementID+4,n,g)];
          psi3 = _bdryFlux[DOFObj(6*elementID+5,n,g)];
        }
        tri.psi0 = psi2;
        tri.psi1 = psi3;
        tri.psi4 = psi0;
        tri.psi5 = psi1;
        tri.sigma = sigma;
        triangleSolveA(tri, theta, _theta[n], _mu[n], elementID, n, g);
      }
      else if (blockID==4) {
        // e12 to v0
        if (tri.edgeNeighbor >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.edgeNeighbor);
          psi2 = _solution[ _dofIndex(edgeIndex,n,g,0) ];
          psi3 = _solution[ _dofIndex(edgeIndex,n,g,1) ];
        }
        else {
          psi2 = 0.0;
          psi3 = 0.0;
          psi2 = _bdryFlux[DOFObj(6*elementID,n,g)];
          psi3 = _bdryFlux[DOFObj(6*elementID+1,n,g)];
        }
        tri.psi0 = psi2;
        tri.psi1 = psi3;
        tri.sigma = sigma;
        triangleSolveB(tri, theta, _theta[n], _mu[n], elementID, n, g);
      }
      else if (blockID==5) {
        // v2 to e01
        if (tri.vertexNeighbor1 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor1);
          psi2 = _solution[ _dofIndex(edgeIndex,n,g,0)];
          psi3 = _solution[ _dofIndex(edgeIndex,n,g,1)];
        }
        else {
          psi2 = 0.0;
          psi3 = 0.0;
          psi2 = _bdryFlux[DOFObj(6*elementID+2,n,g)];
          psi3 = _bdryFlux[DOFObj(6*elementID+3,n,g)];
        }
        if (tri.vertexNeighbor2 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor2);
          psi4 = _solution[ _dofIndex(edgeIndex,n,g,0)];
          psi5 = _solution[ _dofIndex(edgeIndex,n,g,1)];
        }
        else {
          psi4 = 0.0;
          psi5 = 0.0;
          psi4 = _bdryFlux[DOFObj(6*elementID+4,n,g)];
          psi5 = _bdryFlux[DOFObj(6*elementID+5,n,g)];
        }
        tri.psi0 = psi4;
        tri.psi1 = psi5;
        tri.psi4 = psi2;
        tri.psi5 = psi3;
        tri.sigma = sigma;
        triangleSolveA(tri, theta, _theta[n], _mu[n], elementID, n, g);
      }
      else if (blockID==6) {
        // e20 to v1
        if (tri.edgeNeighbor >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,tri.edgeNeighbor);
          psi4 = _solution[ _dofIndex(edgeIndex,n,g,0) ];
          psi5 = _solution[ _dofIndex(edgeIndex,n,g,1) ];
        }
        else {
          psi4 = 0.0;
          psi5 = 0.0;
          psi4 = _bdryFlux[DOFObj(6*elementID,n,g)];
          psi5 = _bdryFlux[DOFObj(6*elementID+1,n,g)];
        }
        tri.psi0 = psi4;
        tri.psi1 = psi5;
        tri.sigma = sigma;
        triangleSolveB(tri, theta, _theta[n], _mu[n], elementID, n, g);
      }
    } // loop over groups
  } // loop over elements
}

void
SolverRegMOC::_getTriangleOrientation(UltraLightElement &element2, int& n,
                                      TriangleDescriptorReg& tri, int& blockID,
                                      double& theta, int &v0, int &v1, int &v2)
{
  long neighbor[3];
  double x[3],y[3],z[3];

  // Get the local connectivity
  for (int v=0; v<3; v++) {
    x[v] = element2.x[v];
    y[v] = element2.y[v];
    z[v] = element2.z[v];
    neighbor[v] = element2.neighborID[v];
  } // v


  // Main vertex loop
  long neighbor_1[3] = {1, 2, 0};
  long neighbor_2[3] = {2, 0, 1};
  v0 = 0;
    /*
                    e1
          v2 ----------------- v1
             \              /
              \            /
               \          /
             e2 \        / e0
                 \      /
                  \    /
                   \  /
                    \/
                    v0

    */
  double dx, dy;

  // v0-v1 coupling
  v1 = neighbor_1[v0];
  long n0 = neighbor[v0];
  dx = x[v1]-x[v0];
  dy = y[v1]-y[v0];
  double d01 = sqrt( pow(dx, 2) + pow(dy, 2) );
  double phi1 = _getAngleFromVector(dx, dy);

  // v0-v2 coupling
  v2 = neighbor_2[v0];
  
  long n2 = neighbor[v2];
  dx = x[v2]-x[v0];
  dy = y[v2]-y[v0];
  double d20 = sqrt( pow(dx, 2) + pow(dy, 2) );
  double phi2 = _getAngleFromVector(dx,dy);

  // v1-v2 coupling
  long n1 = neighbor[v1];
  dx = x[v2]-x[v1];
  dy = y[v2]-y[v1];
  double d12 = sqrt( pow(dx, 2) + pow(dy, 2) );

  if (phi1 < 1.0e-10 && phi2 > pi)
    phi1 = 2*pi;
  if (phi2 < 1.0e-10 && phi1 > pi)
    phi2 = 2*pi;

  // Sort vertices if needed
  double det = (x[v1]-x[v0])*(y[v2]-y[v0]) - (x[v2]-x[v0])*(y[v1]-y[v0]);
  if ( det < 0.0 ) {
    int vtemp = v1;
    v1 = v2;
    v2 = vtemp;

    long ntemp = n0;
    n0 = n2;
    n2 = ntemp;
    
    double dtemp = d01;
    d01 = d20;
    d20 = dtemp;

    double phitemp = phi1;
    phi1 = phi2;
    phi2 = phitemp;
  } // vertex sort

  double pc[2], pm01[2], pm12[2], pm20[2];
  pc[0] = (x[v0] + x[v1] + x[v2])/3.0;
  pc[1] = (y[v0] + y[v1] + y[v2])/3.0;
  pm01[0] = (x[v0] + x[v1])/2.0;
  pm01[1] = (y[v0] + y[v1])/2.0;
  pm12[0] = (x[v1] + x[v2])/2.0;
  pm12[1] = (y[v1] + y[v2])/2.0;
  pm20[0] = (x[v2] + x[v0])/2.0;
  pm20[1] = (y[v2] + y[v0])/2.0;
  
  double edge01[2], edge12[2], edge20[2];
  double edge01b[2], edge12b[2], edge20b[2];
  double edge0d[2], edge1d[2], edge2d[2];
  edge01[0] = x[v1] - x[v0];    edge01[1] = y[v1] - y[v0];
  edge12[0] = x[v2] - x[v1];    edge12[1] = y[v2] - y[v1];
  edge20[0] = x[v0] - x[v2];    edge20[1] = y[v0] - y[v2];
  edge01b[0] = pc[0] - pm01[0]; edge01b[1] = pc[1] - pm01[1];
  edge12b[0] = pc[0] - pm12[0]; edge12b[1] = pc[1] - pm12[1];
  edge20b[0] = pc[0] - pm20[0]; edge20b[1] = pc[1] - pm20[1];
  edge0d[0] = pc[0] - x[v0];    edge0d[1] = pc[1] - y[v0];
  edge1d[0] = pc[0] - x[v1];    edge1d[1] = pc[1] - y[v1];
  edge2d[0] = pc[0] - x[v2];    edge2d[1] = pc[1] - y[v2];
    
  double theta0,theta1,theta2;
  theta = _theta[n] - phi1;
  theta0 = acos((-edge20[0]*edge01[0]-edge20[1]*edge01[1])/norm(edge20,2)/norm(edge01,2));
  theta1 = acos((-edge01[0]*edge12[0]-edge01[1]*edge12[1])/norm(edge01,2)/norm(edge12,2));
  theta2 = acos((-edge12[0]*edge20[0]-edge12[1]*edge20[1])/norm(edge12,2)/norm(edge20,2));

  if (std::abs(edge01[1]*cos(_theta[n]) - edge01[0]*sin(_theta[n])) < 1.0e-8 ||
      std::abs(edge12[1]*cos(_theta[n]) - edge12[0]*sin(_theta[n])) < 1.0e-8 ||
      std::abs(edge20[1]*cos(_theta[n]) - edge20[0]*sin(_theta[n])) < 1.0e-8) {
    theta += 1.0e-8;
  }

  blockID = 0;
  if (theta < 0.0) theta = theta + 2.0*pi;
  if (std::abs(theta-2.0*pi) < 1.0e-9) theta = 0.0;
  if (fmod(theta,pi) <= theta0) {
    if (theta < pi) {
      blockID = 1;
      // v0 to e12
      tri.edgeNeighbor = n1;
      tri.vertexNeighbor1 = n2;
      tri.vertexNeighbor2 = n0;
      tri.i0 = 0;
      tri.i1 = 1;
      tri.i2 = 2;
      tri.i3 = 3;
      tri.edge01[0] = edge01[0];    tri.edge01[1] = edge01[1];
      tri.edge12[0] = edge12[0];    tri.edge12[1] = edge12[1];
      tri.edge20[0] = edge20[0];    tri.edge20[1] = edge20[1];
      tri.theta0 = theta0; tri.theta1 = theta1; tri.theta2 = theta2;
      tri.d01 = d01;       tri.d12 = d12;       tri.d20 = d20;
    }
    else {
      // e12 to v0
      blockID = 4;
      theta -= (theta0 + theta2);
      tri.edgeNeighbor = n1;
      tri.vertexNeighbor1 = n2;
      tri.vertexNeighbor2 = n0;
      tri.i0 = 3;
      tri.i1 = 1;
      tri.i2 = 0;
      tri.i3 = 2;
      tri.edge01[0] = edge12[0];    tri.edge01[1] = edge12[1];
      tri.edge12[0] = edge20[0];    tri.edge12[1] = edge20[1];
      tri.edge20[0] = edge01[0];    tri.edge20[1] = edge01[1];
      tri.theta0 = theta1; tri.theta1 = theta2; tri.theta2 = theta0;
      tri.d01 = d12;       tri.d12 = d20;       tri.d20 = d01;
    }
  }
  else if (fmod(theta,pi) <= pi-theta1) {
    if (theta <= pi) {
      // e01 to v2
      blockID = 2;
      tri.edgeNeighbor = n0;
      tri.vertexNeighbor1 = n1;
      tri.vertexNeighbor2 = n2;
      tri.i0 = 0;
      tri.i1 = 1;
      tri.i2 = 2;
      tri.i3 = 3;
      tri.edge01[0] = edge01[0];    tri.edge01[1] = edge01[1];
      tri.edge12[0] = edge12[0];    tri.edge12[1] = edge12[1];
      tri.edge20[0] = edge20[0];    tri.edge20[1] = edge20[1];
      tri.theta0 = theta0; tri.theta1 = theta1; tri.theta2 = theta2;
      tri.d01 = d01;       tri.d12 = d12;       tri.d20 = d20;
    }
    else {
      // v2 to e01
      blockID = 5;
      theta -= (pi + theta0);
      tri.edgeNeighbor = n0;
      tri.vertexNeighbor1 = n1;
      tri.vertexNeighbor2 = n2;
      tri.i0 = 2;
      tri.i1 = 1;
      tri.i2 = 3;
      tri.i3 = 0;
      tri.edge01[0] = edge20[0];    tri.edge01[1] = edge20[1];
      tri.edge12[0] = edge01[0];    tri.edge12[1] = edge01[1];
      tri.edge20[0] = edge12[0];    tri.edge20[1] = edge12[1];
      tri.theta0 = theta2; tri.theta1 = theta0; tri.theta2 = theta1;
      tri.d01 = d20;       tri.d12 = d01;       tri.d20 = d12;
    }
  }
  else {
    if (theta <= pi) {
      // v1 to e20
      blockID = 3;
      theta -= (theta0 + theta2);
      tri.edgeNeighbor = n2;
      tri.vertexNeighbor1 = n0;
      tri.vertexNeighbor2 = n1;
      tri.i0 = 3;
      tri.i1 = 1;
      tri.i2 = 0;
      tri.i3 = 2;
      tri.edge01[0] = edge12[0];    tri.edge01[1] = edge12[1];
      tri.edge12[0] = edge20[0];    tri.edge12[1] = edge20[1];
      tri.edge20[0] = edge01[0];    tri.edge20[1] = edge01[1];
      tri.theta0 = theta1; tri.theta1 = theta2; tri.theta2 = theta0;
      tri.d01 = d12;       tri.d12 = d20;       tri.d20 = d01;
    }
    else {
      // e20 to v1
      blockID = 6;
      theta -= (pi + theta0);
      tri.edgeNeighbor = n2;
      tri.vertexNeighbor1 = n0;
      tri.vertexNeighbor2 = n1;
      tri.i0 = 2;
      tri.i1 = 1;
      tri.i2 = 3;
      tri.i3 = 0;
      tri.edge01[0] = edge20[0];    tri.edge01[1] = edge20[1];
      tri.edge12[0] = edge01[0];    tri.edge12[1] = edge01[1];
      tri.edge20[0] = edge12[0];    tri.edge20[1] = edge12[1];
      tri.theta0 = theta2; tri.theta1 = theta0; tri.theta2 = theta1;
      tri.d01 = d20;       tri.d12 = d01;       tri.d20 = d12;
    }
  } // triangle orientation selection
}


void
SolverRegMOC::_applyBoundaryConditions()
{
  UltraLightElement element;
  long edgeIndex;
  int nxtNghbr[3] = {1, 2, 0};
  int* n = new int [4];
  // Loop over all boundary elements
  for (long be=0; be < mesh->boundaryElements.size(); be++) {
    long elementID = mesh->boundaryElements[be];
    
    // Get current element
    mesh->getCurrentElementFromID(elementID, element);

    // Loop over all directions
    for (int n=0; n<_theta.size(); n++) {
      TriangleDescriptorReg tri;
      int blockID;
      int v0,v1,v2;
      double theta;
      _getTriangleOrientation(element, n, tri, blockID, theta,v0,v1,v2);
      double OmegaDotN;

      if (tri.edgeNeighbor < 0) {
        // Edge to vertex case
        // blocks 2,4,6
        double nx,ny;
        switch (blockID) {
        case 2:
        case 5:
          nx = element.y[v1]-element.y[v0];  ny = element.x[v0]-element.x[v1];
          break;
        case 4:
        case 1:
          nx = element.y[v2]-element.y[v1];  ny = element.x[v1]-element.x[v2];
          break;
        case 6:
        case 3:
          nx = element.y[v0]-element.y[v2];  ny = element.x[v2]-element.x[v0];
          break;
        }
        double norm;
        norm = sqrt(pow(nx,2)+pow(ny,2));
        nx = nx/norm;
        ny = ny/norm;
        OmegaDotN = sqrt(1.0-pow(_mu[n],2))*cos(_theta[n])*nx +
          sqrt(1.0-pow(_mu[n],2))*sin(_theta[n])*ny;
        if (OmegaDotN < 0.0) {
          if (_prob.globalBC == reflecting) {
            int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
            for (int g=0; g<_prob.numGroups; g++) {
	      edgeIndex = mesh->getEdgeID(elementID, tri.edgeNeighbor);
              _bdryFlux[DOFObj(6*elementID,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,0) ];
              _bdryFlux[DOFObj(6*elementID+1,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,1) ];
            }
          }
          else if (_prob.globalBC == vacuum) {
            for (int g=0; g<_prob.numGroups; g++) {
              _bdryFlux[DOFObj(6*elementID,n,g)] = 0.0;
              _bdryFlux[DOFObj(6*elementID+1,n,g)] = 0.0;
            }
          }
          else if (_prob.globalBC == source) {
            for (int bi=0; bi<_prob.nxnyq.size(); bi+=3) {
              if (std::abs(nx-_prob.nxnyq[bi]) < 1.0e-8 && std::abs(ny-_prob.nxnyq[bi+1]) < 1.0e-8) {
		if (_prob.nxnyq[bi+2]<0.0) {
		  int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
		  for (int g=0; g<_prob.numGroups; g++) {
		    edgeIndex = mesh->getEdgeID(elementID, tri.edgeNeighbor);
		    _bdryFlux[DOFObj(6*elementID,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,0) ];
		    _bdryFlux[DOFObj(6*elementID+1,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,1) ];
		  }
		}
		else {
		  for (int g=0; g<_prob.numGroups; g++) {
		    _bdryFlux[DOFObj(6*elementID,n,g)] = _prob.nxnyq[bi+2];
		    _bdryFlux[DOFObj(6*elementID+1,n,g)] = _prob.nxnyq[bi+2];
		  }
                }
	      }
            }
          }
        }
      }
      //else {
        // Vertex to edge case
        // blocks 1,3,5
        double nx,ny;

        // Edge 1
        if (tri.vertexNeighbor1 < 0.0) {
          switch (blockID) {
          case 5:
          case 2:
            nx = element.y[v2]-element.y[v1];  ny = element.x[v1]-element.x[v2];
            break;
          case 1:
          case 4:
            nx = element.y[v0]-element.y[v2];  ny = element.x[v2]-element.x[v0];
            break;
          case 3:
          case 6:
            nx = element.y[v1]-element.y[v0];  ny = element.x[v0]-element.x[v1];
            break;
          }
          double norm;
          norm = sqrt(pow(nx,2)+pow(ny,2));
          nx = nx/norm;
          ny = ny/norm;
          OmegaDotN = sqrt(1.0-pow(_mu[n],2))*cos(_theta[n])*nx +
            sqrt(1.0-pow(_mu[n],2))*sin(_theta[n])*ny;
          if (OmegaDotN < 0.0) {
            if (_prob.globalBC == reflecting) {
              int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
              for (int g=0; g<_prob.numGroups; g++) {
		edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor1);
                _bdryFlux[DOFObj(6*elementID+2,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,0) ];
                _bdryFlux[DOFObj(6*elementID+3,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,1) ];
              }
            }
            else if (_prob.globalBC == vacuum) {
              for (int g=0; g<_prob.numGroups; g++) {
                _bdryFlux[DOFObj(6*elementID+2,n,g)] = 0.0;
                _bdryFlux[DOFObj(6*elementID+3,n,g)] = 0.0;
              }
            }
            else if (_prob.globalBC == source) {
              for (int bi=0; bi<_prob.nxnyq.size(); bi+=3) {
                if (std::abs(nx-_prob.nxnyq[bi]) < 1.0e-8 && std::abs(ny-_prob.nxnyq[bi+1]) < 1.0e-8) {
		  if (_prob.nxnyq[bi+2]<0.0) {
		    int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
		    for (int g=0; g<_prob.numGroups; g++) {
		      edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor1);
		      _bdryFlux[DOFObj(6*elementID+2,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,0) ];
		      _bdryFlux[DOFObj(6*elementID+3,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,1) ];
		    }
		  }
		  else {
		    for (int g=0; g<_prob.numGroups; g++) {
		      _bdryFlux[DOFObj(6*elementID+2,n,g)] = _prob.nxnyq[bi+2];
		      _bdryFlux[DOFObj(6*elementID+3,n,g)] = _prob.nxnyq[bi+2];
		    }
		  }
		}
              }
            }
          }
        }

        // Edge 2
        if (tri.vertexNeighbor2 < 0.0) {
          switch (blockID) {
          case 5:
          case 2:
            nx = element.y[v0]-element.y[v2];  ny = element.x[v2]-element.x[v0];
            break;
          case 1:
          case 4:
            nx = element.y[v1]-element.y[v0];  ny = element.x[v0]-element.x[v1];
            break;
          case 3:
          case 6:
            nx = element.y[v2]-element.y[v1];  ny = element.x[v1]-element.x[v2];
            break;
          }
          double norm;
          norm = sqrt(pow(nx,2)+pow(ny,2));
          nx = nx/norm;
          ny = ny/norm;
          OmegaDotN = sqrt(1.0-pow(_mu[n],2))*cos(_theta[n])*nx +
            sqrt(1.0-pow(_mu[n],2))*sin(_theta[n])*ny;
          if (OmegaDotN < 0.0) {
            if (_prob.globalBC == reflecting) {
              int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
              for (int g=0; g<_prob.numGroups; g++) {
		edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor2);
                _bdryFlux[DOFObj(6*elementID+4,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,0) ];
                _bdryFlux[DOFObj(6*elementID+5,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,1) ];
              }
            }
            else if (_prob.globalBC == vacuum) {
              for (int g=0; g<_prob.numGroups; g++) {
                _bdryFlux[DOFObj(6*elementID+4,n,g)] = 0.0;
                _bdryFlux[DOFObj(6*elementID+5,n,g)] = 0.0;
              }
            }
            else if (_prob.globalBC == source) {
              for (int bi=0; bi<_prob.nxnyq.size(); bi+=3) {
                if (std::abs(nx-_prob.nxnyq[bi]) < 1.0e-8 && std::abs(ny-_prob.nxnyq[bi+1]) < 1.0e-8) {
		  if (_prob.nxnyq[bi+2]<0.0) {
		    int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
		    for (int g=0; g<_prob.numGroups; g++) {
		      edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor2);
		      _bdryFlux[DOFObj(6*elementID+4,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,0) ];
		      _bdryFlux[DOFObj(6*elementID+5,n,g)] = _solution[ _dofIndex(edgeIndex,np,g,1) ];
		    }
		  }
		  else {
		    for (int g=0; g<_prob.numGroups; g++) {
		      _bdryFlux[DOFObj(6*elementID+4,n,g)] = _prob.nxnyq[bi+2];
		      _bdryFlux[DOFObj(6*elementID+5,n,g)] = _prob.nxnyq[bi+2];
		    }
		  }
		}
              }
            }
          }
        }
	//}
    }
  }
  delete [] n;
  
}


/// Calculate the azimuthal angle of a vector with components dx and dy
double
SolverRegMOC::_getAngleFromVector(double dx, double dy)
{
  double phi;
  if (std::abs(dx) < 1.0e-8) {
    if (dy > 0.0)
      phi = pi/2;
    else
      phi = 3*pi/2;
  }
  else {
    phi = atan(dy/dx);
    if (dx < 0.0)
      phi += pi;
    if (phi < 0.0)
      phi += 2*pi;
  }
  return phi;
}

void
SolverRegMOC::_calculateSource(double* fissionFlux)
{
  PerfStats X("SolverRegMOC::_calculateSource");

  _getExternalSource();
  //  _addFissionSource(fissionFlux);
  _addScatterSource();
}

void
SolverRegMOC::_getExternalSource()
{
  // Get external source
  if (sourceConfig.hasExternalSource) {
    #pragma omp parallel for
    for (long i=0; i<_prob.numCells; i++) {
      for (int n=0; n<_prob.quadOrder; n++) {
        for (int g=0; g<_prob.numGroups; g++) {
          _source[_dofIndexPS(i,n,g,0)] = *_prob.getExtSource(i,n,g) * sourceScaling;
          _source[_dofIndexPS(i,n,g,1)] = *_prob.getExtSource(i,n,g) * sourceScaling;
          _source[_dofIndexPS(i,n,g,2)] = *_prob.getExtSource(i,n,g) * sourceScaling;
          _source[_dofIndexPS(i,n,g,3)] = *_prob.getExtSource(i,n,g) * sourceScaling;
        }
      }
    }
  }
  else {
    for (long i=0; i<_numPhaseSpaceDOF; i++)
      _source[i] = 0.0;
  }
}

void
SolverRegMOC::_addScatterSource()
{
  // Currently only isotropic scattering is supported
  double scattXS, scalFlux;

  #pragma omp parallel for private(scattXS,scalFlux)
  for (long i=0; i<_prob.numCells; i++) {
    for (int subCell=0; subCell<4; subCell++) {
      for (int gp=0; gp<_prob.numGroups; gp++) {
        scalFlux = getSubCellScalarFlux(i, gp, subCell);
        for (int g=0; g<_prob.numGroups; g++) {
          scattXS = *mesh->getElementMat(i)->getSigma_s(gp+1,g+1);
          for (int n=0; n<_prob.quadOrder; n++) 
            _source[_dofIndexPS(i,n,g,subCell)] += scalFlux*scattXS;
        }
      }
    }
  }
}


void
SolverRegMOC::_calculateMatrixAction(double* x, double* y)
{
}

double
SolverRegMOC::getScalarFlux(long space_i, int group_g)
{
  double scalarFlux = 0;
  for (int subCell=0; subCell<4; subCell++) {
    scalarFlux += getSubCellScalarFlux(space_i, group_g, subCell);
  }
  scalarFlux = scalarFlux/4.0;

  return scalarFlux;
}

double
SolverRegMOC::getSubCellScalarFlux(long space_i, int group_g, int subCell)
{
  double scalarFlux = 0;
  double sumOfWeights = 0;
  for (int n=0; n<_prob.quadOrder; n++) {
    scalarFlux += _prob.weights[n]
      * _cellFlux[_dofIndexPS(space_i,n,group_g,subCell)];
    
    sumOfWeights += _prob.weights[n];
  }

  return scalarFlux/sumOfWeights;
}

void
SolverRegMOC::triangleSolveEV(double& q, double& S, double& sigma, double& x,
                              double& psi1, double& psi2,
                              double& psi01, double& psi20, double& cell, double& psiv)
{
  psi01 = psi1*exp(-sigma*S) + (1-exp(-sigma*S))*q/sigma;
  psi01 = 1.0/(sigma*S)*(psi1 - psi01 + S*q);
  psi01 = psi01 + (psi2-psi1)*(1.0 - exp(-sigma*S) - sigma*S*exp(-sigma*S))/pow(sigma*S,2);

  psi20 = psi2*exp(-sigma*S) + (1-exp(-sigma*S))*q/sigma;;
  psi20 = 1.0/(sigma*S)*(psi2 - psi20 + S*q);
  psi20 = psi20 + (psi1-psi2)*(1.0 - exp(-sigma*S) - sigma*S*exp(-sigma*S))/pow(sigma*S,2);

  psiv = ((psi2-psi1)*(x-0.5) + (psi2+psi1)/2.0)*exp(-sigma*S) + (1-exp(-sigma*S))*q/sigma;

  cell = 1.0/(sigma*S)*((psi2-psi1)/2.0 + psi1 - psi01)
    +  1.0/(sigma*S)*((psi1-psi2)/2.0 + psi2 - psi20) + q/sigma;
}

void
SolverRegMOC::triangleSolveVE2(double& q, double& S, double& sigma,
                               double& w1, double& w2, double& w3,
                               double& psi0r, double& psi0l,
                               double& psi1, double& psi2,
                               double& psi12, double& cell)
{
  double psiv = (w1*psi1 + w2*psi2)/w3;
  
  psi12 = psiv*exp(-sigma*S) + (1-exp(-sigma*S))*q/sigma;
  psi12 = (psiv - psi12 + S*q)/(sigma*S);
  psi12 = psi12 + w2/w3*(psi0l-psi2)*(1.0 - exp(-sigma*S) - sigma*S*exp(-sigma*S))/pow(sigma*S,2);
  psi12 = psi12 + w1/w3*(psi0r-psi1)*(1.0 - exp(-sigma*S) - sigma*S*exp(-sigma*S))/pow(sigma*S,2);

  cell = 2.0/(sigma*S)*(psiv - psi12) + q/sigma;
  cell = cell + 1.0/(sigma*S)*(w2/w3*(psi0l-psi2));
  cell = cell + 1.0/(sigma*S)*(w1/w3*(psi0r-psi1));
}

void
SolverRegMOC::triangleSolveA(TriangleDescriptorReg& tri, double& phi, double& theta, double& mu, long& elementID, int& n, int& g)
{
  double S, psiv, q, x;
  double psi0,psi0l,psi0r,psi1,psi2,w1,w2,w3, delta;
  double deriv;
  long edgeIndex;

  S = tri.d01*sin(tri.theta1)/(2.0*sin(pi-tri.theta1-phi));
  x = S*sin(phi)/sin(tri.theta1)/(tri.d12/2.0);
  S = S/sqrt(1.0 - pow(mu,2));

  // Cell 0: vertex to tri.edge
  q = _source[ _dofIndexPS(elementID,n,g,tri.i0) ];
  w1 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
  w2 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
  w3 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
  psi0r = (3.0*tri.psi0-tri.psi1)/2.0;
  psi0l = (3.0*tri.psi5-tri.psi4)/2.0;
  psi1 = (tri.psi0+tri.psi1)/2.0;
  psi2 = (tri.psi4+tri.psi5)/2.0;
  triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi0r,psi0l,psi1,psi2, tri.psi01, tri.cell0);

  // Cell 1: edge to vertex
  q = _source[ _dofIndexPS(elementID,n,g,tri.i1) ];
  psi1 = (tri.psi4+tri.psi5)/2.0; 
  psi2 = (tri.psi0+tri.psi1)/2.0;
  deriv = psi2-psi1;
  psi2 = tri.psi01 + deriv/2.0;
  psi1 = tri.psi01 - deriv/2.0;
  triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi21,tri.psi13,tri.cell1,psiv);

  //psiv = (w1*tri.psi21 + w2*tri.psi13)/w3;
  
  // Cell 2: vertex to edge
  q = _source[ _dofIndexPS(elementID,n,g,tri.i2) ];
  w1 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
  w2 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
  w3 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
  psi0r = psi1;
  psi1 = psiv;
  deriv = psi0r-psi1;
  psi0r = tri.psi21 + deriv/2.0;
  psi1 = tri.psi21 - deriv/2.0;

  psi0l = (tri.psi4+tri.psi5)/2.0;
  psi2 = (3.0*tri.psi4-tri.psi5)/2.0;
  triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi0r,psi0l,psi1,psi2, tri.psi3, tri.cell2);

  psi1 = (tri.psi4+tri.psi5)/2.0; 
  psi2 = (tri.psi0+tri.psi1)/2.0;
  deriv = psi2-psi1;
  psi2 = tri.psi01 + deriv/2.0;
  psi1 = tri.psi01 - deriv/2.0;

  // Cell 3: vertex to edge
  q = _source[ _dofIndexPS(elementID,n,g,tri.i3) ];
  w1 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
  w2 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
  w3 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
  psi0r = (tri.psi0+tri.psi1)/2.0;
  psi1 = (3.0*tri.psi1-tri.psi0)/2.0;

  psi0l = psi2;
  psi2 = psiv;
  deriv = psi0l-psi2;
  psi0l = tri.psi13 + deriv/2.0;
  psi2 = tri.psi13 - deriv/2.0;
  triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi0r,psi0l,psi1,psi2, tri.psi2, tri.cell3);
  
  if ((tri.psi01 < 0.0 || tri.psi21 < 0.0 || tri.psi13 < 0.0 || tri.psi3 < 0.0 || tri.psi2 < 0.0)) {
    LOG_DBG("NEGATIVE FLUX A-- may not be right-- check edge to vertex");
    q = _source[ _dofIndexPS(elementID,n,g,tri.i0) ];
    w1 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
    w2 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
    w3 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
    psi1 = tri.psi0;
    psi2 = tri.psi5;
    triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi1,psi2,psi1,psi2, tri.psi01, tri.cell0);

    // Cell 1: edge to vertex
    q = _source[ _dofIndexPS(elementID,n,g,tri.i1) ];
    psi1 = (tri.psi4+tri.psi5)/2.0; // wrong... no info from incoming bdry
    psi2 = (tri.psi0+tri.psi1)/2.0;
    delta = tri.psi01 - (psi1+psi2)/2.0;
    psi1 = psi1 + 2.0*delta;
    psi2 = psi2 + 2.0*delta;
    psi1 = tri.psi01;
    psi2 = tri.psi01;
    triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi13,tri.psi21,tri.cell1,psiv);
  
    // Cell 2: vertex to edge
    q = _source[ _dofIndexPS(elementID,n,g,tri.i2) ];
    w1 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
    w2 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
    w3 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
    psi1 = tri.psi21;
    psi2 = tri.psi4;
    triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi1,psi2,psi1,psi2, tri.psi3, tri.cell2);

    // Cell 3: vertex to edge
    q = _source[ _dofIndexPS(elementID,n,g,tri.i3) ];
    w1 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
    w2 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
    w3 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
    psi1 = tri.psi1;
    psi2 = tri.psi13;
    triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi1,psi2,psi1,psi2, tri.psi2, tri.cell3);
  }
  

  edgeIndex = mesh->getEdgeID(elementID,tri.edgeNeighbor);
  _solution[ _dofIndex(edgeIndex, n, g, 0) ] = tri.psi3;
  _solution[ _dofIndex(edgeIndex, n, g, 1) ] = tri.psi2;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i0) ] = tri.cell0;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i1) ] = tri.cell1;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i2) ] = tri.cell2;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i3) ] = tri.cell3;
}

void
SolverRegMOC::triangleSolveB(TriangleDescriptorReg& tri, double& phi, double& theta, double& mu, long& elementID, int& n, int& g)
{
  double S, psiv, q, x;
  double psi0,psi0l,psi0r,psi1,psi2,w1,w2,w3,delta;
  double deriv,psiv0,psiv3;
  long edgeIndex;

  S = tri.d20*sin(tri.theta0)/(2.0*sin(pi-phi));
  x = S*sin(phi-tri.theta0)/sin(tri.theta0)/(tri.d01/2.0);
  S = S/sqrt(1.0 - pow(mu,2));

  // Cell 0: e to v
  q = _source[ _dofIndexPS(elementID,n,g,tri.i0) ];
  psi1 = (3.0*tri.psi0-tri.psi1)/2.0;
  psi2 = (tri.psi0+tri.psi1)/2.0;
  triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi5,tri.psi01,tri.cell0,psiv0);

  // Cell 3: e to v
  q = _source[ _dofIndexPS(elementID,n,g,tri.i3) ];
  psi1 = (tri.psi0+tri.psi1)/2.0;
  psi2 = (3.0*tri.psi1-tri.psi0)/2.0;
  triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi13,tri.psi2,tri.cell3,psiv3);

  // Cell 1: v to e
  q = _source[ _dofIndexPS(elementID,n,g,tri.i1) ];
  w1 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
  w2 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
  w3 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
  psi0r = (tri.psi0+tri.psi1)/2.0;//(w1*tri.psi13+w2*tri.psi01)/w3; the commented is bad
  psi1 = psiv3;//(w1*tri.psi13+w2*tri.psi2)/w3;
  deriv = psi0r-psi1;
  psi0r = tri.psi13 + deriv/2.0;
  psi1 = tri.psi13 - deriv/2.0;
  
  psi0l = (tri.psi0+tri.psi1)/2.0;//psi0r; the commented is bad
  psi2 = psiv0;//(w1*tri.psi5+w2*tri.psi01)/w3;
  deriv = psi0l-psi2;
  psi0l = tri.psi01 + deriv/2.0;
  psi2 = tri.psi01 - deriv/2.0;
  triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi0r,psi0l,psi1,psi2, tri.psi21, tri.cell1);
  
  // Cell 2: e to v
  q = _source[ _dofIndexPS(elementID,n,g,tri.i2) ];
  psi1 = (w1*tri.psi5+w2*tri.psi01)/w3;
  psi2 = (w1*tri.psi13+w2*tri.psi2)/w3;  // may need to switch
  psi1 = psiv0;
  psi2 = psiv3;
  deriv = psi2-psi1;
  psi2 = tri.psi21 + deriv/2.0;
  psi1 = tri.psi21 - deriv/2.0;
  triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi4,tri.psi3,tri.cell2,psiv);

  
  if ((tri.psi01<0.0 || tri.psi21<0.0 || tri.psi13<0.0 || tri.psi3<0.0 || tri.psi2<0.0 || tri.psi4<0.0 || tri.psi5<0.0)) {
    LOG_DBG("Going to zero order");
    LOG_DBG(elementID, " ", n, " ", g);
    // Cell 0: e to v
    q = _source[ _dofIndexPS(elementID,n,g,tri.i0) ];
    psi1 = tri.psi0;
    psi2 = tri.psi0;
    triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi01,tri.psi5,tri.cell0,psiv);

    // Cell 3: e to v
    q = _source[ _dofIndexPS(elementID,n,g,tri.i3) ];
    psi1 = tri.psi1;
    psi2 = tri.psi1;
    triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi2,tri.psi13,tri.cell3,psiv);

    // Cell 1: v to e
    q = _source[ _dofIndexPS(elementID,n,g,tri.i1) ];
    w1 = std::abs(-tri.edge20[1]*cos(theta)+tri.edge20[0]*sin(theta))/2;
    w2 = std::abs(-tri.edge12[1]*cos(theta)+tri.edge12[0]*sin(theta))/2;
    w3 = std::abs(-tri.edge01[1]*cos(theta)+tri.edge01[0]*sin(theta))/2;
    psi1 = tri.psi13;
    psi2 = tri.psi01;
    triangleSolveVE2(q,S,tri.sigma,w1,w2,w3,psi1,psi2,psi1,psi2, tri.psi21, tri.cell1);

    // Cell 2: e to v
    q = _source[ _dofIndexPS(elementID,n,g,tri.i2) ];
    psi1 = tri.psi21; // could impose a shape from the boundary points
    psi2 = tri.psi21;
    triangleSolveEV(q,S,tri.sigma,x,psi1,psi2,tri.psi3,tri.psi4,tri.cell2,psiv);
  }
  

  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor1);
  _solution[ _dofIndex(edgeIndex, n, g, 0) ] = tri.psi3;
  _solution[ _dofIndex(edgeIndex, n, g, 1) ] = tri.psi2;
  edgeIndex = mesh->getEdgeID(elementID,tri.vertexNeighbor2);
  _solution[ _dofIndex(edgeIndex, n, g, 0) ] = tri.psi5;
  _solution[ _dofIndex(edgeIndex, n, g, 1) ] = tri.psi4;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i0) ] = tri.cell0;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i1) ] = tri.cell1;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i2) ] = tri.cell2;
  _cellFlux[ _dofIndexPS(elementID,n,g,tri.i3) ] = tri.cell3;
}
