#include "solverdummy.h"
#include "global.h"
#include "sweeper.h"

SolverDummy::SolverDummy(TransportProblem &tp)
  : SolverBase(tp)
{
  // Define number of DOFs
  _numDOF = tp.numCells * tp.quadOrder * tp.numGroups;
  _numSpaceDOF = tp.numCells;
  LOG_DBG("num dof = ",_numDOF);

  mesh = dynamic_cast<MoabMesh*>(_prob.mesh);
  if (!mesh) {
    LOG_ERR("This solver needs a MOAB mesh or a mesh interface with integrated sweep methods.");
    return;
  }

  // Allocate and initialize solution and source arrays
  _mapDOFs();
  _calculateSphericalQuadrature();
}

SolverDummy::~SolverDummy()
{
}

/// Map degrees-of-freedom
void
SolverDummy::_mapDOFs()
{
  DOFlist.reserve(_numDOF);

  for (long i=0; i<_prob.numCells; i++) {
    for (int n=0; n<_prob.quadOrder; n++) {
      for (int g=0; g<_prob.numGroups; g++) {
        DOFlist.push_back( DOFObj(i,n,g) );
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
SolverDummy::_calculateSphericalQuadrature()
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

/// Implementation of the LocalMOC solver
void
SolverDummy::solve()
{
  PerfStats X("SolverDummy::solve");

  #pragma omp parallel for
  for (int n=0; n<_prob.quadOrder; n++)
    _sweep(n);

}

/// Mesh sweep
/**
 *  A MoabMesh must be used.
 */
void
SolverDummy::_sweep(int n)
{
  UltraLightElement element2;
  long elementID;
  double x[3], y[3], z[3];
  long id[3];
  long neighbor[3];
  Sweeper sweep(mesh, n);
  while ((elementID = sweep.getNextElementID()) >=0 ) {
    sweep.getCurrentElementFromID(elementID, element2);

    // Get the local connectivity
    for (int v=0; v<3; v++) {
      x[v] = element2.x[v];
      y[v] = element2.y[v];
      z[v] = element2.z[v];
      id[v] = element2.vertexID[v];
      neighbor[v] = element2.neighborID[v];
    } // v



    // Main vertex loop
    long neighbor_1[3] = {1, 2, 0};
    long neighbor_2[3] = {2, 0, 1};
    int v0 = 0;
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
    int v1 = neighbor_1[v0];
    long n0 = neighbor[v0];
    dx = x[v1]-x[v0];
    dy = y[v1]-y[v0];
    double d01 = sqrt( pow(dx, 2) + pow(dy, 2) );
    double phi1 = _getAngleFromVector(dx, dy);

    // v0-v2 coupling
    int v2 = neighbor_2[v0];
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

    double theta0 = acos(((x[v2]-x[v0])*(x[v1]-x[v0]) + (y[v2]-y[v0])*(y[v1]-y[v0]))/(d01*d20));
    double theta1 = acos(((x[v0]-x[v1])*(x[v2]-x[v1]) + (y[v0]-y[v1])*(y[v2]-y[v1]))/(d12*d01));
    double theta2 = acos(((x[v1]-x[v2])*(x[v0]-x[v2]) + (y[v1]-y[v2])*(y[v0]-y[v2]))/(d20*d12));
      
    // Integrate each characteristic
    double theta,mu01,mu12,mu20,pathDist;
    int evDir,veDir;
    long edgeNeighbor,vertexNeighbor1,vertexNeighbor2;
    
    theta = _theta[n] - phi1;
    if (theta < 0.0) theta = theta + 2.0*pi;

    if (fmod(theta,pi) <= theta0) {
      // v0 to e12
      edgeNeighbor = n1;
      vertexNeighbor1 = n2;
      vertexNeighbor2 = n0;
      pathDist = d01*sin(theta1)/std::abs(sin(theta1+theta));
      mu20 = d20*sin(theta0-theta);
      mu01 = d01*sin(theta);
      mu12 = d12*sin(theta+theta1);
      if (theta < pi) {
        evDir = _negDir[n];
        veDir = n;
      }
      else {
        veDir = _negDir[n];
        evDir = n;
      }
    }
    else if (fmod(theta,pi) <= pi-theta1) {
      // e01 to v2
      theta = theta - theta0;
      edgeNeighbor = n0;
      vertexNeighbor1 = n1;
      vertexNeighbor2 = n2;
      pathDist = d20*sin(theta0)/std::abs(sin(theta0+theta));
      mu20 = d12*sin(theta2-theta);
      mu01 = d20*sin(theta);
      mu12 = d01*sin(theta+theta0);
      if (theta <= pi) {
        veDir = _negDir[n];
        evDir = n;
      }
      else {
        evDir = _negDir[n];
        veDir = n;
      }
    }
    else {
      // v1 to e20
      theta = theta - (pi - theta1);
      edgeNeighbor = n2;
      vertexNeighbor1 = n0;
      vertexNeighbor2 = n1;
      pathDist = d12*sin(theta2)/std::abs(sin(theta2+theta));
      mu20 = d01*sin(theta1-theta);
      mu01 = d12*sin(theta);
      mu12 = d20*sin(theta+theta2);
      if (theta <= pi) {
        evDir = _negDir[n];
        veDir = n;
      }
      else {
        veDir = _negDir[n];
        evDir = n;
      }
    } // triangle orientation selection

    if (evDir==n) {
    }
    else {
    }
    
  } // loop over elements
}


/// Calculate the azimuthal angle of a vector with components dx and dy
double
SolverDummy::_getAngleFromVector(double dx, double dy)
{
  double phi;
  if (std::abs(dx) < 1.0e-8) {
    if (dy > 0.0)
      phi = pi/2;
    else
      phi = 3*pi/2;
  }
  else {
    phi = fmod(2.0*pi + atan2(dy,dx), 2.0*pi);
  }
  return phi;
}


void
SolverDummy::_calculateMatrixAction(double* x, double* y)
{
}

double
SolverDummy::getScalarFlux(long space_i, int group_g)
{
  return 0.0;
}
