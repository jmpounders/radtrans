#include "solverlocalmoc.h"
#include "global.h"
#include "sweeper.h"

#include <iomanip>

SolverLocalMOC::SolverLocalMOC(TransportProblem &tp)
  : SolverBase(tp)
{
  // Define number of DOFs
  _numDOF = 2 * tp.numEdges * tp.quadOrder * tp.numGroups;
  _numPhaseSpaceDOF = 2 * tp.numCells * tp.quadOrder * tp.numGroups;
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
  _surfacePosition = new double [_numDOF];
  _source       = new double [_numPhaseSpaceDOF];
  _cellFlux     = new double [_numPhaseSpaceDOF];
  for (long i=0; i<_numDOF; i++) {
    _solution[i] = 0.0;
    _solutionPrev[i] = 0.0;
    _surfacePosition[i] = 0.0;
  }
  for (long i=0; i<_numPhaseSpaceDOF; i++) {
    _source[i] = 0.0;
    _cellFlux[i] = 1.0;
  }

  _mapDOFs();

  _calculateSphericalQuadrature();
}

SolverLocalMOC::~SolverLocalMOC()
{
  delete [] _solution;
  delete [] _residual;
  delete [] _solutionPrev;
  //delete [] _surfacePosition;
  delete [] _source;
  delete [] _cellFlux;
}

/// Map degrees-of-freedom
void
SolverLocalMOC::_mapDOFs()
{
  DOFlist.reserve(_numDOF);

  for (long i=0; i<_prob.numCells; i++) {
    for (int n=0; n<_prob.quadOrder; n++) {
      for (int g=0; g<_prob.numGroups; g++) {
        DOFlist.push_back( DOFObj(i,n,g) );
        DOFmap.insert( std::pair<DOFObj, double*>(DOFlist[_dofIndex(i,n,g)], &_solution[_dofIndex(i,n,g)]) );
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
SolverLocalMOC::_calculateSphericalQuadrature()
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
SolverLocalMOC::_getReflectedDirection(int i, double nx, double ny, double& OmegaDotn)
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
SolverLocalMOC::solve()
{
  PerfStats X("SolverLocalMOC::solve");

  double convRInf;
  int innerIter;
  int scatterIter;
  int fissionIter;
  double* fissionSrcFlux = NULL;

  sourceConfig.hasFissionSource = false;
  
  convRInf = 1.0e10;

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
      #pragma omp parallel for
      for (int n=0; n<_prob.quadOrder; n++)
        _sweep(n);
      // Test for convergence of inner iterations
      convRInf = calculateRInfSolutionNorm(_solutionPrev);
      if (convRInf < _convRInfTol) break;
    }
    // Test for convergence of fission iterations
    printIterStatus("scatter", scatterIter, convRInf, _convRInfTol);
    convRInf = calculateRInfSolutionNorm(fissionSrcFlux);
    if (convRInf < _convRInfTol) break;
  }
  //printIterStatus("fission", fissionIter, convRInf, _convRInfTol);
  if (sourceConfig.hasFissionSource)
    delete [] fissionSrcFlux;
}

/// Mesh sweep
/**
 *  A MoabMesh must be used.
 */
void
SolverLocalMOC::_sweep(int n)
{
  UltraLightElement element;
  long elementID;
  long id[3];
  long neighbor[3];
  Sweeper sweep(mesh, n);
  while ((elementID = sweep.getNextElementID()) >=0 ) {
    mesh->getCurrentElementFromID(elementID, element);

    double mu01,mu12,mu20;
    double pathDist;
    int evDir,veDir;
    long edgeNeighbor,vertexNeighbor1,vertexNeighbor2;
    int xedge, v0,v1,v2;
    double surfacePosition;
    _getTriangleOrientation(element, n,
                            mu01, mu12, mu20,
                            pathDist, evDir, veDir,
                            edgeNeighbor, vertexNeighbor1, vertexNeighbor2,
                            xedge,v0,v1,v2,surfacePosition);
          
    // Do calculation here
    for (int g=0; g<_prob.numGroups; g++) {
      // Do edge to vertex characteristic
      double psi0,psi1,psi2,psi12,psi01,psi20, q, att,expatt, sigma;
      long edgeIndex;
      sigma = *mesh->getElementMat(elementID)->getSigma_t(g+1);
      att = pathDist/sqrt(1.0 - pow(_mu[n],2));
      expatt = exp(-sigma*att);

      if (evDir==n) {
        // Do edge to vertex characteristic
	double sp,psi12l,psi12r;
        if (edgeNeighbor >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID, edgeNeighbor);
	  sp = _surfacePosition[ _dofIndex(edgeIndex,evDir,g) ];
	  psi12l = _solution[ _dofIndex(edgeIndex,evDir,g,0) ];
	  psi12r = _solution[ _dofIndex(edgeIndex,evDir,g,1) ];
	  if (surfacePosition <= sp) {
	    psi12r = ((sp-surfacePosition)*psi12l + (1-sp)*psi12r)/(1-surfacePosition);
	    psi12l = psi12l;
	  }
	  else {
	    psi12l = (sp*psi12l + (surfacePosition-sp)*psi12r)/surfacePosition;
	    psi12r = psi12r;
	  }
	}
        else {
          psi12l = _bdryFlux[DOFObj(3*elementID,n,g)];
          psi12r = _bdryFlux[DOFObj(3*elementID,n,g)];
        }
        
        q = _source[ _dofIndex(elementID,evDir,g) ];

	edgeIndex = mesh->getEdgeID(vertexNeighbor1,elementID);
        psi0 = expatt*psi12r + (1.0 - expatt)/sigma*q;
        _solution[ _dofIndex(edgeIndex,evDir,g,0) ] = (psi12r - psi0)/(sigma*att) + q/sigma;
        _solution[ _dofIndex(edgeIndex,evDir,g,1) ] = (psi12r - psi0)/(sigma*att) + q/sigma;
	_surfacePosition[ _dofIndex(edgeIndex,evDir,g) ] = 0.5;

	edgeIndex = mesh->getEdgeID(vertexNeighbor2,elementID);
        psi0 = expatt*psi12l + (1.0 - expatt)/sigma*q;
        _solution[ _dofIndex(edgeIndex,evDir,g,0) ] = (psi12l - psi0)/(sigma*att) + q/sigma;
        _solution[ _dofIndex(edgeIndex,evDir,g,1) ] = (psi12l - psi0)/(sigma*att) + q/sigma;
	_surfacePosition[ _dofIndex(edgeIndex,evDir,g) ] = 0.5;

	psi12 = surfacePosition*psi12l + (1-surfacePosition)*psi12r;
        psi0 = expatt*psi12 + (1.0 - expatt)/sigma*q;
        _cellFlux[ _dofIndex(elementID,evDir,g) ] = (psi12 - ((psi12 - psi0)/(sigma*att) + q/sigma))*2/(sigma*att) + q/sigma;
      }
      else {
        // Do vertex to edge characteristic
	double sp;
        if (vertexNeighbor1 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,vertexNeighbor1);
	  sp = _surfacePosition[ _dofIndex(edgeIndex,veDir,g) ];
          psi20 = _solution[ _dofIndex(edgeIndex,veDir,g,0) ]*sp;
	  psi20 += _solution[ _dofIndex(edgeIndex,veDir,g,1) ]*(1.0-sp);
	}
        else {
          psi20 = _bdryFlux[DOFObj(3*elementID+1,n,g)];
        }
        if (vertexNeighbor2 >= 0) {
	  edgeIndex = mesh->getEdgeID(elementID,vertexNeighbor2);
	  sp = _surfacePosition[ _dofIndex(edgeIndex,veDir,g) ];
          psi01 = _solution[ _dofIndex(edgeIndex,veDir,g,0) ]*sp;
          psi01 += _solution[ _dofIndex(edgeIndex,veDir,g,1) ]*(1.0-sp);
	}
        else {
          psi01 = _bdryFlux[DOFObj(3*elementID+2,n,g)];
        }
        
        q = _source[ _dofIndex(elementID,veDir,g) ];

	edgeIndex = mesh->getEdgeID(edgeNeighbor,elementID);
	// Left side
        psi12 = expatt*psi20 + (1.0 - expatt)/sigma*q;
        _solution[ _dofIndex(edgeIndex,veDir,g,0) ] = (psi20 - psi12)/(sigma*att) + q/sigma;
	// Right side
        psi12 = expatt*psi01 + (1.0 - expatt)/sigma*q;
        _solution[ _dofIndex(edgeIndex,veDir,g,1) ] = (psi01 - psi12)/(sigma*att) + q/sigma;
	_surfacePosition[ _dofIndex(edgeIndex,veDir,g) ] = surfacePosition;

        psi0 = (mu20*psi20 + mu01*psi01) / mu12;
        psi12 = expatt*psi0 + (1.0 - expatt)/sigma*q;
        _cellFlux[ _dofIndex(elementID,veDir,g) ] = (psi0 - ((psi0 - psi12)/(sigma*att) + q/sigma))*2/(sigma*att) + q/sigma;
      }
    } // loop over groups
  } // loop over elements

}

void
SolverLocalMOC::_getTriangleOrientation(UltraLightElement &element, int n,
                                        double &mu01, double &mu12, double &mu20,
                                        double &pathDist, int &evDir, int &veDir,
                                        long &edgeNeighbor, long &vertexNeighbor1, long &vertexNeighbor2,
                                        int &xedge, int &v0, int &v1, int &v2, double &surfacePosition)
{
  //long id[3];
  long neighbor[3];
  double x[3],y[3],z[3];

  // Get the local connectivity
  for (int v=0; v<3; v++) {
    x[v] = element.x[v];
    y[v] = element.y[v];
    z[v] = element.z[v];
    //id[v] = element.vertexID[v];
    neighbor[v] = element.neighborID[v];
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

  // Sort vertices if needed to ensure counter-clockwise orientation
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

  // Note: I had to use the min function in the definition of theat1
  // because the argument would sometimes slip maringally above one
  //double theta0 = phi2 - phi1;
  //double theta1 = asin(std::min(1.0,d20*sin(theta0)/d12));
  //double theta2 = pi - theta0 - theta1;
  double theta0 = acos(((x[v2]-x[v0])*(x[v1]-x[v0]) + (y[v2]-y[v0])*(y[v1]-y[v0]))/(d01*d20));
  double theta1 = acos(((x[v0]-x[v1])*(x[v2]-x[v1]) + (y[v0]-y[v1])*(y[v2]-y[v1]))/(d12*d01));
  double theta2 = acos(((x[v1]-x[v2])*(x[v0]-x[v2]) + (y[v1]-y[v2])*(y[v0]-y[v2]))/(d20*d12));

  // Integrate each characteristic
  double theta;
    
  theta = _theta[n] - phi1;
  if (theta < 0.0) theta = theta + 2.0*pi;
  if (std::abs(theta-2.0*pi) < 1.0e-9) theta = 0.0;

  double edge01[2], edge12[2], edge20[2];
  edge01[0] = x[v1] - x[v0];    edge01[1] = y[v1] - y[v0];
  edge12[0] = x[v2] - x[v1];    edge12[1] = y[v2] - y[v1];
  edge20[0] = x[v0] - x[v2];    edge20[1] = y[v0] - y[v2];
  if (std::abs(edge01[1]*cos(_theta[n]) - edge01[0]*sin(_theta[n])) < 1.0e-8 ||
      std::abs(edge12[1]*cos(_theta[n]) - edge12[0]*sin(_theta[n])) < 1.0e-8 ||
      std::abs(edge20[1]*cos(_theta[n]) - edge20[0]*sin(_theta[n])) < 1.0e-8) {
    theta += 1.0e-8;
  }


  if (fmod(theta,pi) <= theta0) {
    // v0 to e12
    xedge = 1;
    edgeNeighbor = n1;
    vertexNeighbor1 = n2;
    vertexNeighbor2 = n0;
    pathDist = d01*sin(theta1)/std::abs(sin(theta1+theta));
    surfacePosition = pathDist*std::abs(sin(theta0-theta)/sin(theta2))/d20;
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
      surfacePosition = 1.0 - surfacePosition;
    }
  }
  else if (fmod(theta,pi) <= pi-theta1) {
    // e01 to v2
    // why do i have to use abs?  are the theta values inconsistent?
    theta = theta - theta0;
    xedge = 0;
    edgeNeighbor = n0;
    vertexNeighbor1 = n1;
    vertexNeighbor2 = n2;
    pathDist = d20*sin(theta0)/std::abs(sin(theta0+theta));
    surfacePosition = pathDist*std::abs(sin(theta1-theta)/sin(theta0))/d01;
    mu20 = d12*sin(theta2-theta);
    mu01 = d20*sin(theta);
    mu12 = d01*sin(theta+theta0);
    if (theta <= pi) {
      veDir = _negDir[n];
      evDir = n;
      surfacePosition = 1.0 - surfacePosition;
    }
    else {
      evDir = _negDir[n];
      veDir = n;
    }
  }
  else {
    // v1 to e20
    theta = theta - (pi - theta1);
    xedge = 2;
    edgeNeighbor = n2;
    vertexNeighbor1 = n0;
    vertexNeighbor2 = n1;
    pathDist = d12*sin(theta2)/std::abs(sin(theta2+theta));
    surfacePosition = pathDist*std::abs(sin(theta2-theta)/sin(theta1))/d12;
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
      surfacePosition = 1.0 - surfacePosition;
    }
  } // triangle orientation selection
}

void
SolverLocalMOC::_applyBoundaryConditions()
{
  UltraLightElement element;
  int nxtNghbr[3] = {1, 2, 0};
  int* n = new int [4];

  // Loop over all boundary elements
  for (long be=0; be < mesh->boundaryElements.size(); be++) {
    long elementID = mesh->boundaryElements[be];
    
    // Get current element
    mesh->getCurrentElementFromID(elementID, element);

    // Loop over all directions
    for (int n=0; n<_theta.size(); n++) {
      double mu01,mu12,mu20;
      double pathDist;
      int evDir,veDir;
      long edgeNeighbor,vertexNeighbor1,vertexNeighbor2;
      int xedge, v0,v1,v2;
      double surfacePosition;
      _getTriangleOrientation(element, n,
                              mu01, mu12, mu20,
                              pathDist, evDir, veDir,
                              edgeNeighbor, vertexNeighbor1, vertexNeighbor2,
                              xedge,v0,v1,v2,surfacePosition);
      double OmegaDotN;
      
      if (edgeNeighbor < 0 && evDir==n) {
        // Edge to vertex case
        double nx,ny;
        switch (xedge) {
        case 0:
          nx = element.y[v1]-element.y[v0];  ny = element.x[v0]-element.x[v1];
          break;
        case 1:
          nx = element.y[v2]-element.y[v1];  ny = element.x[v1]-element.x[v2];
          break;
        case 2:
          nx = element.y[v0]-element.y[v2];  ny = element.x[v2]-element.x[v0];
          break;
        }
        double norm,OmegaDotN;
        norm = sqrt(pow(nx,2)+pow(ny,2));
        nx = nx/norm;
        ny = ny/norm;
        OmegaDotN = sqrt(1.0-pow(_mu[n],2))*cos(_theta[n])*nx +
          sqrt(1.0-pow(_mu[n],2))*sin(_theta[n])*ny;
        if (OmegaDotN < 0.0) {
          if (_prob.globalBC == reflecting) {
            int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
            for (int g=0; g<_prob.numGroups; g++) 
              _bdryFlux[DOFObj(3*elementID,n,g)] = _solution[ _dofIndex(elementID,np,g) ];
          }
          else if (_prob.globalBC == vacuum) {
            for (int g=0; g<_prob.numGroups; g++) 
              _bdryFlux[DOFObj(3*elementID,n,g)] = 0.0;
          }
          else if (_prob.globalBC == source) {
            for (int bi=0; bi<_prob.nxnyq.size(); bi+=3) {
              if (std::abs(nx-_prob.nxnyq[bi]) < 1.0e-8 && std::abs(ny-_prob.nxnyq[bi+1]) < 1.0e-8) {
		if (_prob.nxnyq[bi+2]<0.0) {
		  int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
		  for (int g=0; g<_prob.numGroups; g++) 
		    _bdryFlux[DOFObj(3*elementID,n,g)] = _solution[ _dofIndex(elementID,np,g) ];
		}
		else {
		  for (int g=0; g<_prob.numGroups; g++)
		    _bdryFlux[DOFObj(3*elementID,n,g)] = _prob.nxnyq[bi+2];
		}
	      }
            }
          }
        }
      }
      else {
        // Vertex to edge case
        double nx,ny;

        // Edge 1
        if (vertexNeighbor1 < 0.0) {
          switch (xedge) {
          case 0:
            nx = element.y[v2]-element.y[v1];  ny = element.x[v1]-element.x[v2];
            break;
          case 1:
            nx = element.y[v0]-element.y[v2];  ny = element.x[v2]-element.x[v0];
            break;
          case 2:
            nx = element.y[v1]-element.y[v0];  ny = element.x[v0]-element.x[v1];
            break;
          }
          double norm = sqrt(pow(nx,2)+pow(ny,2));
          nx = nx/norm;
          ny = ny/norm;
          OmegaDotN = sqrt(1.0-pow(_mu[n],2))*cos(_theta[n])*nx +
            sqrt(1.0-pow(_mu[n],2))*sin(_theta[n])*ny;
          if (OmegaDotN < 0.0) {
            if (_prob.globalBC == reflecting) {
              int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
              for (int g=0; g<_prob.numGroups; g++) 
                _bdryFlux[DOFObj(3*elementID+1,n,g)] = _solution[ _dofIndex(elementID,np,g) ];
            }
            else if (_prob.globalBC == vacuum) {
              for (int g=0; g<_prob.numGroups; g++) 
                _bdryFlux[DOFObj(3*elementID+1,n,g)] = 0.0;
            }
            else if (_prob.globalBC == source) {
              for (int bi=0; bi<_prob.nxnyq.size(); bi+=3) {
                if (std::abs(nx-_prob.nxnyq[bi]) < 1.0e-8 && std::abs(ny-_prob.nxnyq[bi+1]) < 1.0e-8) {
		  if (_prob.nxnyq[bi+2]<0.0) {
		    int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
		    for (int g=0; g<_prob.numGroups; g++) 
		      _bdryFlux[DOFObj(3*elementID+1,n,g)] = _solution[ _dofIndex(elementID,np,g) ];
		  }
		  else{
		    for (int g=0; g<_prob.numGroups; g++)
		      _bdryFlux[DOFObj(3*elementID+1,n,g)] = _prob.nxnyq[bi+2];
		  }
		}
              }
            }
          }
        }

        // Edge 2
        if (vertexNeighbor2 < 0.0) {
          switch (xedge) {
          case 0:
            nx = element.y[v0]-element.y[v2];  ny = element.x[v2]-element.x[v0];
            break;
          case 1:
            nx = element.y[v1]-element.y[v0];  ny = element.x[v0]-element.x[v1];
            break;
          case 2:
            nx = element.y[v2]-element.y[v1];  ny = element.x[v1]-element.x[v2];
            break;
          }
          double norm = sqrt(pow(nx,2)+pow(ny,2));
          nx = nx/norm;
          ny = ny/norm;
          OmegaDotN = sqrt(1.0-pow(_mu[n],2))*cos(_theta[n])*nx +
            sqrt(1.0-pow(_mu[n],2))*sin(_theta[n])*ny;
          if (OmegaDotN < 0.0) {
            if (_prob.globalBC == reflecting) {
              int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
              for (int g=0; g<_prob.numGroups; g++)
                _bdryFlux[DOFObj(3*elementID+2,n,g)] = _solution[ _dofIndex(elementID,np,g) ];
            }
            else if (_prob.globalBC == vacuum) {
              for (int g=0; g<_prob.numGroups; g++) 
                _bdryFlux[DOFObj(3*elementID+2,n,g)] = 0.0;
            }
            else if (_prob.globalBC == source) {
              for (int bi=0; bi<_prob.nxnyq.size(); bi+=3) {
                if (std::abs(nx-_prob.nxnyq[bi]) < 1.0e-8 && std::abs(ny-_prob.nxnyq[bi+1]) < 1.0e-8) {
		  if (_prob.nxnyq[bi+2]<0.0) {
		    int np = _getReflectedDirection(n, nx, ny, OmegaDotN);
		    for (int g=0; g<_prob.numGroups; g++)
		      _bdryFlux[DOFObj(3*elementID+2,n,g)] = _solution[ _dofIndex(elementID,np,g) ];
		  }
		  else {
		    for (int g=0; g<_prob.numGroups; g++)
		      _bdryFlux[DOFObj(3*elementID+2,n,g)] = _prob.nxnyq[bi+2];
		  }
		}
              }
            }
          }
        }
      }
    }
  }
  delete [] n;
}


/// Calculate the azimuthal angle of a vector with components dx and dy
double
SolverLocalMOC::_getAngleFromVector(double dx, double dy)
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
SolverLocalMOC::_calculateSource(double* fissionFlux)
{
  PerfStats X("SolverLocalMOC::_calculateSource");

  _getExternalSource();
  //  _addFissionSource(fissionFlux);
  _addScatterSource();
}

void
SolverLocalMOC::_getExternalSource()
{
  // Get external source
  if (sourceConfig.hasExternalSource) {
    #pragma omp parallel for
    for (long i=0; i<_prob.numCells; i++) {
      for (int n=0; n<_prob.quadOrder; n++) {
        for (int g=0; g<_prob.numGroups; g++) {
          _source[_dofIndex(i,n,g)] = *_prob.getExtSource(i,n,g) * sourceScaling;
        }
      }
    }
  }
  else {
    for (long i=0; i<_numDOF; i++)
      _source[i] = 0.0;
  }
}

void
SolverLocalMOC::_addScatterSource()
{
  // Currently only isotropic scattering is supported
  double scattXS, scalFlux;

  #pragma omp parallel for private(scattXS,scalFlux)
  for (long i=0; i<_prob.numCells; i++) {
    for (int gp=0; gp<_prob.numGroups; gp++) {
      scalFlux = getScalarFlux(i, gp);
      for (int g=0; g<_prob.numGroups; g++) {
        scattXS = *mesh->getElementMat(i)->getSigma_s(gp+1,g+1);
        for (int n=0; n<_prob.quadOrder; n++)
          _source[_dofIndex(i,n,g)] += scalFlux*scattXS;
      }
    }
  }
}


void
SolverLocalMOC::_calculateMatrixAction(double* x, double* y)
{
}

double
SolverLocalMOC::getScalarFlux(long space_i, int group_g)
{
  double scalarFlux = 0;
  double sumOfWeights = 0;
  for (int n=0; n<_prob.quadOrder; n++) {
    scalarFlux += _prob.weights[n] * _cellFlux[ _dofIndex(space_i, n, group_g) ];
    sumOfWeights += _prob.weights[n];
  }

  return scalarFlux/sumOfWeights;
}

