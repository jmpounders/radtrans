#include<cmath>
#include<iostream>
#include<iomanip>

#include "solverbase.h"
#include "perfstats.h"
#include "global.h"

/**
 *  Define number of DOF, map DOFs, allocate solution vectors
 */
SolverBase::SolverBase(TransportProblem &tp) :
  _prob(tp), sourceScaling(1), criticalEigenvalue(1), _convRInfTol(1.0e-6), _maxIters(1000)
{
}

/**
 *  Free memory
 */
SolverBase::~SolverBase()
{
}

/**
 *  Calculate the residual (this should be implemented in the base)
 */
void
SolverBase::_calculateResidual()
{
  _calculateMatrixAction( _solution, _residual);
  
  for (int _dof_i=0; _dof_i<_numDOF; _dof_i++) {
      _residual[_dof_i] = _residual[_dof_i] - _source[_dof_i];
  }

}

/**
 *  Calculate the solution convergence norm
 */
double
SolverBase::calculateL2SolutionNorm(double* solutionComp)
{
  int _dof_i = 0;
  double sumSqs = 0.0;

  for (int g=0; g<_prob.numGroups; g++) {
    for (long i=0; i<_prob.numCells; i++) {
      for (int n=0; n<_prob.quadOrder; n++) {
	sumSqs += (_solution[_dof_i] - solutionComp[_dof_i])*
	          (_solution[_dof_i] - solutionComp[_dof_i]);
	_dof_i++;
      }
    }
  }

  return sqrt(sumSqs);
}

/**
 *  Calculate the solution convergence norm
 */
double
SolverBase::calculateLInfSolutionNorm(double* solutionComp)
{
  long _dof_i = 0;
  double maxDiff = 0.0;

  for (int g=0; g<_prob.numGroups; g++) {
    for (long i=0; i<_prob.numCells; i++) {
      for (int n=0; n<_prob.quadOrder; n++) {
	maxDiff = fmax(maxDiff, std::abs(_solution[_dof_i] - solutionComp[_dof_i]));
	_dof_i++;
      }
    }
  }

  return maxDiff;
}


/**
 *  Calculate the solution convergence norm
 */
double
SolverBase::calculateRInfSolutionNorm(double* solutionComp)
{
  long _dof_i = 0;
  double maxDiff = 0.0;

  for (int g=0; g<_prob.numGroups; g++) {
    for (long i=0; i<_prob.numCells; i++) {
      for (int n=0; n<_prob.quadOrder; n++) {
        if (_solution[_dof_i] > 0.0) {
          maxDiff = fmax(maxDiff, std::abs((_solution[_dof_i] - solutionComp[_dof_i])/_solution[_dof_i]));
        }
	_dof_i++;
      }
    }
  }

  return maxDiff;
}

/**
 *  Calculate the solution convergence norm
 */
double
SolverBase::calculateL2SolutionError()
{
  long _dof_i = 0;
  double sumSqs = 0.0;

  for (int g=0; g<_prob.numGroups; g++) {
    for (long i=0; i<_prob.numCells; i++) {
      for (int n=0; n<_prob.quadOrder; n++) {
	sumSqs += (_solution[_dof_i] - *_prob.getRefSolution(i,n,g))*
	          (_solution[_dof_i] - *_prob.getRefSolution(i,n,g));
	_dof_i++;
      }
    }
  }

  return sqrt(sumSqs);
}


/**
 *  Calculate the solution convergence norm
 */
double
SolverBase::calculateRInfSolutionError()
{
  long _dof_i = 0;
  double maxDiff = 0.0;

  for (int g=0; g<_prob.numGroups; g++) {
    for (long i=0; i<_prob.numCells; i++) {
      for (int n=0; n<_prob.quadOrder; n++) {
	maxDiff = fmax(maxDiff, std::abs((_solution[_dof_i] - *_prob.getRefSolution(i,n,g))/_solution[_dof_i]));
	_dof_i++;
      }
    }
  }

  return maxDiff;
}

/**
 *  Save _solution to _solutionPrev
 */
void
SolverBase::_saveOldSolution()
{
  for (long i=0; i<_numDOF; i++) {
    _solutionPrev[i] = _solution[i];
  }
}

double
SolverBase::getNeutronProduction(long space_i)
{
  PerfStats X("SolverBase::getNeutronProduction");
  double productionRate = 0;
  for (int g=0; g<_prob.numGroups; g++) {
    productionRate += *_prob.mesh->getElementMat(space_i)->getNuSigma_f(g+1)
                      * getScalarFlux(space_i, g)*4.0*pi;
  }
  return productionRate;
}


double
SolverBase::getTotalNeutronProduction()
{
  PerfStats X("SolverBase::getTotalNeutronProduction");
  double totalProductionRate = 0;
  double meshVolume;
  for (long i=0; i<_prob.numCells; i++) {
    meshVolume = _prob.mesh->getElementVolume(i);
    for (int g=0; g<_prob.numGroups; g++) {
      totalProductionRate += *_prob.mesh->getElementMat(i)->getNuSigma_f(g+1)
                             * getScalarFlux(i, g)*4.0*pi
                             * meshVolume;
    }
  }
  return totalProductionRate;
}


void
SolverBase::setSolution(double* solution)
{
  for (long i=0; i<_numDOF; i++) {
    _solution[i] = solution[i];
  }
}


double*
SolverBase::copySolution(double* solutionCopy)
{
  if (!solutionCopy) {
    solutionCopy = new double [_numDOF];
  }
  for (long i=0; i<_numDOF; i++) {
    solutionCopy[i] = _solution[i];
  }
  return solutionCopy;
}


void
SolverBase::normalizeSolution(double totalProductionRate)
{
  double norm = totalProductionRate/getTotalNeutronProduction();
  for (long i=0; i<_numDOF; i++) {
    _solution[i] = _solution[i] * norm;
  }  
}


void
SolverBase::zeroSolution()
{
  for (long i=0; i<_numDOF; i++) {
    _solution[i] = 0.0;
  }
}


void
SolverBase::extrapolateSolution(double expExtrapFactor)
{
  for (long i=0; i< _numDOF; i++) {
    _solution[i] = _solution[i]*exp(expExtrapFactor);
  }
}

void
SolverBase::extrapolateSolution(double expExtrapFactor, double* baseSolution)
{
  for (long i=0; i< _numDOF; i++) {
    _solution[i] = baseSolution[i]*exp(expExtrapFactor);
  }
}

void
SolverBase::writeSolution(Output* outputFile)
{
  outputFile->writeData("numDOFs", &_numDOF, 1);
  //outputFile->writeData("solution", _solution, _numDOF);
  //if (outputFile->outputFormat & Output::MESH)
  //  _prob.mesh->writeMesh(outputFile);
}

void
SolverBase::writeScalarFlux(Output* outputFile)
{
  // Should write this to a "Solver" group
  long outputSize = _numSpaceDOF;
  double* outputBuffer = new double [ outputSize ];
  std::ostringstream varNameStream;
  std::string varName;

  LOG_DBG("writing ", outputSize ," output data.");

  // Write flux vector
  for (int g=0; g<_prob.numGroups; g++) {
    for (long i=0; i< _numSpaceDOF; i++) {
      outputBuffer[i] = getScalarFlux( i, g );
    }
    varNameStream.str("");
    varNameStream << "scalarFlux_" << g+1;
    varName = varNameStream.str();
    //outputFile->writeData(varName, outputBuffer, outputSize);
    if (outputFile->outputFormat & Output::MESH) {
      _prob.mesh->tagMesh(varName, outputBuffer, _numSpaceDOF);
      _prob.mesh->writeMesh(outputFile);
    }
  }

  delete [] outputBuffer;
}

void
SolverBase::printIterStatus(std::string name, int iter, double err, double tol)
{
  std::stringstream outs;

  outs << name << ": ";
  outs << std::setw(6) << iter;
  outs << std::scientific;
  outs << "   " << err;
  outs << "   " << tol;
  LOG(outs.str());

}
