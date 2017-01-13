#ifndef TRANSPORTPROBLEM_H
#define TRANSPORTPROBLEM_H

#include<vector>

#include "meshfactory.h"
#include "material.h"
#include "inputparser.h"

/// Boundary condition types
enum TransportBC {vacuum, reflecting, source};

/// Transport problem object
/**
 *  The transport problem object encapsulates (1) a mesh pointer,
 *  (2) high-level descriptive information, (3) boundary conditions,
 *  (4) angular quadrature, and (5) external source information.
 *
 * TODO: generalize boundary conditions
 *       generalize quadrature
 */
class TransportProblem
{
 public:
  TransportProblem(InputParser& input, MeshFactory& meshFactory, int numGroups_in);
  ~TransportProblem();
  MeshInterface *mesh;
  long numCells;
  long numEdges;
  long numNodes;
  int numGroups;
  int quadOrder;
  int scatterAnisotropy;

  bool hasFixedSource;

  TransportBC bcLeft;
  TransportBC bcRight;

  TransportBC globalBC;

  double* extSource;  // needs to become private
  double* refSolution;

  std::vector<double> omega_x;
  std::vector<double> omega_y;
  std::vector<double> omega_z;
  std::vector<double> weights;
  
  std::vector<double> nxnyq;

  void putExtSource(long cell, int quadp, int group, double value)
  { extSource[ si(cell,quadp,group) ] = value; }
  double* getExtSource(long cell, int quadp, int group)
  { return &extSource[ si(cell,quadp,group) ]; }

  void putRefSolution(long cell, int quadp, int group, double value)
  { refSolution[ si(cell,quadp,group) ] = value; }
  double* getRefSolution(long cell, int quadp, int group)
  { return &refSolution[ si(cell,quadp,group) ]; }

 private:
  // Fast, medium, slow
  // space, angle, energy
  long si(long cell, int quadp, int group)
  { return numCells*quadOrder*group + numCells*quadp + cell; }
};
#endif
