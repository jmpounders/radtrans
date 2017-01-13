#ifndef SWEEPER_H
#define SWEEPER_H

#include "moabmesh.h"

/// Mesh sweeper
/**
 *   This object will sweep a (Moab) mesh for a given angle, n.
 *   This was created to make parallel sweeps easier, i.e., make
 *   one object for each thread.
 **/
class Sweeper
{
 public:
  Sweeper(MoabMesh* mesh, int angle);
  void startElementSweep(int angle = -1);
  long getNextElementID();

 private:
  MoabMesh* mesh;
  std::vector<moab::EntityHandle>::const_iterator _nextElementIterator;
  std::vector<moab::EntityHandle> _sweepOrder;

  UltraLightElement _currentElementUL;
  void _setupElementUL();

  moab::ErrorCode _rval;
};


#endif
