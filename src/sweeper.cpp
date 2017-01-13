#include "sweeper.h"

Sweeper::Sweeper(MoabMesh* mesh_in, int angle) : mesh(mesh_in)
{
  _sweepOrder.clear();
  if (angle < 0)
    mesh->_mb->get_entities_by_dimension(0, 2, _sweepOrder);
  else
    mesh->_mb->get_entities_by_dimension(mesh->_orderedMeshSet[angle], 2, _sweepOrder);

  _nextElementIterator = _sweepOrder.begin();
}


long
Sweeper::getNextElementID()
{
  if (_nextElementIterator == _sweepOrder.end()) {
    return -1;
  }
  else {
    long elementID = mesh->_mb->id_from_handle(*_nextElementIterator)-1;
    _nextElementIterator++;
    return elementID;
  }
}




