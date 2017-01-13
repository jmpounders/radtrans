#include "dofobj.h"

DOFObj::DOFObj(long geomID, int qp, int g, int si) : 
  geometricID(geomID), angularQuadraturePoint(qp), energyGroup(g), subIndex(si)
{
}

/**
 *  The less-that comparison keys off space first, then quadrature point, then energy group.
 */
bool
DOFObj::operator<(const DOFObj& dofIn)const
{
  if (geometricID == dofIn.geometricID) {
    if (angularQuadraturePoint == dofIn.angularQuadraturePoint) {
      if (energyGroup == dofIn.energyGroup) {
	return (subIndex < dofIn.subIndex);
      }
      else return (energyGroup < dofIn.energyGroup);
    }
    else return (angularQuadraturePoint < dofIn.angularQuadraturePoint);
  }
  else return (geometricID < dofIn.geometricID);
}
