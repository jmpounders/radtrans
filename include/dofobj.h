#ifndef DOFOBJ_H
#define DOFOBJ_H

/// DOF Object
/**  
 *  A DOF object represents a geometric location, angular quadrature point,
 *  and energy group.
 */
class DOFObj
{
 public:
  DOFObj(long geomID, int angQuadPnt, int eGrp, int si=0);
  long geometricID;
  int angularQuadraturePoint;
  int energyGroup;
  int subIndex;

  bool operator<(const DOFObj& dofIn)const; //!< Less-than operator so that these things can be mapped
};

#endif
