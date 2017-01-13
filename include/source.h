#ifndef SOURCE_H
#define SOURCE_H

/// Source object
/**
 * The source object provides a generic source term for the transport equation.
 */
class Source
{
 public:
  Source(long nDOFS_in);
  ~Source();

 private:
  long _nDOFs;
  double* _sourceTerm;
};

#endif
