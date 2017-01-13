#include "legendre.h"

LegendrePolynomial::LegendrePolynomial(int maxDegree)
{
  _maxDegree = (maxDegree > 1)? maxDegree : 1;
  
  _polyVals = new double [maxDegree+1];
  _polyVals[0] = 1.0;
  for (int n=1; n<=_maxDegree; n++) {
    _polyVals[n] = 0.0;
  }
}

LegendrePolynomial::~LegendrePolynomial()
{
  delete [] _polyVals;
}

void
LegendrePolynomial::evaluate(double x)
{
  _polyVals[1] = x;
  
  for (int n=2; n<=_maxDegree; n++) {
    _polyVals[n] = ( (2.0*n - 1.0)*x*_polyVals[n-1] - (n-1.0)*_polyVals[n-2] )/n;
  }
}
