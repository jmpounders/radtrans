#include "associatedlegendre.h"
#include "mathematics.h"

namespace math{

  AssociatedLegendrePolynomial::AssociatedLegendrePolynomial(int maxDegree)
    : _maxDegree(maxDegree)
  {
    _numTerms = (_maxDegree + 1)*(_maxDegree + 1);

    _polyVals = new double [_numTerms];
    _polyVals[0] = 1.0;
    for (int i=1; i<=_numTerms; i++)
      _polyVals[i] = 0.0;
  }

  
  AssociatedLegendrePolynomial::~AssociatedLegendrePolynomial()
  {
    delete [] _polyVals;
  }

  
  void
  AssociatedLegendrePolynomial::evaluate(double x)
  {
    if (_maxDegree==0) return;

    _polyVals[1] = 0.5*(std::sqrt(1.0 - std::pow(x,2)));
    _polyVals[2] = x;
    _polyVals[3] = -std::sqrt(1.0 - std::pow(x,2));
    
    for (int ell=2; ell<_maxDegree; ell++) {
      for (int m=2-ell; m<ell-2; m++) {
        _polyVals[ _ind(ell,m) ] = ((2*ell-1)*x*_polyVals[ _ind(ell-1,m) ] - (ell+m-1)*_polyVals[ _ind(ell-2,m) ])/(ell-m);
      }
      _polyVals[ _ind(ell,ell) ] = -(2*ell - 1)*std::sqrt(1 - std::pow(x,2))*_polyVals[ _ind(ell-1,ell-1) ];
      _polyVals[ _ind(ell,ell-1) ] = x*(2*ell - 1)*_polyVals[ _ind(ell-1,ell-1) ];
      _polyVals[ _ind(ell,1-ell) ] = std::pow(-1,ell-1)/factorial(2*ell-1)*_polyVals[ _ind(ell,ell-1) ];
      _polyVals[ _ind(ell,-ell) ] = std::pow(-1,ell)/factorial(2*ell)*_polyVals[ _ind(ell,ell) ];
    }
  }
  
}
