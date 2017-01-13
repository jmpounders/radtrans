#ifndef ASSOCIATEDLEGENDRE_H
#define ASSOCIATEDLEGENDRE_H

namespace math{

  /// Evaluate associated Legendre polynomials
  class AssociatedLegendrePolynomial
  {
   public:
    AssociatedLegendrePolynomial(int maxDegree);
    ~AssociatedLegendrePolynomial();
    
    void evaluate(double x);

    double operator()(int degree, int order) { return _polyVals[degree]; };
  
   private:
    int _maxDegree;
    int _numTerms;
    double* _polyVals;

    int _ind(int degree, int order) { return degree*degree + degree + order; };
  };
 
}
 
#endif
