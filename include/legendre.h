#ifndef LEGENDRE_H
#define LEGENDRE_H

/// Evaluate Legendre polynomials
/**
 *  Each object represents a set of Legendre polynomials from degree
 *  0 to _maxDegree evaluated a given number.  The evaluations can
 *  be accessed using square brackets, i.e. polynomial[degree].
 **/
class LegendrePolynomial
{
 public:
  LegendrePolynomial(int maxDegree);
  ~LegendrePolynomial();

  void evaluate(double x);
  
  double operator[](int degree) { return _polyVals[degree]; };

 private:
  int _maxDegree;
  double* _polyVals;
};

#endif
