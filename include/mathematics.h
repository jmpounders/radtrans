#ifndef MATHEMATICS_H
#define MATHEMATICS_H

#include <cmath>

namespace math
{
  inline int factorial(int n) { return (int) std::tgamma(n+1); };
  int doubleFactorial(int n);
}

#endif
