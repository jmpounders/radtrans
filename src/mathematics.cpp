#include "mathematics.h"

namespace math
{
  /// Computes the double factorial
  int doubleFactorial(int n)
  {
    int k;
    if (n%2 == 0) {
      k = n/2;
      return std::pow(2,k)*factorial(k);
    }
    else {
      k = (n+1)/2;
      return factorial(2*k)/(std::pow(2,k)*factorial(k));
    }
        
  }
}
