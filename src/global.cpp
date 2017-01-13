#include "global.h"

float
norm(double* vec, int length)
{
  float sum = 0.0;
  for (int i=0; i<length; i++) {
    sum += pow(vec[i], 2);
  }
  return sqrt(sum);
}
