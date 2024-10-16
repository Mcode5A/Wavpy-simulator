/*#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <memory>
*/

#include <gsl/gsl_complex.h>
#include "ComplexWaveform.hpp"
 
int main(int argc, char* argv[]) {

  int N = 5;
  double re[N] = {1, 2, 3,  4,  5};
  double im[N] = {7, 8, 9, 10, 11}; 
  ComplexWaveform a(N, re, im);
  ComplexWaveform b(N, re, im);
  gsl_complex cmplx = gsl_complex_rect(1, -2);
  double dbl = 6.432;

  ComplexWaveform d;
  d = a+b;
  d = a+cmplx;

  d = a-b;
  d = a-cmplx; 
  
  d = a*b;
  d = a*dbl;

  d = a.Conj();

  d = a/b;
  d = a/dbl;

  return 0;
}
