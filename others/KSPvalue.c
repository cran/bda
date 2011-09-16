#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void RKSPvalue(double *t0){
  double ksp=0.,t=t0[0];
  int i;
  for(i=1;i<100;i++){
    ksp += exp(-2.*pow(i*t,2.));
    i++;
    ksp -= exp(-2.*pow(i*t,2.));
  }
  t0[0] = 2.0 * ksp;
}
