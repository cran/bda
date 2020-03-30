/*
// ISNA(x): true for R's NA only
// ISNAN(x): true for R's NA and IEEE NaN
// R_FINITE(x): false for Inf,-Inf,NA,NaN
// R_IsNaN(x): true for NaN but not NA

// Rprintf:  printing from a C routine compiled into R

 */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"

void KSPvalue(double *x0){
  double ksp=0., t=x0[0];
  int i;
  for(i=1;i<100;i++){
    ksp += exp(-2.*pow(i*t,2.));
    i++;
    ksp -= exp(-2.*pow(i*t,2.));
  }
  x0[0] = 2.*ksp;
}


void pks2(double *x, int *size1, int *size2){
  double md, nd, q, *u, w;
  int i, j, m=size1[0], n=size2[0];
  if(m > n) {
    i = n; n = m; m = i;
  }
  md = (double) m;
  nd = (double) n;
  q = (0.5 + floor(*x * md * nd - 1e-7)) / (md * nd);
  u = (double *) R_alloc(n + 1, sizeof(double));
  for(j = 0; j <= n; j++) {
    u[j] = ((j / nd) > q) ? 0 : 1;
  }
  for(i = 1; i <= m; i++) {
      w = (double)(i) / ((double)(i + n));
      if((i / md) > q)
	u[0] = 0;
      else
	u[0] = w * u[0];
      for(j = 1; j <= n; j++) {
	if(fabs(i / md - j / nd) > q) 
	  u[j] = 0;
	else
	  u[j] = w * u[j] + u[j - 1];
      }
  }
  x[0] = fabs(1.0 - u[n]);
}



void KSP2x(double *D, int *size){
  /* copied from R ks.c file */
  
  /* Compute Kolmogorov's distribution.
     Code published in
     George Marsaglia and Wai Wan Tsang and Jingbo Wang (2003),
     "Evaluating Kolmogorov's distribution".
     Journal of Statistical Software, Volume 8, 2003, Issue 18.
     URL: http://www.jstatsoft.org/v08/i18/.
  */

  int n=size[0];
  double d=D[0];
  
  double p=0.0, s;
   
  /* 
     The faster right-tail approximation.
  */
  s = d*d*n; 
  if(s > 7.24 || (s > 3.76 && n > 99)){ 
    p = 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
  }
  D[0] = p;
}
