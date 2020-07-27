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

void subdKDE(double y0[], double x0[], int n0, 
	  double x[], double h[], double f[], int n)
{
  int i,j;
  double xi,nsum=0.0;
  for(j=0; j<n;j++){
    nsum += f[j];
  }

  for(i=0; i<n0;i++){
    y0[i] = 0.0;
    for(j=0; j<n;j++){
      xi = (x0[i] - x[j])/h[j];
      y0[i] += dnorm(xi,0.,1.,0)/h[j]*f[j];
    }
    y0[i] /= nsum;
  }
}

void subpKDE(double y0[], double x0[], int n0, 
	  double x[], double h[], double f[], int n)
{
  int i,j;
  double xi,nsum=0.0;
  for(j=0; j<n;j++){
    nsum += f[j];
  }

  for(i=0; i<n0;i++){
    y0[i] = 0.0;
    for(j=0; j<n;j++){
      xi = (x0[i] - x[j])/h[j];
      y0[i] += pnorm(xi,0.,1.,1,0)/h[j]*f[j]; 
    }
    y0[i] /= nsum;
  }
}

void dKDE(double *x0, int *n0, 
	  double *x, double *h, double *f, int *n)
{
  int i,k=n0[0];
  double y[k];
  for(i=0; i<k;i++){
    y[i] = 0.0;
  }

  subdKDE(y,x0,k,x,h,f,n[0]);

  for(i=0; i<k;i++){
    x0[i] = y[i];
  }
}

void pKDE(double *x0, int *n0, 
	  double *x, double *h, double *f, int *n)
{
  int i,k=n0[0];
  double y[k];
  for(i=0; i<k;i++){
    y[i] = 0.0;
  }

  subpKDE(y,x0,k,x,h,f,n[0]);

  for(i=0; i<k;i++){
    x0[i] = y[i];
  }
}
