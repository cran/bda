#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*
c  Part of R package BDA
c  Copyright (C) 2009-2010 Bin Wang
c
c  Unlimited use and distribution (see LICENCE).

c     This program is to compute the MLE of Gamma(k,theta). (1) compute
c     s=log(mean(x))-mean(log(x)).  (2) compute initial value k =
c     (3-s+sqrt((s-3)^2+24s))/(12s) (3) use Newton method find a
c     numerical solution k1 = k - (log(k)-psi(k)-s)/(1/k-psi'(k)), where
c     psi(.) is the digamma function, and psi'(.) is the trigamma
c     function. [S. C. Choi and R. Wette. (1969) Maximum Likelihood
c     Estimation of the Parameters of the Gamma Distribution and Their
c     Bias, Technometrics, 11(4) 683–690]

c     Last changed: 17 Feb 2011


*/

void remlegamma(double *x,double *f, double *b,int *size, double *theta)
{
  int n=size[0],i,iter;
  double s0,s,k0,k1,t1,t2,xbar,llk;
  
  llk = -999.0;
  k0 = theta[0];
  s = theta[1];
  iter=0;
  t2 = 100.;
  while(iter<100 && t2 > 0.000001){
      k1 = k0;
      k0 -= (log(k0)-digamma(k0)-s)/(1./k0-trigamma(k0));
      iter = iter + 1;
      t1 = fabs(k1-k0);
      t2 = fmax(t1,fabs(t1/fmin(k1,k0)));
  }
  size[0] = iter; // store the iteration number in computing kappa.
  theta[0] = k0;//k0;// first parameter kappa
}

