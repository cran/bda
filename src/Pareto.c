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

void qmPareto(double *p, double *q, int *npar,
	      double *xm, double *alpha)
{
  int i,j,n=npar[0], k;

  k = 0;
  for(i=0; i<n-1; i++){
    for(j=i+1; j<n; j++){
      alpha[k] = log((1.0-p[i])/(1.0-p[j]))/log(q[j]/q[i]);
      //Rprintf("\nalpha=%f, ",alpha[k]);
      if(alpha[k] <= 0.0){
	xm[k] = -99.;
      }else{
	xm[k] = pow(1.0-p[i], 1.0/alpha[k]) * q[i];
      }
      //Rprintf("\nxm=%f",xm[k]);
      k++;
    }
  }
}

double binParetoLLK(double f[], double b[], int n,
		    double xm, double a){
  int i;
  double p,q, llk=0.0;

  p = 1.0 - pow(xm/b[0], a);
  if(p>0.0){
    llk += log(p)*f[0];
  }else{
    llk -= 999. * f[0];
  }
  for(i=1; i<n-1; i++){
    q = 1.0 - pow(xm/b[i], a);
    if(q > p){
      llk += log(q - p)*f[i];
      p = q;
    }else{
      llk -= 999. * f[0];
    }
  }
  if(p<1.0){
    llk += log(1.0 - p)*f[n-1];
  }else{
    llk -= 999. * f[0];
  }

  return llk;
}

void mle1Pareto(double *cnts, double *b, int *nclass,
	      double *xm, double *alpha)
{
  int i,n=nclass[0];
  double a=0.,s=0.,Fn[n],p=0.,q=0.;
  double llk,LLK,mle;
  //first find an initial value of alpha: Pareto index
  s = 0.0;
  for(i=0; i<n; i++){
    s += cnts[i];
    Fn[i] = s;
  }
  s = 1.0;
  for(i=0; i<n-1; i++){
    Fn[i] /= Fn[n-1];
    if(fabs(Fn[i] - 0.5) < s){
      p = Fn[i];
      q = b[i];
      s = fabs(Fn[i] - 0.5);
    }
  }
  a = log(1.0-p)/log(xm[0]/q);
  s = 0.05*a; //increment over [a/10, 10a]
  a = s; mle = a;
  LLK = binParetoLLK(cnts,b,n,xm[0],a);

  for(i=0; i<100; i++){
    a += s;
    llk = binParetoLLK(cnts,b,n,xm[0],a);
    if(llk > LLK){
      LLK = llk;
      mle = a;
    }
  }
  alpha[0] = mle;
  xm[0] = LLK;
}

void mle2Pareto(double *cnts, double *b, int *nclass,
	      double *xm, double *alpha)
{
  int i,j,n=nclass[0];
  double a=0.,s=0.,Fn[n],p=0.,q=0.,xm1;
  double llk,LLK,mle1,mle2;
  //first find an initial value of alpha: Pareto index
  s = 0.0;
  for(i=0; i<n; i++){
    s += cnts[i];
    Fn[i] = s;
  }
  s = 1.0;
  for(i=0; i<n-1; i++){
    Fn[i] /= Fn[n-1];
    if(fabs(Fn[i] - 0.5) < s){
      p = Fn[i];
      q = b[i];
      s = fabs(Fn[i] - 0.5);
    }
  }

  xm1 = xm[0]; mle2 = xm1;
  a = log(1.0-p)/log(xm1/q);
  a = s; mle1 = a;  
  LLK = binParetoLLK(cnts,b,n,xm1,a);

  s = 0.05*a; 
  for(j=0; j<99; j++){
    a = log(1.0-p)/log(xm1/q);
    s = a*0.1; a = s; 
    for(i=0; i<100; i++){
      llk = binParetoLLK(cnts,b,n,xm1,a);
      if(llk > LLK){
	LLK = llk;
	mle1 = a;
	mle2 = xm1;
      }
      a += s;
    }
    xm1 += xm[0];
  }
  alpha[0] = mle1;
  xm[0] = mle2;
  b[0] = LLK;
}

