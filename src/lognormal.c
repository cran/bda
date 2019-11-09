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

void mclnorm(double *x, double *fn,
	     int *size, double *mu, double *s)
{
  int i, j, k, iter = 50;
  double D0, D1, Dk, mu0, s0, step1, step2, mu1, s1, Fx;
  
  step1 = s[0] * 4.0/iter;
  step2 = s[0] * 10.0/iter;
  D1 = 2.0;
  mu0 = mu[0] - 2.0 * s[0]; //mu-2s to mu+2s
  mu1 = mu0; 
  s0 = s[0] * 5.0/iter; //s0 yo 5*s0
  s1 = s0;
  for(i=0; i<iter; i++){
    s0 = s[0] * 0.01;
    for(j=0; j<iter; j++){
      D0 = 0.0;
      for(k=0; k<size[0]; k++){
	Fx = plnorm(x[k], mu0, s0,1,0);
	Dk = fabs(Fx-fn[k]);
	if(Dk > D0){
	  D0 = Dk;
	}
      }
      if(D0 < D1){
	D1 = D0;
	mu1 = mu0;
	s1 = s0;
      }
      s0 += step2;
    }
    mu0 += step1;
  }
  mu[0] = mu1;
  s[0] = s1;
}

double LlkTN(double x[], double d[], double f[], int n,
	      double xp, double qp, double sig)
{
  int i;
  double mu, llk,F2,F1,tmp;

  mu = xp - sig * qp;
  
  llk = 0.0;
  for(i=0; i<n; i++){
    F1 = pnorm(x[i], mu, sig, 1, 0);
    F2 = pnorm(x[i]+d[i], mu, sig, 1, 0);
    tmp = fabs(F2-F1);
    if(tmp < 0.00000001){
      llk += f[i] * (-10.00);
    }else{
      llk += f[i] * log(tmp);
    }
  }
  return llk;
}

double dfLlkTN(double x[], double d[], double f[], int n,
	      double xp, double qp, double sig)
{
  int i;
  double mu, llk,F2,F1,f1,f2;
  
  mu = xp - sig * qp;
  
  llk = 0.0;
  for(i=0; i<n; i++){
    F1 = pnorm(x[i], mu, sig, 1, 0);
    F2 = pnorm(x[i]+d[i], mu, sig, 1, 0);
    f1 = dnorm(x[i], mu, sig, 0);
    f2 = dnorm(x[i]+d[i], mu, sig, 0);
    llk += f[i] * (f2-f1)/(F2-F1);
  }
  return llk;
}

void mleTN(double *x, double *d, double *f, int *n,
	   double *xp, double *qp, double *sig)
{
  int i, iter=1000;
  double sa, sb, sc, fa,fb,fc,tol=1.0e-5;
  double xa, fxa;
  sa = sig[0] * 0.01;
  sb = sig[0] * 10;
  sc = (sa + sb)*0.5;
  
  fa = LlkTN(x,d,f,n[0],xp[0],qp[0],sa);
  fb = LlkTN(x,d,f,n[0],xp[0],qp[0],sb);
  fc = LlkTN(x,d,f,n[0],xp[0],qp[0],sc);
  if(fa > fc){
    //    Rprintf("Solution is smaller than %10.3f.\n", sa);
    n[0] = -1;
  }else if(fb > fc){
    //    Rprintf("Solution is greater than %10.3f.\n", sb);
    n[0] = -2;
  }else{
    for (i=0; i < iter; i++) {
      xa = (sa + sc) * 0.5;
      fxa = LlkTN(x,d,f,n[0],xp[0],qp[0],xa);
      if(fxa < fc){
	sa = xa; fa = fxa;
      }else{
	sb = sc; fb = fc;
	sc = xa; fc = fxa;
      }
      if(fabs(sb-sa) < tol){
	n[0] = 0;
	sig[0] = sc;
	break;
      }
    }
  }
}

double CompSig(double x[], double d[], double f[], int n,
	       double xp, double qp, double sig)
{
  int i, iter=1000;
  double sa, sb, sc, fa,fb,fc,tol=1.0e-5;
  double xa, fxa, s=0.0;
  sa = sig * 0.01;
  sb = sig * 10.0;
  sc = (sa + sb)*0.5;
  
  fa = LlkTN(x,d,f,n,xp,qp,sa);
  fb = LlkTN(x,d,f,n,xp,qp,sb);
  fc = LlkTN(x,d,f,n,xp,qp,sc);

  //  Rprintf("fa= %10.3f, fb= %10.3f, fc= %10.3f.\n", fa, fb, fc);

  if(fa == fc){
    s = sc;
  }else if(fb == fc){
    s = sc;
  }else if((fb-fc)*(fc-fa)>0.0){
    for (i=0; i < iter; i++) {
      xa = (sa + sc) * 0.5;
      fxa = LlkTN(x,d,f,n,xp,qp,xa);
      if(fxa < fc){
	sa = xa; fa = fxa;
      }else{
	sb = sc; fb = fc;
	sc = xa; fc = fxa;
      }
      if(fabs(sb-sa) < tol){
	s = sc;
	break;
      }
    }
  }
  
  return s;
}

void mlemixTN(double *x, double *d, double *f, int *n,
	      double *xp, double *qp, double *sig,
	      double *pmix, int *k)
{
  int i,j,iter=100,m;
  double F1, F2, s, mu[k[0]], fk[n[0]], w[n[0]][k[0]];
  double wsum, fsum=0.0, delta=1.0, pnew=0.0,tol=1.0e-5;
  double s2[k[0]],mu2[k[0]],p2[k[0]];
  //initialize fk and w;
  for (i=0; i < n[0]; i++) {
    fk[i] = 0.0;
    fsum += f[i];
    for (j=0; j < k[0]; j++) {
      w[i][j] = 0.0;
    }
  }
  
  
  s = CompSig(x,d,f,n[0],xp[0],qp[0],sig[0]);
  //initialize the parameters
  
  for (i=0; i < k[0]; i++) {
    sig[i] = s * pow(1.1,1.0*i);
    mu[i] = xp[0] - sig[i] * qp[0];
    s2[i] = sig[i];
    mu2[i] = mu[i];
    p2[i] = pmix[i];
  }

  for (m=0; m < iter; m++) {

    for (j=0; j < k[0]; j++) {
      pmix[j] = 0.0;
    }
    
    for (i=0; i < n[0]; i++) {
      wsum = 0.0;
      for (j=0; j < k[0]; j++) {
	F1 = pnorm(x[i], mu[j], sig[j], 1, 0);
	F2 = pnorm(x[i]+d[i], mu[j], sig[j], 1, 0);
	w[i][j] = F2-F1;
	wsum += F2-F1;
      }
      for (j=0; j < k[0]; j++) {
	w[i][j] = w[i][j] / wsum * f[i];
	pmix[j] += w[i][j];
      }
    }
    for (j=0; j < k[0]; j++) {
      pmix[j] /= fsum;
    }

    //update estimate of sigma and mu
    delta = 0.0;
    for (j=0; j < k[0]; j++) {
      // get the new frequencies
      for (i=0; i < n[0]; i++) {
	fk[i] = w[i][j];
      }
      sig[j] = CompSig(x,d,fk,n[0],xp[0],qp[0],sig[j]);

      //      Rprintf("sig=%10.3f,j=%3d\n",sig[j],j);
      
      pnew = fabs(sig[j] - s2[j]);
      if(pnew > delta){
	delta = pnew;
      }
      s2[j] = sig[j];
      
      mu[j] = xp[0] - sig[j] * qp[0];
      pnew = fabs(mu[j] - mu2[j]);
      if(pnew > delta){
	delta = pnew;
      }
      mu2[j] = mu[j];
      
      pnew = fabs(pmix[j] - p2[j]);
      if(pnew > delta){
	delta = pnew;
      }
      p2[j] = pmix[j];
    }
    if(delta > tol){
      n[0] = 0;
      break;
    }
  }
}


void qtlmlnorm(double *q, int *k,
	       double *p, double *mu, double *s)
{
  int i, j, iter=100000;
  double x0=10.0, dx, Fx, fx, tol=1.0e-5;

  for (i=0; i < iter; i++) {
    fx = 0.0;
    Fx = 0.0;
    for(j=0; j < k[0]; j++){
      fx = p[j] * dlnorm(x0, mu[j], s[j],0);
      Fx = p[j] * plnorm(x0, mu[j], s[j],1,0);
    }
    if(fx < 1.e-5){
      dx = -x0;
      x0 *= 2.0;
    }else{
      dx = (Fx-q[0]) / fx;
      if(x0 < dx){
	x0 *= 0.5;
      }else{
	x0 -= dx;
      }
    }

    if (fabs(dx) < tol){
      q[0] = x0;
      k[0] = 0;
      break;
    }
  }
}


/* 
The following algorithm fit a lognormal distribution to NGS data
truncated at 1. We find the best fit by matching the tail shapes --
min(distance(fn, E(fn))).
 */

void mclnorm2(double *x, double *fn, double *delta,
	      int *size, double *mu, double *s)
{
  int i, j, k, iter = 50;
  double D0, D1, mu0, s0, step1, step2, mu1, s1;
  double Fx[size[0]-1], n;

  n = 0.0;
  for(i=0; i<size[0]; i++){
    n += fn[i];
  }
  for(i=1; i<size[0]; i++){
    Fx[i-1] = 0.0;
  }
  
  step1 = s[0] * 4.0/iter;
  step2 = s[0] * 10.0/iter;
  D1 = 99999999999.0;
  mu0 = mu[0] - 2.0 * s[0]; //mu-2s to mu+2s
  mu1 = mu0; 
  s0 = s[0] * 5.0/iter; //s0 yo 5*s0
  s1 = s0;
  
  for(i=0; i<iter; i++){
    s0 = s[0] * 0.01;
    for(j=0; j<iter; j++){
      D0 = 0.0;
      for(k=1; k<size[0]; k++){
	Fx[k-1] = plnorm(x[k]+delta[k], mu0, s0,1,0) -
	  plnorm(x[k], mu0, s0,1,0);
	Fx[k-1] *= n;
	D0 += fabs(fn[k] - Fx[k-1]);
      }

      if(D0 < D1){
	D1 = D0;
	mu1 = mu0;
	s1 = s0;
      }
      s0 += step2;
    }
    mu0 += step1;
  }
  mu[0] = mu1;
  s[0] = s1;
  fn[0] = plnorm(x[1], mu1, s1,1,0);
}

// lognormal based on binned data.  Using MLE by search over [mu0,mu1]
// and [0.5s, 3s]

void lnormBinMLE(int *size, double *x, double *fn,
		 double *mu, double *s)
{
  int i, j, k, iter = 1000;
  double D0, D1,Dmax, mu0, s0, step1, step2, mu1, s1, F0, Fx;

  step1 = 3.0 * mu[0] / iter;
  step2 = 2.5 * s[0] / iter;
  D1 = -999999999999999.0;
  Dmax = D1;
  mu0 = mu[0] * 0.01;
  mu1 = mu0;   //best solution
  s0 = s[0] * 0.5 / iter; // 0.5s0 to 3*s0
  s1 = s0;     //best solution
  //  printf("mu0=%10.2f, s0=%10.2f.\n", mu0, s0);
  for(i=0; i<iter; i++){
    s0 = s[0] * 0.5 / iter; // 0.5s0 to 3*s0
    for(j=0; j<iter; j++){
      D0 = 0.0;
      F0 = 0.0;
      for(k=0; k<size[0]-1; k++){
	Fx = plnorm(x[k], mu0, s0,1,0);
	//	printf("Fx %10.2f.\n", Fx);
	if(Fx - F0 > 0.0){
	  D0 += fn[k] * log(fabs(Fx - F0));
	}else{
	  D0 += Dmax;
	}
	F0 = Fx;
      }
      if(F0 < 1.0){
	D0 += fn[size[0]] * log(1.0-F0);
      }else{
	D0 += Dmax;
      }
      //      printf("LLK %10.2f, mu=%10.2f, s=%10.2f.\n", D0, mu0, s0);
      if(D0 > D1){
	D1 = D0;
	mu1 = mu0;
	s1 = s0;
      }
      
      s0 += step2;
    }
    mu0 += step1;
  }
  mu[0] = mu1; 
  s[0] = s1;
}



void lnormBinChisq(int *size, double *x, double *fn,
		 double *mu, double *s)
{
  int i, j, k, iter = 1000;
  double D0, D1, mu0, s0, step1, step2, mu1, s1, F0, Fx;
  double ne, nsum;

  nsum = 0.0;
  for(i=0; i<=iter; i++){
    nsum += fn[i];
  }

  step1 = (5.0*mu[1] - mu[0])/iter;
  step2 = s[0] * 2.5/iter;
  D1 = 999999999999999.0;
  mu0 = mu[0]; //mu0 -- mu1
  mu1 = mu0;   //best solution
  s0 = s[0] * 0.5/iter; // 0.5s0 to 3*s0
  s1 = s0;     //best solution
  for(i=0; i<iter; i++){
    for(j=0; j<iter; j++){
      D0 = 0.0;
      F0 = 0.0;
      for(k=0; k<size[0]; k++){
	Fx = plnorm(x[k], mu0, s0,1,0);
	ne = (Fx - F0) * nsum;
	D0 += pow(fn[k] - ne, 2.0)/ne;
	F0 = Fx;
      }
      ne = (1.0 - F0) * nsum;
      D0 += pow(fn[k] - ne, 2.0)/ne;
      if(D0 < D1){
	D1 = D0;
	mu1 = mu0;
	s1 = s0;
      }
      s0 += step2;
    }
    mu0 += step1;
  }
  mu[0] = mu1; mu[1] = D1;
  s[0] = s1;
}
