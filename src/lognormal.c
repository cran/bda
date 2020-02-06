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

// lognormal based on binned data.  Using MLE by search over [0.5mu0,1.5mu1]
// and [0.5s, 1.5s]

double lnormDist(double a[], double b[], double w[],
		 int n, double mu0, double s0){
  int i;
  double d1=0.0, d2=0.0, x1,x2, f1, f2, F1, F2;
  
  x1 = 0.0;
  x2 = log(b[0]) - mu0;
  f1 = 0.0;
  f2 = dlnorm(b[0], mu0, s0,0);
  F1 = 0.0;
  F2 = plnorm(b[0], mu0, s0,1,0);
  d1 += w[0] * pow(f2, 2.0) * x2/F2;
  d2 += w[0]*pow(f2, 2.0)*(pow(x2/s0,2.0)-1.0)/F2;
  for(i=1; i<n-1; i++){
    x1 = log(a[i]) - mu0;
    x2 = log(b[i]) - mu0;
    f1 = dlnorm(a[i], mu0, s0,0);
    f2 = dlnorm(b[i], mu0, s0,0);
    F1 = plnorm(a[i], mu0, s0,1,0);
    F2 = plnorm(b[i], mu0, s0,1,0);
    d1 += w[i] * (pow(f2, 2.0)*x2-pow(f1,2.0)*x1)/(F2-F1);
    d2 += w[i]*(pow(f2,2.0)*(pow(x2/s0,2.0)-1.0)-pow(f1,2.0)*(pow(x1/s0,2.0)-1.0))/(F2-F1);
  }
  x1 = log(a[n-1]) - mu0;
  x2 = 0.0;
  f1 = dlnorm(a[n-1], mu0, s0,0);
  f2 = 0.0;
  F1 = plnorm(a[n-1], mu0, s0,1,0);
  F2 = 1.0;
  d1 += - w[n-1] * pow(f1,2.0)*x1/(F2-F1);
  d2 += - w[n-1]*pow(f1,2.0)*(pow(x1/s0, 2.0)-1.0)/(F2-F1);

  F1 = pow(d1, 2.0) + pow(d2, 2.0);

  return F1;
}

void lnormBinMLE2(double *a, double *b, double *w,
		  int *size, double *mu, double *s)
{
  int i, j, iter = 1000, n=size[0];
  double mu0, s0, mu1, s1, step1, step2;
  double D0, Dmin;

  //range of searching: 
  // mu: 0.5mu0 to 1.5mu
  // sigma: 0.5s0 to 1.5s0
  step1 = mu[0]*2.0 / iter;
  step2 = s[0]*2.0 / iter;

  mu1 = mu[0];
  s1 = s[0];
  Dmin = lnormDist(a,b,w,n,mu[0],s[0]);
  Rprintf("Distance= %10.3f.\n", Dmin);

  s0 = 0.05 * s[0];
  for(i=0; i<iter; i++){
    s0 += step2; 
    mu0 = 0.05 * mu[0];
    for(j=0; j<iter; j++){
      mu0 += step1;
      D0 = lnormDist(a,b,w,n,mu0,s0);
      if(D0 < Dmin){
	mu1 = mu0;
	s1 = s0;
	Dmin = D0;
      }
    }
  }
  Rprintf("Distance= %10.3f.\n", Dmin);
  mu[0] = mu1; 
  s[0] = s1;
}

void slr(double *y, double *x, int *n,
	       double *a, double *b)
{
  int i;
  double mux=0.0, muy=0.0,ssxx=0.0, ssxy=0.0;
  
  for(i=0; i<n[0]; i++){
    mux += x[i];
    muy += y[i];
  }
  mux /= 1.0 * n[0];
  muy /= 1.0 * n[0];
  
  for(i=0; i<n[0]; i++){
    ssxx += pow(x[i]-mux, 2.0);
    ssxy += (x[i]-mux) * (y[i]-muy);
  }
  b[0] = ssxy/ssxx;
  a[0] = muy - b[0] * mux;
}

static double llkbinlnorm(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  
  double *tmp= (double*)ex, res=0.0;
  int i,n = (int)tmp[0]; //first element is the length of x;
  double mu0 = pars[0], s0= pars[1]; 
  double a[n], b[n], w[n], F1, F2, penalty=0.0;
  
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    a[i] = tmp[i+1]; 
    b[i] = tmp[i+n+1]; 
    w[i] = tmp[i+2*n+1]; 
  }
  
  for(i=0;i<n;i++) {
    F1 = plnorm(a[i], mu0, s0,1,0);
    F2 = plnorm(b[i], mu0, s0,1,0);
    if(w[i] > 0){
      if(F2>F1){
	res += w[i] * log(F2 - F1);
      }else{
	penalty += w[i];
      }
    }
  }
  
  return(-res/(1.0 - penalty/n));
}

void lnormBinMLE3(double *a, double *b, double *w,
		  int *size, double *mu, double *s)
{
  //maximize llk directly using Nelder'sMead simplex method
  int i,nx=size[0],npar=2;
  double dpar[npar],opar[npar]; 
  dpar[0] = mu[0]; dpar[1] = s[0]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[3*nx+1];
  yaux[0] = nx; //sample size
  for(i=0;i<nx;i++){
    yaux[i+1] = a[i];
    yaux[i+nx+1] = b[i];
    yaux[i+2*nx+1] = w[i];
  }
  nmmin(npar,dpar,opar,&val,llkbinlnorm,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  //Rprintf("mu = %10.3f.\n", opar[0]);
  //Rprintf("sigma = %10.3f.\n", opar[1]);

  mu[0] = opar[0]; s[0] = opar[1];
}


static double llkmixlnorm(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  
  double *tmp= (double*)ex, res=0.0;
  int i,n = (int)tmp[0];
  double n0 = tmp[1]; 
  //first element is the length of x;
  double mu0 = pars[0], s0= pars[1],xlmt=tmp[2]; 
  double x[n], w[n], fx, penalty=0.0;
  
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    x[i] = tmp[i+3]; 
    w[i] = tmp[i+n+3]; 
  }
  
  res = plnorm(xlmt, mu0, s0,1,1)*n0; //log(F(xmin))
  for(i=0;i<n;i++) {
    fx = dlnorm(x[i], mu0, s0, 0);
    if(fx > 0){
      res += w[i] * log(fx);
    }else{
      penalty += w[i];
    }
  }
  
  return(-res/(1.0 - penalty/n));
}

void lnormMix1(double *x, double *w, double *xlmt,
	       int *size, double *pars)
{
  //maximize llk directly using Nelder'sMead simplex method
  int i,nx=size[0],n0=size[1],npar=2;
  double dpar[npar],opar[npar]; 
  dpar[0] = pars[0]; dpar[1] = pars[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[2*nx+3];
  yaux[0] = nx; //sample size 1
  yaux[1] = n0; //sample size 2
  yaux[2] = xlmt[0]; //sample size 2
  for(i=0;i<nx;i++){
    yaux[i+3] = x[i];
    yaux[i+nx+3] = w[i];
  }
  nmmin(npar,dpar,opar,&val,llkmixlnorm,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  //Rprintf("mu = %10.3f.\n", opar[0]);
  //Rprintf("sigma = %10.3f.\n", opar[1]);
  
  pars[0] = opar[0]; pars[1] = opar[1];
}

static double llkmixnorm(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  
  double *tmp= (double*)ex, res=0.0;
  int i,n = (int)tmp[0];
  double n0 = tmp[1]; 
  //first element is the length of x;
  double mu0 = pars[0], s0= pars[1],xlmt=tmp[2]; 
  double x[n], w[n], fx, penalty=0.0;
  
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    x[i] = tmp[i+3]; 
    w[i] = tmp[i+n+3]; 
  }
  
  res = pnorm(xlmt, mu0, s0,1,1)*n0; //log(F(xmin))
  for(i=0;i<n;i++) {
    fx = dnorm(x[i], mu0, s0, 0);
    if(fx > 0){
      res += w[i] * log(fx);
    }else{
      penalty += w[i];
    }
  }
  
  return(-res/(1.0 - penalty/n));
}

void mleNorm0(double x[], double w[], int nx,
	      double xlmt, double wt0,
	      double mle[]){
  int i,npar=2;
  double dpar[npar],opar[npar]; 
  dpar[0] = mle[0]; dpar[1] = mle[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[2*nx+3];
  yaux[0] = nx; //sample size 1
  yaux[1] = wt0; //sample size 2
  yaux[2] = xlmt; 
  for(i=0;i<nx;i++){
    yaux[i+3] = x[i];
    yaux[i+nx+3] = w[i];
  }
  nmmin(npar,dpar,opar,&val,llkmixnorm,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);

  //Rprintf("subroutine: Mean=%10.3f, Std= %10.3f.\n", 
  //opar[0], opar[1]);
  mle[0] = opar[0]; mle[1] = opar[1];
}

void mleNorm1(double x[],double w[],int n,double mle[]){
  int i;
  double mu,s,xsum, wsum=0.0;
  xsum = 0.0; wsum = 0.0;
  for(i=0;i<n;i++){
    xsum += x[i] * w[i];
    wsum += w[i];
  }
  mu = xsum/wsum;
  xsum = 0.0; 
  for(i=0;i<n;i++){
    xsum += pow(x[i]-mu, 2.0) * w[i];
  }
  s = sqrt(xsum/(wsum-1.0));
  mle[0] = mu; mle[1] = s;
}

void lnormMixK(double *x, double *w, int *size, double *w0,
	       double *ps, double *mus, double *sigs)
{
  int i,j,k,ncomp=size[1],n=size[0];
  double n0=w0[0], xlmt=w0[1];
  double wts[n][ncomp], mu[ncomp], s[ncomp], p[ncomp];
  double llk1=-4.0e-9, llk2=-2.0e-9, tol=1.0e-2,delta;
  double tmp[ncomp],ttmp;
  double mle[2];
  double w2[n];

  for(i=0; i<ncomp;i++){
    mu[i] = mus[i];
    s[i] = sigs[i];
    p[i] = ps[i];
  }
  // iterate until converge: llk
  delta = 1.0;
  j = 0;
  while(delta > tol && j < 100){
    // compute weight of x=0 ====================
    p[0] = n0;
    if(ncomp > 1){
      for(i=1; i<ncomp;i++){
	p[i] = 0.0; // all 0's to comp#1
      }
    }
      
    /*ttmp = 0.0;
    for(i=0; i<ncomp;i++){
      tmp[i] = pnorm(xlmt, mu[i], s[i], 1, 0);
      ttmp += tmp[i];
      p[i] = 0.0; //to compute mixing coefficients
    }
    if(ttmp > 0){
      for(i=0; i<ncomp;i++){
	wt0[i] = tmp[i]/ttmp*n0;
	p[i] = wt0[i];
      }
      }*/
    // compute weight of x>0 ====================
    for(k=0; k<n;k++){
      ttmp = 0.0;
      for(i=0; i<ncomp;i++){
	tmp[i] = dnorm(x[k], mu[i], s[i], 0);
	ttmp += tmp[i];
      }
      if(ttmp > 0){
	for(i=0; i<ncomp;i++){
	  wts[k][i] = tmp[i]/ttmp*w[k];
	  p[i] += wts[k][i];
	}
      }
    }
    // compute mixing coefficients
    ttmp = 0.0;
    for(i=0; i<ncomp;i++){
      ttmp += p[i];
    }
    if(ttmp > 0){
      for(i=0; i<ncomp;i++){
	p[i] /= ttmp;
      }
    }
    // update the parameter estimates of mu&s
    for(i=0; i<ncomp;i++){
      for(k=0; k<n;k++){
	w2[k] = wts[k][i];
      }
      //Rprintf("Before: K=%d, Mean=%10.3f, Std= %10.3f.\n", 
      //i, mu[i], s[i]);
      mle[0] = mu[i];
      mle[1] = s[i];
      if(i==0){
	mleNorm0(x,w2,n,xlmt,n0,mle);
      }else{
	mleNorm1(x,w2,n,mle);
      }
      mu[i] = mle[0];
      s[i] = mle[1];
      //Rprintf("After: K=%d, Mean=%10.3f, Std= %10.3f.\n",
      // i, mu[i], s[i]);
    }

    // renew the log-likelihood ==================
    llk1 = llk2;
    llk2 = 0.0;
    ttmp = 0.0;
    for(i=0; i<ncomp;i++){
      ttmp += p[i] * pnorm(xlmt, mu[i], s[i], 1, 0);
    }
    if(ttmp > 0){
      llk2 += n0 * log(ttmp);
    }
    for(k=0; k<n;k++){
      ttmp = 0.0;
      for(i=0; i<ncomp;i++){
	ttmp += p[i] * dnorm(x[k], mu[i], s[i], 0);
      }
      if(ttmp > 0){
	llk2 += w[k] * log(ttmp);
      }
    }
    //update variables for convergence check
    delta = fabs(llk2 - llk1);
    j++;
  }
  for(i=0; i<ncomp;i++){
    ps[i] = p[i];
    mus[i] = mu[i];
    sigs[i] = s[i];
  }
  w0[0] = llk2;
  Rprintf("Iteration=%d, delta= %10.3f.\n", j, delta);

}

void lseNorm(double x[], double w[], int n, double mle[]){
  int i;
  double mi[n], qi[n];
  double y,xmu,ymu,ssxx,ssxy;

  y = 0.0;
  for(i=0;i<n;i++){
    mi[i] = x[i];
    y += w[i];
    qi[i] = y;
  }
  y += 0.5; 
  for(i=0;i<n;i++){
    qi[i] = qnorm(qi[i]/y,0.0, 1.0, 1,0);
  }
  //compute LSE
  xmu = 0.0; ymu = 0.0;
  for(i=0;i<n;i++){
    xmu += qi[i];
    ymu += mi[i];
  }
  xmu = xmu/y;
  ymu = ymu/y;
  ssxx = 0.0; ssxy = 0.0;
  for(i=0;i<n;i++){
    ssxx += (qi[i]-xmu)*(mi[i]-ymu)*w[i];
    ssxy += pow(qi[i]-xmu, 2.0)*w[i];
  }

  mle[1] = ssxy/ssxx;
  mle[0] = ymu - xmu*mle[1];
}

void lnormLSEK(double *x, double *w, int *size, 
	       double *ps, double *mus, double *sigs)
{
  int i,j,k,ncomp=size[1],n=size[0];
  double wts[n][ncomp], mu[ncomp], s[ncomp];
  double p0[ncomp], p[ncomp];
  double tol=1.0e-9,delta;
  double tmp[ncomp],ttmp;
  double mle[2];
  double w2[n];

  //stopping criteria 1) sum(abs(para-change))
  delta = 0.0;
  for(i=0; i<ncomp;i++){
    mu[i] = mus[i];
    s[i] = sigs[i];
    p[i] = ps[i];
    p0[i] = ps[i];
    delta += fabs(mu[i]);
    delta += fabs(s[i]);
    delta += fabs(p[i]);
  }
  // iterate until converge: absolute changes of LSEs
  j = 0;
  while(delta > tol && j < 500){
    delta = 0.0; //initialize delta
    // compute weight of x ====================
    for(k=0; k<n;k++){
      ttmp = 0.0;
      for(i=0; i<ncomp;i++){
	tmp[i] = dnorm(x[k], mu[i], s[i], 0);
	ttmp += tmp[i];
      }
      for(i=0; i<ncomp;i++){
	wts[k][i] = tmp[i]/ttmp*w[k];
	p[i] += wts[k][i];
      }
    }
    // compute mixing coefficients
    ttmp = 0.0;
    for(i=0; i<ncomp;i++){
      ttmp += p[i];
    }
    for(i=0; i<ncomp;i++){
      p[i] /= ttmp;
      delta += fabs(p[i]-p0[i]);
      p0[i] = p[i];
    }
    // update the parameter estimates of mu&s
    for(i=0; i<ncomp;i++){
      for(k=0; k<n;k++){
	w2[k] = wts[k][i];
      }
      //Rprintf("Before: K=%d, Mean=%10.3f, Std= %10.3f.\n", 
      //i, mu[i], s[i]);
      mle[0] = mu[i];
      mle[1] = s[i];
      lseNorm(x,w2,n,mle); //compute LSE
      delta += fabs(mu[i]-mle[0]);
      mu[i] = mle[0];
      delta += fabs(s[i]-mle[1]);
      s[i] = mle[1];
      //Rprintf("After: K=%d, Mean=%10.3f, Std= %10.3f.\n",
      // i, mu[i], s[i]);
    }

    //update variables for convergence check
    j++;
  }
  for(i=0; i<ncomp;i++){
    ps[i] = p[i];
    mus[i] = mu[i];
    sigs[i] = s[i];
  }
  //Rprintf("Iteration=%d, delta= %10.3f.\n", j, delta);
}


static double lsmixlnorm(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  double *tmp= (double*)ex, res=0.0;
  int i,j,n = (int)tmp[0];
  double mu,s,psum;
  
  for(i=0;i<n;i++) {
    psum = 0.0;
    for(j=0;j<npar;j++){
      mu = pars[j];
      s = pars[j+npar];
      if(s>0.0){
	psum += pnorm(tmp[i+1], mu, s,1,0);
      }else{
	psum += 1.0;
      }
      res += fabs(tmp[i+1+npar]-psum);
    }
  }
    
  return(res);
}

void lnormMixNM(double *x, double *F, int *size, 
		int *ncomp, double *mus, double *sigs)
{
  //minimize LS directly using Nelder'sMead simplex method
  int i,nx=size[0],m=ncomp[0], npar=2*m;
  double dpar[npar],opar[npar]; 
  //initial values
  for(i=0;i<m;i++){
    dpar[i] = mus[i];
    dpar[i+m] = sigs[i];
    opar[i] = mus[i];
    opar[i+m] = sigs[i];
  }

  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[2*nx+1];
  yaux[0] = nx; //sample size 1
  for(i=0;i<nx;i++){
    yaux[i+1] = x[i];
    yaux[i+nx+1] = F[i];
  }
  nmmin(npar,dpar,opar,&val,lsmixlnorm,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  
  for(i=0;i<m;i++){
    mus[i] = opar[i];
    sigs[i] = opar[i+m];
  }
}
