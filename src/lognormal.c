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


/*
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation (a copy of the GNU
 *  General Public License is available at
 *  http://www.r-project.org/Licenses/
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/*
 * Grid.Binning is designed to redistributed the weight along a
 * defined grid.
 *
 *  ngrid is the number of grid points, which equals to the bin number
 *  plus one.
 *
 *  Last updated: Nov 2, 2012
 */
void GridBinning(double *x, double *w, int *nx,
		 double *xlo, double *bw, int *ngrid,
		 int *truncate, int *linbin, double *gcounts)
{
  int i, li, m=ngrid[0], n=nx[0];
  double binwidth = bw[0], lxi, rem, a=xlo[0];
  
  for(i=0; i<m; i++) gcounts[i] = 0.0;
  
  for(i=0; i<n; i++){
    lxi = (x[i] - a)/binwidth;
    li = (int) lxi;
    if(linbin[0] == 1)
      rem = lxi - li;
    else
      rem = 0.0;

    if((li>0) && (li<m-1)){
      gcounts[li] += (1.0-rem) * w[i];
      gcounts[li+1] += rem * w[i];
    }
    
    if((li<=0)&&(truncate==0))
      gcounts[0] += w[i];
    if((li>=m-1)&&(truncate==0)&&(linbin[0] == 1))
      gcounts[m-1] += w[i];
    if((li>=m-1)&&(truncate==0)&&(linbin[0] == 0))
      gcounts[m-2] += w[i];
  }
}

static double rcllkweibull(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{

  double *tmp= (double*)ex, res=0.0;
  int i,n = (int)tmp[0]; //first element is the length of x;
  double kappa = pars[0], lambda= pars[1]; 
  double x[n], w[n];
  
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    x[i] = tmp[i+1]; 
    w[i] = tmp[i+n+1]; 
  }
  
  for(i=0;i<n;i++) {
    res += w[i]*(log(kappa) + (kappa-1.0)*log(x[i]) - kappa*log(lambda))
      -  pow(x[i]/lambda, kappa);
  }
  
  return(-res);
}

void RcMleWeibull(double *x,double *w,int *size,double *pars)
{
  int i,nx=size[0],npar=2;
  double dpar[npar],opar[npar]; 
  dpar[0] = pars[0]; dpar[1] = pars[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[2*nx+1];
  yaux[0] = nx; //sample size
  for(i=0;i<nx;i++){
    yaux[i+1] = x[i];
    yaux[i+nx+1] = w[i];
  }
  nmmin(npar,dpar,opar,&val,rcllkweibull,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  pars[0] = opar[0]; pars[1] = opar[1];
}

/*  
 * product limit estimate for data with right censoring
 * Call: rcple(x,y,n[0],...);
 *
 */
void myrcple(double x[], double w[], int n, double y[], double h[], int m) 
{
  int i,j;
  double xprod=1.0;
  for(i=0;i<m;i++){
    h[i] = 1.0; 
  }
  i = 0; j = 0;
  while(j < m){
    if(y[j] <= x[i]){
      h[j] = xprod;
      j++;
    }else{
      i++;
      if(i < n)
	xprod *= pow((n-i)/(n-i+1.0), 1.0-w[i]);
      else xprod = 0.0;
    }
  }
}



void wkdemae(double *x,double *w,int *size,double *y,int *ny)
{
  int i, j, n=size[0], m=ny[0];
  double lambda=0.0,delta=0.0;
  double Hx[m], HX[n], hx[m], x0;
  for(i=0; i<n;i++){
    lambda += x[i];
    delta  += w[i];
  }
  lambda /= delta; //mle of lambda
  myrcple(x,w,n,y,Hx,m);
  myrcple(x,w,n,x,HX,n);
  double t1,t2;
  t1 = 0.7644174 * lambda * pow(n,-.2);
  t2 = 0.2/lambda;

  for(i=0; i<m;i++){
    hx[i] = t1 * exp(t2) * pow(Hx[i],-.2);
  }

  for(i=0; i<m;i++){
    x0 = y[i]; y[i] = 0.0; // reuse y[]
    for(j=0; j<n; j++){
      t1 = (x0 - x[j])/hx[i];
      y[i] += w[j]/(HX[j]*hx[i])*dnorm(t1,0.,1.0,0);
    }
  }
  
  for(i=0; i<m;i++){
    y[i] /= n;
  }

}


void BDMLE(double *f,double *a,double *b,int *nbin,
	   double *pars, int *npars, int *dist)
{
  int i,nx=nbin[0],npar=2; //2-parameter distribution only
  double dpar[npar],opar[npar]; 
  dpar[0] = pars[0]; dpar[1] = pars[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[3*nx+1];
  yaux[0] = nx; //sample size
  for(i=0;i<nx;i++){
    yaux[i+1] = f[i];
    yaux[i+nx+1] = a[i];
    yaux[i+2*nx+1] = b[i];
  }

  if(dist[0]==0) npars[0] = 2;  //reserved for later

  nmmin(npar,dpar,opar,&val,rcllkweibull,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  pars[0] = opar[0]; pars[1] = opar[1];
}

double bllkWeibull(double x[], double counts[], double kappa, 
		   double lambda, double alpha, int n, int nu)
{
  int i;
  double res=0.0, tmp=0.0;
  tmp = counts[0] * pow(1.0 - exp(-pow(x[0]/lambda,kappa)), alpha);
  if(tmp>0){
    res = log(tmp);
  }

  for(i=1;i<n; i++){
    tmp = counts[i] * (pow(1.0-exp(-pow(x[i]/lambda,kappa)),alpha) - 
		       pow(1.0-exp(-pow(x[i-1]/lambda,kappa)), alpha));
    if(tmp>0){
      res += log(tmp);
    }
  }
  tmp = nu * (1.0 - pow(1.0 - exp(-pow(x[0]/lambda,kappa)), alpha));
  if(tmp>0){
    res += log(tmp);
  }

  return(res);
}

double bllkDagum(double x[], double counts[], double kappa, 
		 double lambda, double alpha, int n, int nu)
{
  int i;
  double res=0.0, tmp=0.0;
  tmp = counts[0] * pow(1.0 + pow(x[0]/lambda,-kappa), -alpha);
  if(tmp>0.){
    res = log(tmp);
  }

  for(i=1;i<n; i++){
    tmp = counts[i] * (pow(1.0 + pow(x[i]/lambda,-kappa), -alpha) - 
		       pow(1.0 + pow(x[i-1]/lambda,-kappa), -alpha));
    if(tmp>0.){
      res += log(tmp);
    }
  }
  tmp = nu * (1.0 - pow(1.0 + pow(x[n-1]/lambda,-kappa), -alpha));
  if(tmp>0.){
    res += log(tmp);
  }
  return(res);
}

void bdrWeibull(double F[], double X[], double counts[], int n, int nu, double pars[]) 
{
  int i,j;
  double y[n], x[n],xbar,ybar,ssxy, ssxx, tmp, llk=0.0,alpha;
  alpha = pars[2];  //passed to here: cannot be zero or negative.
  xbar = 0.0;
  ybar = 0.0;
  for(i=0; i<n; i++){
    y[i] = log(-log(1.0-exp(log(F[i])/alpha)));
    x[i] = log(X[i]);
    xbar += x[i];
    ybar += y[i];
  }
  xbar /= n;
  ybar /= n;

  ssxy = 0.0; ssxx = 0.0;
  for(i=0; i<n; i++){
    tmp = x[i]-xbar;
    ssxy += tmp * (y[i]-ybar);
    ssxx += tmp * tmp;
  }
  tmp = ssxy/ssxx;
  pars[0] = tmp;
  pars[1] =  exp(xbar - ybar/tmp);

  llk = bllkWeibull(X, counts, pars[0], pars[1], alpha, n,nu);
  pars[2] = llk;

  double lstep, kstep,lambda0,kappa0;
  lstep = 0.01 * pars[1];
  kstep = 0.01 * pars[0];
  lambda0 = 0.8 * pars[1];
  kappa0 = 0.8 * pars[0];
  for(i=0; i<50; i++){
    for(j=0;j<50;j++){
      tmp = bllkWeibull(X, counts, kappa0, lambda0,alpha, n,nu);
      if(tmp > pars[2]){
	pars[0] = kappa0;
	pars[1] = lambda0;
	pars[2] = tmp;
      }
      kappa0 += kstep;
    }
    lambda0 += lstep;
  }
}


void bdrDagum(double F[], double X[], double counts[], int n, int nu, double pars[]) 
{
  int i,j;
  double a,b;
  double y[n], x[n],xbar,ybar,ssxy, ssxx, tmp, llk=0.0,alpha;
  alpha = pars[2];  //passed to here: cannot be zero or negative.
  xbar = 0.0;
  ybar = 0.0;
  for(i=0; i<n; i++){
    y[i] = log(exp(-log(F[i])/alpha)-1.0);
    x[i] = log(X[i]);
    xbar += x[i];
    ybar += y[i];
  }
  xbar /= n;
  ybar /= n;

  ssxy = 0.0; ssxx = 0.0;
  for(i=0; i<n; i++){
    tmp = x[i]-xbar;
    ssxy += tmp * (y[i]-ybar);
    ssxx += tmp * tmp;
  }
  b = ssxy/ssxx; a = ybar - b * xbar;
  pars[0] = -b;  //parameter: a
  pars[1] =  exp(-a/b); //parameter: b

  llk = bllkDagum(X, counts, pars[0], pars[1],alpha,n,nu);
  pars[2] = llk;

  double lstep, kstep,lambda0,kappa0;
  lstep = 0.01 * pars[1];
  kstep = 0.01 * pars[0];
  lambda0 = 0.8 * pars[1];
  kappa0 = 0.8 * pars[0];
  for(i=0; i<40; i++){
    for(j=0;j<40;j++){
      tmp = bllkWeibull(X, counts, kappa0, lambda0,alpha,n,nu);
      if(tmp > pars[2]){
	pars[0] = kappa0;
	pars[1] = lambda0;
	pars[2] = tmp;
      }
      kappa0 += kstep;
    }
    lambda0 += lstep;
  }
}

void bdregmle(double *F, double *x, double *counts,
	      int *nusize, int *size, int *dist, double *pars)
{
  int i,n=size[0], nu = nusize[0];
  double llk,lambda=0.0, kappa=0.0,alpha=0.0,tmp;

  switch(dist[0]){
  case 1: //EWD
    alpha = 1.;
    pars[2] = alpha;
    bdrWeibull(F, x, counts, n, nu, pars);
    llk = pars[2];
    tmp = 0.5;
    for(i=0; i<40;i++){
      tmp += 0.02;
      pars[2] = tmp;
      bdrWeibull(F, x, counts, n, nu, pars);
      if(pars[2] > llk && R_FINITE(pars[2])){
	llk = pars[2];
	alpha = tmp;
	kappa = pars[0];
	lambda = pars[1];
      }
    }
    pars[0] = kappa;
    pars[1] = lambda;
    pars[2] = alpha;
    break;
  case 2: //Dagum
    alpha = 0.0001;
    pars[2] = alpha;
    bdrDagum(F, x, counts, n, nu, pars);
    llk = pars[2];
    tmp = alpha;
    for(i=0; i<1000;i++){
      if(tmp < 1.5){
	tmp += 0.002;
      }else{
	tmp += 0.1;
      }
      pars[2] = tmp;
      bdrDagum(F, x, counts, n, nu, pars);
      if(pars[2] > llk && R_FINITE(pars[2])){
	llk = pars[2];
	alpha = tmp;
	kappa = pars[0];
	lambda = pars[1];
      }
    }
    pars[0] = kappa;
    pars[1] = lambda;
    pars[2] = alpha;
    break;
  default:
    pars[2] = 1.0;
    bdrWeibull(F, x, counts, n, nu, pars);
  }
}

double qGldFmkl(double u, double lambdas[])
{
  double l1=lambdas[0],l2=lambdas[1],l3=lambdas[2],l4=lambdas[3];
  return(l1+((pow(u,l3)-1.)/l3-(pow(1.-u,l4)-1.)/l4)/l2);
}

double gRootGldFmkl(double u, double lambdas[], double q)
{
  return(qGldFmkl(u, lambdas) - q);
}

void rootGldFmklBisection(double *q, double *lambdas)
{
  int  iter=100, ctr=1;
  double delta=0.5, tol=1.e-8;
  double l1=0.0, l2=1.0, r=0.0, f1, f2, f3;
  if(!isfinite(q[0])){
    if(q[0]>0.)
      r = 0.0;
    else
      r = 1.0;
  }else{
    f1 = gRootGldFmkl(l1,lambdas,q[0]);
    f2 = gRootGldFmkl(l2,lambdas,q[0]);
    f3 = f2; 
    if(f1 == 0.0)
      r = l1;
    else if(f2 == 0.0)
      r = l2;
    else if(f1 * f2 > 0.0){
      if(f1 > 0.0)
	r = 0;
      else
	r = 1;
    }else{
      while(ctr <= iter && fabs(delta) > tol){
	r = 0.5 * (l1 + l2);
	f3 = gRootGldFmkl(r,lambdas,q[0]);
	delta = f2 - f3;
	if(f3 == 0) break;
	if(f1 * f3 < 0.0){
	  l2 = r;
	  f2 = f3;
	}
	else{
	  l1 = r;
	  f1 = f3;
	}
	ctr++;
      }
    }
  }
  q[0] = r;
}


void bootsd(int *size, double *x, double *y, double *sx, double *sy,
	    int *iter, double *sig, double *rho)
{
  int  n=size[0],i,j;
  double xr, yr, dyx, d[n],dbar,dsum,xbar,ybar,xysum,xsum,ysum;

  for(i=0; i<iter[0]; i++){
    dsum = 0.0;
    xsum = 0.0;
    ysum = 0.0;
    xbar = 0.0;
    ybar = 0.0;
    xysum = 0.0;
    for(j=0; j<n; j++){
      xr = x[j] + rnorm(0.0, sx[j]);
      yr = y[j] + rnorm(0.0, sy[j]);
      dyx = yr - xr;
      d[j] = dyx;
      dsum += dyx;
      xbar += xr;
      ybar += yr;
      xsum += pow(xr,2.0);
      ysum += pow(yr,2.0);
      xysum += xr * yr;
    }
    dbar = dsum/n;
    xbar = xbar/n;
    ybar = ybar/n;
    dsum = 0.0;
    for(j=0; j<n; j++){
      dsum += pow(d[j]-dbar,2.0);
    }
    rho[i] = (xysum - n*xbar*ybar)/sqrt((xsum-n*xbar*xbar)*(ysum-n*ybar*ybar));
    sig[i] = sqrt(dsum/(n-1.0));
  }
}

void fitdpro1(double *ll, double *ul, int *n, double *mu, double *s)
{
  int  i,j,k;
  double dmu,ds,mu0,s0,mu1,s1,maxllk,llk,Fx1,Fx0,dFx;
  dmu = mu[0] * 0.01;
  ds = s[0] * 0.063;
  mu0 = 0.8 * mu[0];
  s0 = 0.9 * s[0];
  mu1 = mu0; s1 = s0;

  maxllk = -999.99;
  for(i = 0; i < 50; i++){
    for(j = 0; j < 50; j++){
      llk = 0.0;
      for(k = 0; k < n[0]; k++){
	if(fabs(ll[k])>100){
	  Fx0 = 0.0;
	}else{
	  Fx0 = pnorm(ll[k], mu0, s0, 1, 0);
	}
	
	if(fabs(ll[k])>100){
	  Fx1 = 1.0;
	}else{
	  Fx1 = pnorm(ul[k], mu0, s0, 1, 0);
	}
	
	dFx = fabs(Fx1 - Fx0);
	if(dFx > 0.00000001){
	  llk += log(dFx);
	}
      }
      if(llk > maxllk){
	maxllk = llk;
	mu1 = mu0;
	s1 = s0;
      }
      s0 += ds;
    }
    mu0 += dmu;
  }
  //printf("Mean=%10.2f,SD=%10.2f\n",mu1,s1);
  mu[0] = mu1;
  s[0] = s1;
}

void fitdpro2(double *ll, double *ul, int *n2,
	      double *x, int *n1,
	      double *mu, double *s)
{
  int  i,j,k;
  double dmu,ds,mu0,s0,mu1,s1,maxllk,llk,Fx1,Fx0,dFx;
  dmu = mu[0] * 0.005;
  ds = s[0] * 0.033;
  mu0 = 0.8 * mu[0];
  s0 = 0.9 * s[0];
  mu1 = mu0; s1 = s0;

  maxllk = -999.99;
  for(i = 0; i < 100; i++){
    for(j = 0; j < 100; j++){
      llk = 0.0;
      for(k = 0; k < n2[0]; k++){
	if(fabs(ll[k])>100){
	  Fx0 = 0.0;
	}else{
	  Fx0 = pnorm(ll[k], mu0, s0, 1, 0);
	}
	
	if(fabs(ll[k])>100){
	  Fx1 = 1.0;
	}else{
	  Fx1 = pnorm(ul[k], mu0, s0, 1, 0);
	}

	dFx = fabs(Fx1 - Fx0);
	if(dFx > 0.00000001){
	  llk += log(dFx);
	}
      }
      for(k = 0; k < n1[0]; k++){
	Fx0 = dnorm(x[k], mu0, s0, 0);
	if(Fx0 > 0.00000001){
	  llk += log(Fx0);
	}
      }
      if(llk > maxllk){
	maxllk = llk;
	mu1 = mu0;
	s1 = s0;
      }
      s0 += ds;
    }
    mu0 += dmu;
  }
  mu[0] = mu1;
  s[0] = s1;
}


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


//to compute the D=max(Fx-Fy)
void compFnx(double X[], int n, double a, double b, int M, double Fx[])
{
  int i,li;
  double delta, lxi;
  for(i=0; i<M; i++){
    Fx[i] = 0.0;
  }
  delta = (b-a)/(M-1.0);
  for(i=0; i<n; i++){
    lxi = (X[i]-a)/delta;
    li = (int) floor(lxi);
    Fx[li] = Fx[li] + 1.0;
  }
  lxi = 0.0;
  for(i=0; i<M; i++){
    lxi = lxi + Fx[i]/n;
    Fx[i] =  lxi;
  }
}

double compD(int nx,double x[],int ny,double y[],double a,double b,int ngrid)
{
  int i;
  double Fx[ngrid], Fy[ngrid], xmax,xtmp;
  for(i=0; i<ngrid; i++){
    Fx[i] = 0.0;
    Fy[i] = 0.0;
  }
  compFnx(x,nx,a,b,ngrid,Fx);
  compFnx(y,ny,a,b,ngrid,Fy);
  xmax = fabs(Fx[0]-Fy[0]);
  for(i=1; i<ngrid; i++){
    xtmp = fabs(Fx[i] - Fy[i]);
    if(xtmp > xmax)
      xmax = xtmp;
  }
  return(xmax);
}


void permtest(double *x, int *nx, double *y, int *ny, double *a,
	      double *b, double *D, double *pv, int *iter)
{
  int i, j, k,l;
  int ngrid=1001;
  double xy1[nx[0]+ny[0]], xy2[nx[0]+ny[0]];
  double xmin=a[0], xmax=b[0],D1=0.0;

  GetRNGstate();
  
  //get min and max of positive values and define grid
  for(i=0; i<nx[0]; i++){
    xy1[i] = x[i];
  }
  for(i=0; i<ny[0]; i++){
    xy1[i+nx[0]] = y[i];
  }
  D[0] = compD(nx[0],x,ny[0],y,xmin,xmax,ngrid);
  pv[0] = 1.0;
  for(i=0; i<iter[0]; i++){
    for(j=0; j<nx[0]+ny[0]; j++){
      xy2[j] = xy1[j];
    }
    //resample
    k = nx[0] + ny[0];
    for(j=0; j<nx[0]; j++){
      l = (int) runif(0.0, k-1.0);
      x[j] = xy2[l];
      k--;
      xy2[l] = xy2[k];
    }
    for(j=0; j<ny[0]; j++){
      y[j] = xy2[j];
    }
    D1 = compD(nx[0],x,ny[0],y,xmin,xmax,ngrid);
    if(D[0] <= D1){
      pv[0] = pv[0] + 1.0;
    }
  }
  pv[0] = pv[0]/(iter[0] + 1.0);
  PutRNGstate();
 
}



void permtest2(double *D, int *M, int *nx, int *ny,
	       int *F, int *iter)
{
  int i, j, k, l, pv=1,n=nx[0]+ny[0];
  int x[M[0]],F2[n];
  double Dmax,D1,D2,Dif;

  GetRNGstate();

  
  for(i=0; i<iter[0]; i++){
    //initialize F2
    l = 0;
    for(j=0; j<M[0]; j++){
      for(k=0; k<F[j]; k++){
	F2[l] = j;
	l++;
      }
    }
    //resample
    for(j=0; j<M[0]; j++){
      x[j] = 0;
    }
    k = n;
    for(j=0; j<nx[0]; j++){
      l = (int) runif(0.0, k-1.0);
      //      printf("k=%d,l=%d,F2[l]=%d,x[before]=%d\n",k,l,F2[l],x[F2[l]]);
      x[F2[l]]++;
      //      printf("x[after]=%d\n",x[F2[l]]);
      k--;
      F2[l] = F2[k];
    }
    //Find max(|Fn(x)-Fn(y)|)
    k = 0; l = 0;
    Dmax = 0.0; D2 = 0.0; D1 = 0.0; Dif = 0.0;
    for(j=0; j<M[0]; j++){
      k += x[j];
      l += F[j] - x[j];
      //      printf("k=%d,l=%d,n.x=%d,n.y=%d,",k,l,nx[0],ny[0]);
      D1 = (double) k/nx[0];
      D2 = (double) l/ny[0];
      Dif = fabs(D1 - D2);
      if(Dif > Dmax) Dmax = Dif;
      //      printf("D1=%f,D2=%f,Diff=%f,\nDmax=%f,D0=%f\n",D1, D2, Dif,Dmax,D[0]);
    }
    if(D[0] <= Dmax) pv++;
    //    printf("counts=%d\n\n",pv);
  }
  D[0] = pv/(iter[0] + 1.0);
  PutRNGstate();
}



void permtest3(double *xy, int *nx, int *ny,
	       double *pv, int *iter)
{
  int i, j, k, l, n=nx[0]+ny[0],plt=0,pne=0,pgt=0;
  double D=pv[0],xysum=pv[1],D1=0.0,D2,xsum,xy0[n];
  D2 = fabs(D);
  
  GetRNGstate();
  
  for(i=0; i<iter[0]; i++){
    //initialize the pooled values
    for(j=0; j<n; j++){
      xy0[j] = xy[j];
    }
    
    //resampling
    k = n-1; xsum = 0.0;
    for(j=0; j<nx[0]; j++){
      l = (int) runif(0.0, k);
      xsum += xy0[l];
      xy0[l] = xy0[k];
      k--;
    }
    //Compute diff(mu1-mu2)
    D1 = xsum/nx[0] - (xysum-xsum)/ny[0];
    if(D > D1) plt++;
    if(D < D1) pgt++;
    if(D2 < fabs(D1)) pne++;
  }
  pv[0] = (double) plt/(iter[0] + 1.0);
  pv[1] = (double) pne/(iter[0] + 1.0);
  pv[2] = (double) pgt/(iter[0] + 1.0);
  PutRNGstate();
}


//recursive version to compute factorial of n := n!
int factorial(int n){
  return n>=1 ? n * factorial(n-1) : 1;
}

double g1(double p, int m1, int n11, double a[], double alpha){
  int i;
  double p1=0.0, p2=0.0;
  for(i=0;i<n11;i++)
    p1 += a[i] * pow(p, i);
  for(i=n11;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
    p2 += a[i] * pow(p, i);
  }
  return p2/p1 - 0.5 * alpha;
}

double g2(double p, int m1, int n11, double a[], double alpha){
  int i;
  double p1=0.0, p2=0.0;
  for(i=0;i<n11+1;i++){
    p1 += a[i] * pow(p, i);
    p2 += a[i] * pow(p, i);
  }
  for(i=n11+1;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
  }
  return p2/p1 - 0.5 * alpha;
}

double dg1(double p, int m1, int n11, double a[]){
  int i;
  double p1, p2=0.0, dp1=0.0, dp2=0.0;

  p1 = a[0];
  for(i=1;i<n11;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
  }
  for(i=n11;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
    p2 += a[i] * pow(p, i);
    dp2 += i * a[i] * pow(p, i-1);
  }
  return (dp2 * p1 - p2 * dp1)/(p1 * p1);
}

double dg2(double p, int m1, int n11, double a[]){
  int i;
  double p1, p2=0.0, dp1=0.0, dp2=0.0;

  p1 = a[0];
  p2 = a[0];
  for(i=1;i<n11+1;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
    dp2 += i * a[i] * pow(p, i-1);
  }
  for(i=n11+1;i<m1+1;i++){
    p1 += a[i] * pow(p, i);
    dp1 += i * a[i] * pow(p, i-1);
  }
  return (dp2 * p1 - p2 * dp1)/(p1 * p1);
}
/* use out to pass the initial value and the final estimate

   The Newton method was used, fast but not stable.  We change to the
   bisection method instead.  We take small value a = 1e-16 and b=1e6.
   If not in the range, simply use 0 or 100000+
 */

void orexactl(int *counts, double *alpha, double *out)
{
  int i,n11,n12,n21,n22, n1, n2, m1;
  n11 = counts[0];
  n12 = counts[1];
  n21 = counts[2];
  n22 = counts[3];
  n1 = n11+n12;
  n2 = n21+n22;
  m1 = n11+n21;
  double delta = 1.0, p0 = out[0],f0, fa, fb,pa,pb;
  double a[m1 + 1]; // to store the coefficients
  

  for(i=0;i < m1+1;i++)
    a[i] = choose(n1, i) * choose(n2, m1 - i);
  
  i = 0;
  pa = 1.e-16; pb = 1.e6;
  f0 = g1(p0, m1, n11, a, alpha[0]);
  if(f0>0.) pb = p0; else pa = p0;
  fa = g1(pa, m1, n11, a, alpha[0]);
  fb = g1(pb, m1, n11, a, alpha[0]);   

  while(i<10000 && fabs(delta) >0.00001){
    if(fa >=0.){
      p0 = pa; break; //exit
    }else if(fb<=0.){
      p0 = pb; break;
    }else{
      p0 = 0.5 * (pa + pb);      
      f0 = g1(p0, m1, n11, a, alpha[0]);
      if(f0>0.){
	pb = p0; 
	fb = g1(pb, m1, n11, a, alpha[0]);   
      }else{
	pa = p0;
	fa = g1(pa, m1, n11, a, alpha[0]);
      }
      i++;
    }
  }

  out[0] = p0;
}

void orexactu(int *counts, double *alpha, double *out)
{
  int i,n11,n12,n21,n22, n1, n2, m1;
  n11 = counts[0];
  n12 = counts[1];
  n21 = counts[2];
  n22 = counts[3];
  n1 = n11+n12;
  n2 = n21+n22;
  m1 = n11+n21;
  double delta = 1.0, p0 = out[0],f0, fa, fb,pa,pb;
  double a[m1 + 1]; // to store the coefficients
  

  for(i=0;i < m1+1;i++)
    a[i] = choose(n1, i) * choose(n2, m1 - i);
  
  i = 0;
  pa = 1.e-16; pb = 1.e6;
  f0 = g2(p0, m1, n11, a, alpha[0]);
  if(f0 < 0.) pb = p0; else pa = p0;
  fa = g2(pa, m1, n11, a, alpha[0]);
  fb = g2(pb, m1, n11, a, alpha[0]);   

  while(i<10000 && fabs(delta) >0.00001){
    if(fa <=0.){
      p0 = pa; break; //exit
    }else if(fb>=0.){
      p0 = pb; break;
    }else{
      p0 = 0.5 * (pa + pb);      
      f0 = g2(p0, m1, n11, a, alpha[0]);
      if(f0<0.){
	pb = p0; 
	fb = g2(pb, m1, n11, a, alpha[0]);   
      }else{
	pa = p0;
	fa = g2(pa, m1, n11, a, alpha[0]);
      }
      i++;
    }
  }
  //  p0 = g2(10., m1, n11, a, alpha[0]);

  out[0] = p0;
}

