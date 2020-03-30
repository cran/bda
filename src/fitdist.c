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
  double l1=0.0, l2=1.0, r, f1, f2, f3;
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

