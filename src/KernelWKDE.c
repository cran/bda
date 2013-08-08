#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include <Rinternals.h>
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

/*  binned data analysis*/

static double bdcdf(int npar, double *pars, void *ex)
/* loglikelihood function to be called by the Nelder-Mead simplex
method */
{

  double *tmp= (double*)ex, res=0.0;
  int i,n = (int)tmp[0]; //first element is the length of x;
  double kappa = pars[1], lambda= pars[0]; 
  double f[n], a[n], b[n], llk1,llk2;

  if(lambda > 0.0 && kappa > 0.0){
    res = 0.0;
    for(i=0;i<n;i++) {//restore auxiliary information from ex;
      f[i] = tmp[i+1]; 
      a[i] = tmp[i+n+1]; 
      b[i] = tmp[i+2*n+1]; 
    }
    llk1 = exp(-pow(a[0]/lambda, kappa));
    for(i=0;i<n;i++) {
      if(finite(b[i])) {
	llk2 = exp(-pow(b[i]/lambda, kappa));
      }else{ 
	llk2 = 0.0;
      }
      res += f[i]*log(llk1-llk2);
      llk1 = llk2;
    }
    res = -res;
  }else{
    res = 999999999999.99;
  }
  return(res);
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
