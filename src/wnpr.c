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

double Knorm(double z){
  return dnorm(z,0.0,1.0,0);
}

/*Instead of choosing the optimal bandwidth by minimizing the
  leave-one-out cross-validation score, we minimize the generalize
  cross-validation score by replacing Lii with nu/n=average of
  (Lii)=tr(L).
 */
double wnprgcv(double x[],double y[],double w[],int n,double h,double s){
  int i,j;
  double K[n], l[n], rhat[n];
  double nu=0.0,S1, S2, t1, t2, bsum,tmp,z;


  for(i=0; i<n; i++){
    // to compute S1, and S2
    S1 = 0.0;
    S2 = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - x[i];
      z = t1/h;
      tmp = 1. + pow(s/h,2.0)*(1.0-z*z);
      K[j] = w[j] * Knorm(z)*tmp;

      t2 = K[j] * t1;
      S1 += t2;
      S2 += t2 * t1;
    }
    // compute li: use li[i][j] to store b[i][j]
    bsum = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - x[i];
      l[j] = w[j] * K[j]*(S2-S1*t1); //bi(x)
      bsum += l[j];
    }
    rhat[i] = 0.0;
    for(j=0; j<n; j++){
      l[j] /= bsum;
      rhat[i] += l[j]*y[j];
    }
    nu += l[i];
  }
  
  bsum = 0.0; // lscv score
  for(j=0; j<n; j++){
    //generalized cross-validation (an approximation avoid Inf, NaN)
    t1 = (y[j]-rhat[j]);
    bsum += t1 * t1;
  }
  //  return bsum/n;  //cross-validation score
  return bsum/n/(1.0 - nu/n)/(1.0 - nu/n);
}

//the minimum may not exist. check the scatter plot to make sure hopt is valid
double hgcv(double x[],double y[],double w[],int n,double h, double s)
{  
  int i;
  double delta=0.03*h, h0, hopt, Rhmin=1.0e7, Rhat=0.0;
  h0 = 0.3 * h; hopt = h0;
  for(i=0; i<101; i++){
    Rhat = wnprgcv(x,y,w,n,h0,s);
    if(Rhat <= Rhmin && R_FINITE(Rhat)){
      hopt = h0;
      Rhmin = Rhat;
    }
    h0 += delta;
  }  
  return hopt;
}

void wnpreg(double xgrid[], int m, double x[], double y[],
	    double w[],int n,double h,double rhat[],double s)
{
  int i,j;
  double K[n], l[n];
  double t1,z, bsum,tmp;
  
  for(i=0; i<m; i++){
    bsum = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - xgrid[i];
      z = t1/h;
      tmp = 1. + pow(s/h,2.0)*(1.0-z*z);
      K[j] = w[j] * Knorm(z)*tmp;
      bsum += K[j];
    }
    rhat[i] = 0.0;
    for(j=0; j<n; j++){
      l[j] = K[j]/bsum;
      rhat[i] += l[j]*y[j];
    }
  }
}

void wlpreg(double xgrid[], int m, double x[], double y[],
	    double w[],int n,double h,double rhat[],double s)
{
  int i,j;
  double K[n], l[n];
  double S1, S2, t1, t2,z, bsum,tmp;
  
  for(i=0; i<m; i++){
    // to compute K, S1, and S2
    S1 = 0.0;
    S2 = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - xgrid[i];
      z = t1/h;
      tmp = 1. + pow(s/h,2.0)*(1.0-z*z);
      K[j] = w[j] * Knorm(z)*tmp;
      t2 = K[j] * t1;
      S1 += t2;
      S2 += t2 * t1;
    }
    // compute li: use li[i][j] to store b[i][j]
    bsum = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - xgrid[i];
      l[j] = w[j] * K[j]*(S2-S1*t1); //bi(x)
      bsum += l[j];
    }
    rhat[i] = 0.0;
    for(j=0; j<n; j++){
      l[j] /= bsum;
      rhat[i] += l[j]*y[j];
    }
  }
}

void wnpr(double *xgrid, int *ngrid, double *x, double *y, 
	  double *w, int *size, double *bw, double *sx)
{
  int i;
  double rhat[ngrid[0]],hopt=bw[0];

  for(i=0; i<ngrid[0]; i++){
    rhat[i] = 0.0;
  }

  hopt = hgcv(x,y,w,size[0],hopt,sx[0]);
  bw[0] = hopt;

  if(sx[0] <= 0.0){
    wnpreg(xgrid, ngrid[0],x,y,w,size[0],hopt,rhat,sx[0]);
  }else{
    wlpreg(xgrid, ngrid[0],x,y,w,size[0],hopt,rhat,sx[0]);
  }

  for(i=0; i<ngrid[0]; i++){
    xgrid[i] = rhat[i];
  }
}

void wdekde(double *x, double *w, int *n,
	    double *xgrid, int *ngrid,  
	    double *bw, double *sx)
{
  int i,j;
  double fx[ngrid[0]],xh=0.0,tmp=0.0;

  for(i=0; i<ngrid[0]; i++){
    fx[i] = 0.0;
  }

  for(i=0; i<ngrid[0]; i++){
    for(j=0; j<n[0]; j++){
      xh = (xgrid[i]-x[j])/bw[0];
      tmp = 1. + pow(sx[0]/bw[0],2.0)*(1.0-xh*xh);
      fx[i] += w[j] * dnorm(xh,0.0,1.0,0)*tmp;
    }
    fx[i] /= bw[0];
  }
  
  
  for(i=0; i<ngrid[0]; i++){
    xgrid[i] = fx[i];
  }
}
