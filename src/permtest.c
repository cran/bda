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
