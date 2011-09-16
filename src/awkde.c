#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"

// ISNA(x): true for R's NA only
// ISNAN(x): true for R's NA and IEEE NaN
// R_FINITE(x): false for Inf,-Inf,NA,NaN
// R_IsNaN(x): true for NaN but not NA

// Rprintf:  printing from a C routine compiled into R

void awcdf(double x[],int nx,double w[],double h,double hs[],
	   double y[],int ny,double Fy[])
{
  int i,j;
  for(i=0;i<ny;i++) Fy[i]=0.0;
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      Fy[j] += w[i]*pnorm(y[j]-x[i],0.,h*hs[i],1,0);
    }
  }
}

 

void awpdf(double x[], int nx, double w[], double h, double hs[],
	       double y[], int ny, double fy[], int rr)
{
  int i,j;
  double kappa = 0.0;
  //  double wsum = 1.,kappa = 0.0;
  for(i=0;i<ny;i++) fy[i]=0.0;
  //  for(i=0;i<nx;i++){
  //    if(x[i]<4.*h){wsum += w[i];}
  //  }
  for(i=0;i<ny;i++){
    for(j=0;j<nx;j++){
      if(rr!=1){
	fy[i] += w[j]*dnorm(y[i]-x[j],0.,h*hs[j],0);
      }else{
	if(x[j]<4.*h){
	  fy[i] += w[j]*(dnorm(y[i]-x[j],0.,h*hs[j],0) + 
			    dnorm(y[i]+x[j],0.,h*hs[j],0));
	}else{
	  fy[i] += w[j]*dnorm(y[i]-x[j],0.,h*hs[j],0);
	}
      }
    }
    //    fy[i] /=wsum;
    kappa += fy[i];
  }
  kappa *= fabs(y[ny-1]-y[0])/(ny-1.);
  for(i=0;i<ny;i++){
    fy[i] /= kappa;
  }

}
 
void UpdateBwfactor(double fx[], int n, double alpha,double hs[])
{
  double g=0.0;
  int i;
  for(i=0;i<n;i++){g += log(fx[i]);}
  g = exp(g/n);
  for(i=0;i<n;i++) {
    hs[i] = pow(fx[i]/g,-alpha);
  }
}


static double wise2(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  double res2=0.0, pi2 = M_2_SQRTPI/4., res=0.0, *tmp= (double*)ex;
  int i,j,n = (int)tmp[0]; //first element is the length of x;
  double x[n],hs[n],w[n],fx[n],h0;
  h0 = tmp[1]; //second element is h0
  double sp=pars[1], h= pars[0]; 
  if(sp<0.|sp>1.){
    GetRNGstate();sp=runif(0.,1.);PutRNGstate();
  }
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    x[i]=tmp[i+2]; 
    w[i] = tmp[i+n+2]; 
    fx[i]=0.0;
  }
  if(h<0.||h>5.*h0){
    GetRNGstate();
    h = runif(0.001,5.0*h0);
    PutRNGstate();
  }
  for(i=0;i<n;i++) hs[i]=1.0;
  awpdf(x,n,w,h,hs,x,n,fx,0);
  UpdateBwfactor(fx,n,sp,hs);
  
  for(i=0;i<n;i++)
    res += w[i]*(fx[i]-w[i]*M_1_SQRT_2PI/h/hs[i])/(1.-w[i]);
  res *= -2.;
  for(i=0;i<n;i++){
    res2 += pi2*pow(w[i],2.0)/h/hs[i];
    for(j=i+1;j<n;j++)
      res2 += 2.*w[i]*w[j]*dnorm(x[i]-x[j],0.,sqrt(hs[i]*hs[i]+hs[j]*hs[j])*h,0);
  }
  return(res+res2);
}


static double wise1(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  double res2=0.0, pi2 = M_2_SQRTPI/4., res=0.0, *tmp= (double*)ex;
  int i,j,n = (int)tmp[0]; //first element is the length of x;
  double x[n],hs[n],w[n],fx[n];
  double h= pars[0],h0; 
  h0 = tmp[1]; //second element is h0
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    x[i]=tmp[i+2]; 
    w[i] = tmp[i+n+2]; 
    fx[i]=0.0;
    hs[i]=1.;
  }
  if(h<0.||h>5.*h0){
    GetRNGstate();
    h = runif(0.001,5*h0);
    PutRNGstate();
  }
  awpdf(x,n,w,h,hs,x,n,fx,0);
  for(i=0;i<n;i++){ res += w[i]*(fx[i]-w[i]*M_1_SQRT_2PI/h)/(1.-w[i]);}
  res *= -2.;
  for(i=0;i<n;i++){
    res2 += pi2*pow(w[i],2.0)/h;
    for(j=i+1;j<n;j++){
      res2 += 2.*w[i]*w[j]*dnorm(x[i]-x[j],0.,M_SQRT2*h,0);
    }
  }
  return(res+res2);
}

void awkde(double *x,double *w,int *xsize,
	   double *y, double *fy, double *Fy, 
	   int *ysize, double *pars, int *RR)
{
  int i,j,k,rr=RR[0],nx=xsize[0],ny=ysize[0],npar=2;
  double h0,h,hs[nx],fx[nx],sp=0.8; //sp=sensitivity parameter;
  double dpar[npar],opar[npar]; 
  h0 = pars[0]; h=h0;
  dpar[0]=h; dpar[1] = sp;
  double abstol=0.000000001,reltol=0.0000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  for(i=0;i<nx;i++){hs[i]=1.;}
 
  double yaux[2*nx+2];
  yaux[0] = nx;yaux[1]=h0;
  for(i=0;i<nx;i++){
    yaux[i+2] = x[i];
    yaux[i+nx+2] = w[i];
  }
  nmmin(npar,dpar,opar,&val,wise2,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  h=opar[0]; sp=opar[1];
  if(sp<0.||sp>1.||ifail != 0){
    sp=0;
    npar=1;dpar[0]=h0;
    nmmin(npar,dpar,opar,&val,wise1,&ifail,abstol,reltol, 
	  (void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
    h=opar[0]; 
    pars[0]=h; pars[1]=0.;
  }else{
    awpdf(x,nx,w,h,hs,x,nx,fx,0);
    UpdateBwfactor(fx,nx,sp,hs);
    pars[0]=h; pars[1]=sp;
  }
  
  awpdf(x,nx,w,h,hs,y,ny,fy,rr);
  awcdf(x,nx,w,h,hs,y,ny,Fy);
}

void wkde(double *x,double *w,int *xsize,
	   double *y, double *fy, double *Fy, 
	   int *ysize, double *pars, int *RR)
{
  int i,j,k,rr=RR[0],nx=xsize[0],ny=ysize[0],npar=2;
  double h0,h,hs[nx],fx[nx],sp=0.8; //sp=sensitivity parameter;
  double dpar[npar],opar[npar]; 
  h0 = pars[0]; h=h0;
  dpar[0]=h; dpar[1] = sp;
  double abstol=0.000000001,reltol=0.0000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  for(i=0;i<nx;i++){hs[i]=1.;}
 
  double yaux[2*nx+2];
  yaux[0] = nx;yaux[1]=h0;
  for(i=0;i<nx;i++){
    yaux[i+2] = x[i];
    yaux[i+nx+2] = w[i];
  }

  sp=0;
  npar=1;dpar[0]=h0;
  nmmin(npar,dpar,opar,&val,wise1,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  h=opar[0]; 
  pars[0]=h; pars[1]=0.;
  
  awpdf(x,nx,w,h,hs,y,ny,fy,rr);
  awcdf(x,nx,w,h,hs,y,ny,Fy);
}

void ckde(double *x,double *w,int *xsize,
	   double *y, double *fy, double *Fy, 
	   int *ysize, double *pars, int *RR)
{
  int i,j,k,rr=RR[0],nx=xsize[0],ny=ysize[0];
  double h,hs[nx],fx[nx];

  h = pars[0]; 
  for(i=0;i<nx;i++){hs[i]=1.;}
  
  awpdf(x,nx,w,h,hs,y,ny,fy,rr);
  awcdf(x,nx,w,h,hs,y,ny,Fy);
}

void wmise(double *x,double *w,int *xsize,
	   double *hs, double *mises, int *hsize)
{
  int i,j,k;
  double xh1, xh2, g1,g2,g3,tmp;

  for(k=0; k<hsize[0]; k++){
    g1 = 0.0; g2 = 0.0; g3 = 0.0;
    for(i=0;i<xsize[0];i++){
      tmp = 0.0;
      for(j=0;j<xsize[0];j++){
	xh1 = (x[i] - x[j])/hs[k];
	xh2 = xh1/1.414214;
	g1 += w[i]*w[j]*dnorm(xh2,0.0,1.0,0);
	tmp += w[j]*dnorm(xh1,0.0,1.0,0);
      }
      g2 += tmp*w[i]/(1.0-w[i]);
      g3 += w[i]*w[i]/(1.0-w[i]);
    }
    mises[k] = g1 - 2.0 * (g2 - g3*M_1_SQRT_2PI); 
  }
}

