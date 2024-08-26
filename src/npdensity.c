#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"

typedef double (*Fun1p)(double);
typedef double (*Fun2p)(double,double);
typedef double (*Fun3p)(double,double,double);
typedef double (*Fun3d)(double,double,double*);
typedef double (*Fun4p)(double,double*,double,int);
typedef double (*Fun5p)(double,double,double,double*,int,int);
typedef double (*Fun6p)(double,double,double,double*,double*,int);

/////////// Defining Kernels ////////////

double KernHLap(double t,double z,double sig){
  double res=0.0;
  if(fabs(z) < 1.0){
    res = cos(-t * z) * exp(-0.5*t*t) * sig;
  }
  return res;
}

//////////////////////////////////////////////////////////////////////////

void BootSample(double x[], double y[], int n)
/* 
To draw a bootstrap sample from y and the results will be stored in x.
 */
{
  int i,j;
  GetRNGstate();
  for(i=0;i<n;i++){
    j = runif(0.0,1.0)*(1.0+n);
    x[i] = y[(int)j];
  }
  PutRNGstate();
}

void BS(double *x,double *y, int *n)
//Test code
{
  BootSample(x,y,n[0]);
}

void rlaplace(double x[], int n, double sig[])
//to draw a random sample from double exponential (Laplacian distribution)
{
  int i; 
  GetRNGstate();
  for(i=0;i<n;i++){
    x[i] = rexp(1./sig[i]);
    if(runif(0.0,1.0)<0.5) x[i]=-x[i];
  }
  PutRNGstate();
}

void rdexp(double *x,int *n, double *sig)
//Test code
{
  rlaplace(x,n[0],sig);
}

//////////////////////////////////////////////////////////////////////////    
// double Gauss_Legendre_Integration_2pts( double a, double b, double (*f)(double) ) 
// void   Gauss_Legendre_Zeros_2pts( double nodes[] )                  
//    void   Gauss_Legendre_Coefs_2pts( double wght[] )                   
//////////////////////////////////////////////////////////////////////////

//  16 pts
static const double B16[] = {
    9.50125098376374401877e-02,    2.81603550779258913231e-01,
    4.58016777657227386350e-01,    6.17876244402643748452e-01,
    7.55404408355003033891e-01,    8.65631202387831743866e-01,
    9.44575023073232576090e-01,    9.89400934991649932601e-01
};

static const double A16[] = {
    1.89450610455068496287e-01,    1.82603415044923588872e-01,
    1.69156519395002538183e-01,    1.49595988816576732080e-01,
    1.24628971255533872056e-01,    9.51585116824927848073e-02,
    6.22535239386478928628e-02,    2.71524594117540948514e-02
};

#define NOPZ16  sizeof(B16) / sizeof(double)
#define NOZ16   NOPZ16+NOPZ16


double GLInteg(double a, double b, double (*f)(double))
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B16[NOPZ16 - 1];
   const double *pA = &A16[NOPZ16 - 1];

   for (; pB >= B16; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*f)(d - dum) + (*f)(d + dum) );
   }

   return c * integral;
}

double GLInt3p(double a, double b, double (*f)(double,double,double),
	       double sig2,double h)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B16[NOPZ16 - 1];
   const double *pA = &A16[NOPZ16 - 1];

   for (; pB >= B16; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*f)(d - dum,sig2,h) + (*f)(d + dum,sig2,h) );
   }

   return c * integral;
}

double GLInt3d(double a, double b, 
	       double (*f)(double,double,double),
	       double x,double *psig, double *nsig)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B16[NOPZ16 - 1];
   const double *pA = &A16[NOPZ16 - 1];
   int k=NOPZ16;

   for (; pB >= B16; pA--, pB--,k--) {
      dum = c * *pB;
      integral += *pA * ((*f)(d-dum,x,nsig[k]) 
			 + (*f)(d + dum,x,psig[k]));
   }

   return c * integral;
}


double GLInt4p(double a, double b, double (*f)(double,double*,double,int),
	       double *ss,double h,int n)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B16[NOPZ16 - 1];
   const double *pA = &A16[NOPZ16 - 1];

   for (; pB >= B16; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*f)(d - dum,ss,h,n) + (*f)(d + dum,ss,h,n) );
   }

   return c * integral;
}

double GLInt6p(double a, double b, 
	       double (*f)(double,double,double,double*,double*,int),
	       double h, double g, double *sig, double *x,int n)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B16[NOPZ16 - 1];
   const double *pA = &A16[NOPZ16 - 1];

   for (; pB >= B16; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*f)(d - dum,h,g,sig,x,n) + (*f)(d + dum,h,g,sig,x,n) );
   }

   return c * integral;
}

//////////////////////////////////////////////////////////////////////////
void DKEGauss(double *y,int *ny,double *x, int *nx, 
		 double *bw,double *sig, int *type)
//type=0 for PDF and type=1 for CDF
{  
  double z[ny[0]], sigh;
  sigh = 1.0 - pow(sig[0]/bw[0],2.0);
  int i,j;
  switch(type[0]){
  case 0:
    for(i=0;i<nx[0];i++){
      for(j=0;j<ny[0];j++){
	z[j] = pow((x[i]-y[j])/bw[0],2.0)/sigh*0.5;
      }
      x[i]=0.0;
      for(j=0;j<ny[0];j++){
	x[i] += exp(-z[j]);
      }
      x[i] = x[i]/ny[0]/bw[0]/sqrt(sigh)*M_1_SQRT_2PI;
      if(x[i]<0.0) x[i]=0.0;
    }
    break;
  case 1:
    sigh = sqrt(sigh);
    for(i=0;i<nx[0];i++){
      for(j=0;j<ny[0];j++){
	z[j] = (x[i]-y[j])/bw[0];
      }
      x[i]=0.0;
      for(j=0;j<ny[0];j++){
	x[i] += pnorm(z[j],0.0,sigh,1,0);
      }
      x[i] = x[i]/ny[0];
      if(x[i]<0.0) x[i]=0.0;
    }
    break;
  default:
    Rprintf("No type is specified!");
  }
}

void DKELaplace(double *y,int *ny,double *x, int *nx, 
		double *bw,double *sig, int *type)
//type=0 for PDF and type=1 for CDF
{  
  double z[ny[0]], sigh2,tmp1,tmp2,sum1,sum2;
  sigh2 = pow(sig[0]/bw[0],2.0);
  int i,j,k=0;
  switch(type[0]){
  case 0:
    for(i=0;i<nx[0];i++){
      for(j=0;j<ny[0];j++){
	z[j] = pow((x[i]-y[j])/bw[0],2.0);
      }
      x[i]=0.0;sum1=0.0;sum2=0.0;
      for(j=0;j<ny[0];j++){
	tmp1 = exp(-z[j]*0.5);
	tmp2 = tmp1*z[j];
	if(sigh2*(1.0 - z[j]) > -1.0){
	  sum1 += tmp1;
	  sum2 += tmp2;
	}
      }
      x[i] = ((1.0+sigh2)*sum1-sigh2*sum2)/ny[0]/bw[0]*M_1_SQRT_2PI;
      if(x[i]<0.0) x[i]=0.0;
    }
    break;
  case 1:
    for(i=0;i<nx[0];i++){
      if(k==0){
	for(j=0;j<ny[0];j++){
	  z[j] = (x[i]-y[j])/bw[0];
	}
	sum1=0.0;sum2=0.0;
	for(j=0;j<ny[0];j++){
	  tmp1 = pnorm(z[j],0.0,1.0,1,0);
	  tmp2 = dnorm(z[j],0.0,1.0,0)*z[j];
	  sum1 += tmp1;
	  sum2 += tmp2;
	}
	x[i] =(sum1+sigh2*sum2)/ny[0];
	if(x[i]<0.0) x[i]=0.0;
	if(x[i]>1.0) {x[i]=1.0;k=1;}
      }else{
	x[i]=1.0;
      }
    }
    break;
  default:
    Rprintf("No type is specified!");
  }
}


void DKESupport(double *y,int *ny,double *x, int *nx, 
		double *bw,double *sig, int *type)
/*
  type=0 for PDF and type=1 for CDF; We use K=10pts Legendre-Gauss
  Quadrature integration method to compute the kernel.
*/
{  
  double z[ny[0]];
  int i,j,k;
  double integral = 0.0; 
  double a=0.0,b=1.0;
  double c = 0.5 * (b - a);
  double d = 0.5 * (b + a);
  double dum, sb2,t,nt,pt;
  int K=NOPZ16;  //K is changable
  double ntexp[K],ptexp[K];  //K is changable
  const double *pB = &B16[K - 1];  //K is changable 
  const double *pA = &A16[K - 1];  //K is changable
  double ppart1,npart1;

   switch(type[0]){
   case 0:
     sb2 = pow(sig[0]/bw[0],2.0)*0.5;
     k=K-1;
     for (; pB >= B16; pA--, pB--,k--) {
       dum = c * *pB;
       nt = pow(d-dum,2.0);
       ntexp[k] = pow(1.0-nt,3.0)*exp(sb2*nt);
       pt = pow(d+dum,2.0);
       ptexp[k] = pow(1.0-pt,3.0)*exp(sb2*pt);
     }
     for(i=0;i<nx[0];i++){
       for(j=0;j<ny[0];j++){
	 z[j] = (x[i]-y[j])/bw[0];
       }
       k = K-1;
       integral = 0.0;
       pB = &B16[K - 1];   
       pA = &A16[K - 1];  
       for (; pB >= B16; pA--, pB--,k--) {
	 dum = c * *pB;
	 nt = d-dum;
	 pt = d+dum;
	 npart1=0.0;ppart1=0.0;
	 for(j=0;j<ny[0];j++){
	   npart1 += cos(nt*z[j]);
	   ppart1 += cos(pt*z[j]);
	 }
	 integral +=  *pA * (npart1*ntexp[k] + ppart1*ptexp[k]);
       }
       x[i] = c*integral/bw[0]/ny[0]*M_1_PI;
       if(x[i]<0.0) x[i]=0.0;
     }
     break;
   case 1:
     sb2 = pow(sig[0]/bw[0],2.0)*0.5;
     k=K-1;
     for (; pB >= B16; pA--, pB--,k--) {
       dum = c * *pB;
       t= d-dum;
       nt = pow(t,2.0);
       ntexp[k] = pow(1.0-nt,3.0)*exp(sb2*nt)/t;
       t=d+dum;
       pt = pow(t,2.0);
       ptexp[k] = pow(1.0-pt,3.0)*exp(sb2*pt)/t;
     }
     for(i=0;i<nx[0];i++){
       for(j=0;j<ny[0];j++){
	 z[j] = (x[i]-y[j])/bw[0];
       }
       k = K-1;
       integral = 0.0;
       pB = &B16[K - 1];   
       pA = &A16[K - 1];  
       for (; pB >= B16; pA--, pB--,k--) {
	 dum = c * *pB;
	 nt = d-dum;
	 pt = d+dum;
	 npart1=0.0;ppart1=0.0;
	 for(j=0;j<ny[0];j++){
	   npart1 += sin(nt*z[j]);
	   ppart1 += sin(pt*z[j]);
	 }
	 integral +=  *pA * (npart1*ntexp[k] + ppart1*ptexp[k]);
       }
       x[i] = 0.5 + c*integral/ny[0]*M_1_PI;
       if(x[i]<0.0) x[i]=0.0;
     }
     break;
  default:
    Rprintf("No type is specified!");
  }
}

void densHSupport(double y[],int ny,double x[], int nx, 
		double bw,double sig[])
{  
  int i,j,k;
  double integral = 0.0; 
  double a=0.0,b=1.0;
  double c = 0.5 * (b - a);
  double d = 0.5 * (b + a);
  double dum;
  int K=NOPZ16;  //K is changable
  double nt[K], pt[K], nt2[K], pt2[K],nt3[K],pt3[K];
  double sigh[ny], z[ny];
  double nsig[K][ny], psig[K][ny];
  const double *pB = &B16[K - 1];  //K is changable 
  const double *pA = &A16[K - 1];  //K is changable
  double nsum, psum;

  nsum = bw*bw;//reuse valuable
  for(j=0;j<ny;j++){
    sigh[j] = -0.5 * pow(sig[j],2.0)/nsum;
  }
  k=K-1;
  for (; pB >= B16; pA--, pB--,k--) {
    dum = c * *pB;
    nt[k] = d-dum;
    nt2[k] = pow(nt[k],2.0);
    nt3[k] = pow(1.0-nt2[k],3.0);
    pt[k] = d+dum;
    pt2[k] = pow(pt[k],2.0);
    pt3[k] = pow(1.0-pt2[k],3.0);
    nsum = 0.0;psum=0.0;//sum of exp part;
    for(j=0;j<ny;j++){
      nsig[k][j] = exp(sigh[j]*nt2[k]);
      psig[k][j] = exp(sigh[j]*pt2[k]);
      nsum += pow(nsig[k][j],2.0);
      psum += pow(psig[k][j],2.0);
    }
    for(j=0;j<ny;j++){
      nsig[k][j] /= nsum;
      psig[k][j] /= psum;
    }
  }
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      z[j] = (x[i]-y[j])/bw;
    }
    k = K-1;
    integral = 0.0;
    pB = &B16[K - 1];   
    pA = &A16[K - 1];  
    for (; pB >= B16; pA--, pB--,k--) {
      nsum=0.0;psum=0.0;
      for(j=0;j<ny;j++){
	nsum += cos(nt[k]*z[j]) * nsig[k][j];
	psum += cos(pt[k]*z[j]) * psig[k][j];
      }
      integral +=  *pA * (nsum*nt3[k] + psum*pt3[k]);
    }
    x[i] = c*integral/bw*M_1_PI;
    if(x[i]<0.0) x[i]=0.0;
  }
}


void DKEHSupport(double *y,int *ny,double *x, int *nx, 
		double *bw,double *sig, int *type)
/*
  type=0 for PDF and type=1 for CDF; We use K=10pts Legendre-Gauss
  Quadrature integration method to compute the kernel.
*/
{  
  int i,j,k;
  double integral = 0.0; 
  double a=0.0,b=1.0;
  double c = 0.5 * (b - a);
  double d = 0.5 * (b + a);
  double dum;
  int K=NOPZ16;  //K is changable
  double nt[K], pt[K], nt2[K], pt2[K],nt3[K],pt3[K];
  double sigh[ny[0]], z[ny[0]];
  double nsig[K][ny[0]], psig[K][ny[0]];
  const double *pB = &B16[K - 1];  //K is changable 
  const double *pA = &A16[K - 1];  //K is changable
  double nsum, psum;

  nsum = bw[0]*bw[0];//reuse valuable
  for(j=0;j<ny[0];j++){
    sigh[j] = -0.5 * pow(sig[j],2.0)/nsum;
  }
  k=K-1;
  for (; pB >= B16; pA--, pB--,k--) {
    dum = c * *pB;
    nt[k] = d-dum;
    nt2[k] = pow(nt[k],2.0);
    nt3[k] = pow(1.0-nt2[k],3.0);
    pt[k] = d+dum;
    pt2[k] = pow(pt[k],2.0);
    pt3[k] = pow(1.0-pt2[k],3.0);
    nsum = 0.0;psum=0.0;//sum of exp part;
    for(j=0;j<ny[0];j++){
      nsig[k][j] = exp(sigh[j]*nt2[k]);
      psig[k][j] = exp(sigh[j]*pt2[k]);
      nsum += pow(nsig[k][j],2.0);
      psum += pow(psig[k][j],2.0);
    }
    for(j=0;j<ny[0];j++){
      nsig[k][j] /= nsum;
      psig[k][j] /= psum;
    }
  }
   switch(type[0]){
   case 0:
     for(i=0;i<nx[0];i++){
       for(j=0;j<ny[0];j++){
	 z[j] = (x[i]-y[j])/bw[0];
       }
       k = K-1;
       integral = 0.0;
       pB = &B16[K - 1];   
       pA = &A16[K - 1];  
       for (; pB >= B16; pA--, pB--,k--) {
	 nsum=0.0;psum=0.0;
	 for(j=0;j<ny[0];j++){
	   nsum += cos(nt[k]*z[j]) * nsig[k][j];
	   psum += cos(pt[k]*z[j]) * psig[k][j];
	 }
	 integral +=  *pA * (nsum*nt3[k] + psum*pt3[k]);
       }
       x[i] = c*integral/bw[0]*M_1_PI;
       if(x[i]<0.0) x[i]=0.0;
     }
     break;
   case 1:
     for(i=0;i<nx[0];i++){
       for(j=0;j<ny[0];j++){
	 z[j] = (x[i]-y[j])/bw[0];
       }
       k = K-1;
       integral = 0.0;
       pB = &B16[K - 1];   
       pA = &A16[K - 1];  
       for (; pB >= B16; pA--, pB--,k--) {
	 nsum=0.0;psum=0.0;
	 for(j=0;j<ny[0];j++){
	   nsum += sin(nt[k]*z[j]) * nsig[k][j];
	   psum += sin(pt[k]*z[j]) * psig[k][j];
	 }
	 integral +=  *pA * (nsum*nt3[k]/nt[k] + psum*pt3[k]/pt[k]);
       }
       x[i] = 0.5 + c*integral*M_1_PI;
       if(x[i]<0.0) x[i]=0.0;
     }
     break;
  default:
    Rprintf("No type is specified!");
  }
}

// Bandwidth selection
void NormNorm1(int *n, double *Rfx,double *s2, double *h1,double *grid,double *ub)
//Kernel=normal; Error=normal; homoscedastic=yes;
{
  double h=fmax(h1[0]/ub[0],sqrt(s2[0])), hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999.,hdiff, hopt=0.0;
  int i;
  for(i=0;i<30;i++){
    h += hstep;
    hdiff = h*h-s2[0];
    mise = 0.5/n[0]/pow(M_PI*hdiff,-0.5)+Rfx[0]*pow(h,4.0);
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}

void NormLap1(int *n, double *Rfx,double *s2, double *h1,double *grid,double *ub)
//Kernel=normal; Error=normal; homoscedastic=yes;
{
  double h=h1[0]/ub[0], hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999., hopt=0.0;
  int i;
  for(i=0;i<(int)grid[0];i++){
    h += hstep;
    mise = M_1_SQRT_2PI/n[0]/h*(1.+2.*s2[0]/pow(h,2.0)+3.*pow(s2[0],2.0)/pow(h,4.0))
      +Rfx[0]*pow(h,4.0);
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}

double funSuppNorm1(double t,double sig2,double h){
  return pow(1.0-pow(t,2.0),6.0)*exp(sig2*pow(t/h,2.0));
}

void SuppNorm1(int *n, double *Rfx,double *s2, double *h1,double *grid,double *ub)
//Kernel=normal; Error=normal; homoscedastic=yes;
{
  double h=h1[0]/ub[0], hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999.,fint, hopt=0.0;
  Fun3p f[1];
  f[0] = funSuppNorm1;
  int i;
  for(i=0;i<(int)grid[0];i++){
    h += hstep;
    fint = GLInt3p(0.,1.0,f[0],s2[0],h);
    mise = M_1_PI/n[0]/h*fint + 6.0*Rfx[0]*pow(h,4.0);
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}

double funSuppLap1(double t,double sig2,double h){
  double t2=pow(t,2.0);
  return pow(1.0-t2,6.0)*pow(1.+sig2*t2/pow(h,2.0),2.0);
}

void SuppLap1(int *n, double *Rfx,double *s2, double *h1,double *grid,double *ub)
//Kernel=normal; Error=normal; homoscedastic=yes;
{
  double h=h1[0]/ub[0], hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999.,fint, hopt=0.0;
  Fun3p f[1];
  f[0] = funSuppLap1;
  int i;
  for(i=0;i<(int)grid[0];i++){
    h += hstep;
    fint = GLInt3p(0.,1.0,f[0],s2[0],h);
    mise = M_1_PI/n[0]/h*fint + 6.0*Rfx[0]*pow(h,4.0);
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}


double funSuppNorm2(double t,double *ss,double h,int n){
  double t2=pow(t,2.0), h2=pow(h,2.0);
  double fsum=0.0;
  int i;
  double *ptrss;
  ptrss = &ss[0];

  for(i=0;i<n;i++){
    fsum += exp(-(*ptrss)*t2/h2);
  }
  return pow(1.-t2,6.0)/fsum;
}

void SuppNorm2(int *n, double *Rfx,double *ss, double *h1,double *grid,double *ub)
//Kernel=normal; Error=normal; homoscedastic=yes;
{
  double h=h1[0]/ub[0], hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999.,fint, hopt=0.0;
  Fun4p f[1];
  f[0] = funSuppNorm2;
  int i;
  for(i=0;i<(int)grid[0];i++){
    h += hstep;
    fint = GLInt4p(0.,1.0,f[0],ss,h,n[0]);
    mise = M_1_PI/h*fint + 6.0*Rfx[0]*pow(h,4.0);
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}


double funSuppLap2(double t,double *ss,double h,int n){
  double t2=pow(t,2.0), h2=pow(h,2.0);
  double fsum=0.0;
  int i;
  double *ptrss;
  ptrss = &ss[0];

  for(i=0;i<n;i++){
    fsum += pow(1.0+(*ptrss)*t2/h2,-2.0);
  }
  return pow(1.-t2,6.0)/fsum;
}

void SuppLap2(int *n, double *Rfx,double *ss, double *h1,double *grid,double *ub)
//Kernel=normal; Error=normal; homoscedastic=yes;
{
  double h=h1[0]/ub[0], hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999.,fint, hopt=0.0;
  Fun4p f[1];
  f[0] = funSuppLap2;
  int i;
  for(i=0;i<(int)grid[0];i++){
    h += hstep;
    fint = GLInt4p(0.,1.0,f[0],ss,h,n[0]);
    mise = M_1_PI/h*fint + 6.0*Rfx[0]*pow(h,4.0);
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}


double funNormNorm2(double t,double *ss,double h,int n){
  double t2=pow(t,2.0), h2=pow(h,2.0);
  double fsum=0.0;
  int i;
  double *ptrss;
  ptrss = &ss[0];

  for(i=0;i<n;i++){
    fsum += exp(-(*ptrss)*t2/h2);
  }
  return exp(-t2)/fsum;
}

void NormNorm2(int *n, double *Rfx,double *ss, double *h1,double *grid,double *ub)
//Kernel=normal; Error=normal; homoscedastic=yes;
{
  double h=h1[0]/ub[0], hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999.,fint, hopt=0.0;
  Fun4p f[1];
  f[0] = funNormNorm2;
  int i;
  for(i=0;i<(int)grid[0];i++){
    h += hstep;
    fint = GLInt4p(0.,5.,f[0],ss,h,n[0]);
    mise = M_1_PI/h*fint + Rfx[0]*pow(h,4.0);
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}

double funNormLap2(double t,double *ss,double h,int n){
  double t2=pow(t,2.0), h2=pow(h,2.0);
  double fsum=0.0;
  int i;
  double *ptrss;
  ptrss = &ss[0];

  for(i=0;i<n;i++){
    fsum += pow(1.0+(*ptrss)*t2/h2,-2.0);
  }
  return exp(-t2)/fsum;
}

void NormLap2(int *n, double *Rfx,double *ss, double *h1,
	      double *grid,double *ub)
//Kernel=normal; Error=laplace; homoscedastic=no;
{
  double h=h1[0]/ub[0], hstep=(ub[0]-1./ub[0])*h1[0]/grid[0];
  double mise,mmin=99999999999.,fint, hopt=0.0;
  Fun4p f[1];
  f[0] = funNormLap2;
  int i;
  for(i=0;i<(int)grid[0];i++){
    h += hstep;
    fint = GLInt4p(0.,5.,f[0],ss,h,n[0]);
    mise = M_1_PI/h*fint + Rfx[0]*pow(h,4.0);
    //printf("MISE = %f10.2, bw = %f10.2\n", mise , h );
    if(mise<mmin) {
      hopt=h; mmin=mise;
    }
  }
  h1[0]=hopt;
}

double KernGL(double x, double sigh){
  return(dnorm(x,0.0,1.0,0)*(1.0+pow(sigh,2.0)*(1.0-pow(x,2.0))));
}

double dknpreg(double x, double Z[], double Y[], double S[], int n, double h){
  int i;
  double usum=0.,dsum=0.,tmp=0.0;
  for(i=0;i<n;i++){
    tmp = KernGL((x-Z[i])/h, S[i]);
    usum += Y[i] * tmp;
    dsum += tmp;
  }
  return(usum/dsum);
}

double dknpreg2(int m, double Z[], double Y[], double S[], int n, double h){
  int i;
  double usum=0.,dsum=0.,tmp=0.0,x;

  x = Z[m];
  for(i=0;i<n;i++){
    if(i != m){
      tmp = KernGL((x-Z[i])/h, S[i]);
      usum += Y[i] * tmp;
      dsum += tmp;
    }
  }
  return(usum/dsum);
}

void DkNpReg(double *Z,double *Y, double *S, int *size, 
	     double *bandwidth, double *X, int *nx, 
	     double *loo, double *opt){
  int i,j,n=size[0];
  double hopt=bandwidth[0];
  for(i=0;i<n;i++){
    S[i] /= hopt;
  }  

  // choose optimal bandwidth using the risk as the leave-one-out
  // cross-validation score:
  double h,hstep,rh,t,Rh=1.0e9;
  if(opt[0] > 0.0){
    h = bandwidth[0]*0.8;
    hstep = 0.0035*bandwidth[0];
    for(j=0;j<400;j++){
      rh = 0.0;
      for(i=0;i<size[0];i++){
	if(loo[0] > 0.0){
	  t = dknpreg2(i,Z,Y,S,n,h);
	}else{
	  t = dknpreg(Z[i],Z,Y,S,n,h);
	}
	rh += pow(t-Y[i],2.0);
      }
      rh = rh / size[0];
      if(rh < Rh){
	Rh = rh;
	hopt = h;
      }
      h += hstep;
    }  
    bandwidth[0] = hopt;
    opt[0] = Rh;
  }else{
    rh = 0.0;
    for(i=0;i<size[0];i++){
      if(loo[0] > 0.0){
	t = dknpreg2(i,Z,Y,S,n,hopt);
      }else{
	t = dknpreg(Z[i],Z,Y,S,n,hopt);
      }
      //t = dknpreg2(i,Z,Y,S,n,hopt);
      //t = dknpreg(Z[i],Z,Y,S,n,hopt);
      rh += pow(t-Y[i],2.0);
    }
    rh = rh / size[0];
    opt[0] = rh;
  }

  for(i=0;i<nx[0];i++){
    X[i] = dknpreg(X[i],Z,Y,S,n,hopt);
  }
}

/*    --------------------------------------------------------------
The following codes are for nonparametric regression
------------------------------------------------------------------ */


void NPRGauss(double *y,int *ny,double *zs,double *x, int *nx, 
		 double *bw,double *sig)
{  
  double z[ny[0]], sigh,tmp0,tmp1;
  sigh = 1.0 - pow(sig[0]/bw[0],2.0);
  int i,j;
  for(i=0;i<nx[0];i++){
    for(j=0;j<ny[0];j++){
      z[j] = pow((x[i]-y[j])/bw[0],2.0)/sigh*0.5;
    }
    x[i]=0.0; tmp0=0.;
    for(j=0;j<ny[0];j++){
      tmp1 = exp(-z[j]);
      x[i] += tmp1;
      tmp0 += tmp1*zs[j];
    }
    x[i] = tmp0/x[i];
  }
}

void NPRLaplace(double *y,int *ny,double *zs,double *x, int *nx, 
		double *bw,double *sig)
{  
  double z[ny[0]], sigh2,tmp1,tmp2,sum1,sum2;
  sigh2 = pow(sig[0]/bw[0],2.0);
  int i,j;
  for(i=0;i<nx[0];i++){
    for(j=0;j<ny[0];j++){
      z[j] = pow((x[i]-y[j])/bw[0],2.0);
    }
    x[i]=0.0;sum1=0.0;sum2=0.0;
    for(j=0;j<ny[0];j++){
      tmp1 = exp(-z[j]*0.5)*(1.+sigh2*(1.-z[j]));
      tmp2 = tmp1*zs[j];
      sum1 += tmp1;
      sum2 += tmp2;
    }
    x[i] = sum2/sum1;
  }
}

void NPRSupport(double *y,int *ny,double *zs,double *x, int *nx, 
		double *bw,double *sig)
{  
  //*zs is the dependent variable
  double z[ny[0]];
  int i,j,k;
  double integral = 0.0; 
  double a=0.0,b=1.0;
  double c = 0.5 * (b - a);
  double d = 0.5 * (b + a);
  double dum, sb2,nt,pt;
  int K=NOPZ16;  //K is changable
  double ntexp[K],ptexp[K];  //K is changable
  const double *pB = &B16[K - 1];  //K is changable 
  const double *pA = &A16[K - 1];  //K is changable
  double ppart1,npart1,denom,tmp0,ppart2,npart2;
  sb2 = pow(sig[0]/bw[0],2.0)*0.5;
  k=K-1;
  for (; pB >= B16; pA--, pB--,k--) {
    dum = c * *pB;
    nt = pow(d-dum,2.0);
    ntexp[k] = pow(1.0-nt,3.0)*exp(sb2*nt);
    pt = pow(d+dum,2.0);
    ptexp[k] = pow(1.0-pt,3.0)*exp(sb2*pt);
  }
  for(i=0;i<nx[0];i++){
    for(j=0;j<ny[0];j++){
      z[j] = (x[i]-y[j])/bw[0];
    }
    k = K-1;
    integral = 0.0;denom=0.0;
    pB = &B16[K - 1];   
    pA = &A16[K - 1];  
    for (; pB >= B16; pA--, pB--,k--) {
      dum = c * *pB;
      nt = d-dum;
      pt = d+dum;
      npart1=0.0;ppart1=0.0;npart2=0.;ppart2=0.;
      for(j=0;j<ny[0];j++){
	npart1 += cos(nt*z[j]);
	ppart1 += cos(pt*z[j]);
	npart2 += cos(nt*z[j]) * zs[j];
	ppart2 += cos(pt*z[j]) * zs[j];
      }
      tmp0 = *pA * (npart1*ntexp[k] + ppart1*ptexp[k]);
      denom += tmp0;
      tmp0 = *pA * (npart2*ntexp[k] + ppart2*ptexp[k]);
      integral +=  tmp0;
    }
    x[i] = integral/denom;
  }
}

void NPRHSupport(double *y,int *ny,double *zs,double *x, int *nx, 
		double *bw,double *sig)
{  
  int i,j,k;
  double integral = 0.0; 
  double a=0.0,b=1.0;
  double c = 0.5 * (b - a);
  double d = 0.5 * (b + a);
  double dum;
  int K=NOPZ16;  //K is changable
  double nt[K], pt[K], nt2[K], pt2[K],nt3[K],pt3[K];
  double sigh[ny[0]], z[ny[0]];
  double nsig[K][ny[0]], psig[K][ny[0]];
  const double *pB = &B16[K - 1];  //K is changable 
  const double *pA = &A16[K - 1];  //K is changable
  double nsum, psum,denom,tmp0,nsum2,psum2;

  nsum = bw[0]*bw[0];//reuse valuable
  for(j=0;j<ny[0];j++){
    sigh[j] = -0.5 * pow(sig[j],2.0)/nsum;
  }
  k=K-1;
  for (; pB >= B16; pA--, pB--,k--) {
    dum = c * *pB;
    nt[k] = d-dum;
    nt2[k] = pow(nt[k],2.0);
    nt3[k] = pow(1.0-nt2[k],3.0);
    pt[k] = d+dum;
    pt2[k] = pow(pt[k],2.0);
    pt3[k] = pow(1.0-pt2[k],3.0);
    nsum = 0.0;psum=0.0;//sum of exp part;
    for(j=0;j<ny[0];j++){
      nsig[k][j] = exp(sigh[j]*nt2[k]);
      psig[k][j] = exp(sigh[j]*pt2[k]);
      nsum += pow(nsig[k][j],2.0);
      psum += pow(psig[k][j],2.0);
    }
    for(j=0;j<ny[0];j++){
      nsig[k][j] /= nsum;
      psig[k][j] /= psum;
    }
  }

  for(i=0;i<nx[0];i++){
    for(j=0;j<ny[0];j++){
      z[j] = (x[i]-y[j])/bw[0];
    }
    k = K-1;
    integral = 0.0;denom=0.0;
    pB = &B16[K - 1];   
    pA = &A16[K - 1];  
    for (; pB >= B16; pA--, pB--,k--) {
      nsum=0.0;psum=0.0;nsum2=0.;psum2=0.;
      for(j=0;j<ny[0];j++){
	nsum += cos(nt[k]*z[j]) * nsig[k][j];
	psum += cos(pt[k]*z[j]) * psig[k][j];
	nsum2 += cos(nt[k]*z[j]) * nsig[k][j] * zs[j];
	psum2 += cos(pt[k]*z[j]) * psig[k][j] * zs[j];
      }
      tmp0 = *pA * (nsum*nt3[k] + psum*pt3[k]);
      denom += tmp0;
      tmp0 = *pA * (nsum2*nt3[k] + psum2*pt3[k]);
      integral +=  tmp0;
    }
    x[i] = integral/denom;
  }
}

//Bootstrap-type bandwidth selector

double BootHomoSupp(double t,double h, double g, double *sig, double *x,int n){
  int i;
  double mucos=0.,phi2,ht2,t2;
  ht2 = 1.0-pow(h*t,2.0);
  t2 = pow(t,2.0);
  for(i=0;i<n;i++){
    mucos += cos(x[i]*t);
  }
  mucos = pow(mucos/n,2.0);
  phi2 = pow(1.0-g*g*t*t,6.0)*mucos*exp(pow(sig[0]*t,2.0));
  return pow(1.-t2,6.0)*exp(pow(sig[0]*t/h,2.0))/n/h*0.5 
    - phi2*pow(ht2,3.0) +(n-1.0)*0.5/n*phi2*pow(ht2,6.0);
}

double BootHeteroSupp(double t,double h, double g, double *sig, double *x,int n){
  int i;
  double mucos=0.,phi2,ht2,t2,tsig2,fsum1,fsum2,rn,rd,ratio,sumexp;
  ht2 = 1.0-pow(h*t,2.0);
  t2 = pow(t,2.0);
  fsum1 = 0.0;fsum2=0.;sumexp=0.0;rn=0.0;rd=0.0;ratio=0.0;
  for(i=0;i<n;i++){
    tsig2 = pow(t*sig[i],2.0);
    mucos = exp(-tsig2*0.5); //used as temporary variable
    fsum2 += mucos*mucos;
    fsum1 += mucos*cos(t*x[i]);
    sumexp += exp(-tsig2/h/h);
    rn += pow(mucos,4.0); 
    rd += pow(mucos,2.0);
  }
  mucos = pow(fsum1/fsum2,2.0);
  ratio = rn/pow(rd,2.0);
  phi2 = pow(1.0-g*g*t*t,6.0)*mucos;
  return pow(1.-t2,6.0)/sumexp/h*0.5 
    - phi2*pow(ht2,3.0) +(n-1.0)*0.5*phi2*pow(ht2,6.0)*ratio;
}

double BootHomoNorm(double t,double h, double g, double *sig, double *x,int n){
  int i;
  double mucos=0.,phi2,ht2,tsig2,gt2;
  ht2 = pow(h*t,2.0);
  tsig2 = pow(t*sig[0],2.0);
  gt2 = pow(g*t,2.0);
  for(i=0;i<n;i++){
    mucos = cos(t*x[i]); //used as temporary variable
  }
  mucos = mucos*(1.+tsig2)*exp(-0.5*gt2);
  phi2 = pow(mucos,2.0);
  return exp(-t*t)/n/h*pow(1.+tsig2/h/h,2.0)-2.*phi2*exp(-0.5*ht2)
    +(n-1)/n*phi2*exp(-ht2);
}

double BootHeteroNorm(double t,double h, double g, double *sig, double *x,int n){
  int i;
  double mucos=0.,phi2,ht2,tsig2,gt2,fsum1,fsum2,rn,rd,ratio,sumexp;
  ht2 = pow(h*t,2.0);
  gt2 = pow(g*t,2.0);
  fsum1 = 0.0;fsum2=0.;sumexp=0.0;;rn=0.0;rd=0.0;ratio=0.0;
  for(i=0;i<n;i++){
    tsig2 = pow(t*sig[i],2.0);
    sumexp += pow(1.0+tsig2/h/h,-2.0);
    mucos = 1./(1.0+tsig2);
    fsum1 += mucos*cos(t*x[i]);
    fsum2 += mucos*mucos;
    rn += pow(mucos,4.0); 
    rd += pow(mucos,2.0);
  }
  mucos = fsum1/fsum2;
  ratio = rn/pow(rd,2.0);
  phi2 = pow(mucos*exp(-0.5*gt2),2.0);
  return exp(-t*t)/h/sumexp-2.*phi2*exp(-0.5*ht2)
    +(n-1)*phi2*exp(-ht2)*ratio;
}

void bwBoot1(double *h0,int *size,int *type,double *y,double *sig, 
	     int *grid,double *ub)
{
  double g=h0[0],h=h0[0]/ub[0], hstep=(ub[0]-1./ub[0])*h0[0]/grid[0];
  double mise,mmin=99999999999., fint, hopt=0.0;
  Fun6p f[4];
  f[0] = BootHomoSupp;
  f[1] = BootHeteroSupp;
  f[2] = BootHomoNorm;
  f[3] = BootHeteroNorm;
  int i,n=size[0];
  switch(type[0]){
  case 1:
    for(i=0;i<grid[0];i++){
      h += hstep;
      fint = GLInt6p(-1.0,1.0,f[0],h,g,sig,y,n);
      mise=fint/M_PI;
      if(mise<mmin) {
	hopt=h; mmin=mise;
      }
    }
    break;
  case 2:
    for(i=0;i<grid[0];i++){
      h += hstep;
      fint = GLInt6p(0,4.0,f[2],h,g,sig,y,n);//approximate 
      mise=fint/M_PI;
      if(mise<mmin) {
	hopt=h; mmin=mise;
      }
    }
    break;
  case 3:
    for(i=0;i<grid[0];i++){
      h += hstep;
      fint = GLInt6p(-1.0,1.0,f[1],h,g,sig,y,n);
      mise=fint/M_PI;
      if(mise<mmin) {
	hopt=h; mmin=mise;
      }
    }
    break;
  case 4:
    for(i=0;i<grid[0];i++){
      h += hstep;
      fint = GLInt6p(0,4.0,f[3],h,g,sig,y,n);//approximate 
      mise=fint/M_PI;
      if(mise<mmin) {
	hopt=h; mmin=mise;
      }
    }
    break;
  default:
    Rprintf("The bandwidth selector for this type has not been implemented yet!");
  }
  h0[0]=hopt;
}
 
// 10/11/2013: SCB estimation
void dkdenest(double *y, double *sig, int *ny,
	      double *x, double *Sn, int *nx, 
	      double *bw, int *ktype)
/*
  type=0 for PDF and type=1 for CDF; We use K=10pts Legendre-Gauss
  Quadrature integration method to compute the kernel.
*/
{  
  double z[ny[0]];
  int i,j,k;
  double integral = 0.0; 
  double a=0.0,b=1.0;
  double c = 0.5 * (b - a);
  double d = 0.5 * (b + a);
  double dum, sb2,nt,pt;
  int K=NOPZ16;  //K is changable
  double ntexp[K],ptexp[K];  //K is changable
  const double *pB = &B16[K - 1];  //K is changable 
  const double *pA = &A16[K - 1];  //K is changable
  double ppart1,npart1;

  //Fun4p f[1];
  //f[0] = funSuppNorm2;

  switch(ktype[0]){
  case 0:
    sb2 = pow(sig[0]/bw[0],2.0)*0.5;
    k=K-1;
    for (; pB >= B16; pA--, pB--,k--) {
      dum = c * *pB;
      nt = pow(d-dum,2.0);
      ntexp[k] = pow(1.0-nt,3.0)*exp(sb2*nt);
      pt = pow(d+dum,2.0);
      ptexp[k] = pow(1.0-pt,3.0)*exp(sb2*pt);
    }
    for(i=0;i<nx[0];i++){
      for(j=0;j<ny[0];j++){
	z[j] = (x[i]-y[j])/bw[0];
      }
      k = K-1;
      integral = 0.0;
      pB = &B16[K - 1];   
      pA = &A16[K - 1];  
      for (; pB >= B16; pA--, pB--,k--) {
	dum = c * *pB;
	nt = d-dum;
	pt = d+dum;
	npart1=0.0;ppart1=0.0;
	for(j=0;j<ny[0];j++){
	  npart1 += cos(nt*z[j]);
	  ppart1 += cos(pt*z[j]);
	}
	integral +=  *pA * (npart1*ntexp[k] + ppart1*ptexp[k]);
      }
      x[i] = c*integral/bw[0]/ny[0]*M_1_PI;
      if(x[i]<0.0) x[i]=0.0;
    }
    Sn[0] = 1.0;
    break;
  default:
    Rprintf("No type is specified!");
  }
}

/*  The following codes are added on 10/12/2013: to compute bootstrap
    simutaneous confidence band for density estimates. */

double funSupport(double t,double x,double sig){
  double t2=t*t, res;
  res = cos(t*x)*pow(1.0-t2,3.0)*exp(0.5*sig*sig*t2);
  return res;
}

/* To compute the integral of a deconvoluting kernel for
   heteroscedastic errors. */
double funHKernel(double t,double z,double sig){
  return cos(t * z) * sig;
}

void pdfHSupport(double *x0, double *Sn, int *nx0, 
		 double *x, double *sig, 
		 double *bw, int *size)
{
  Fun3p f[1];
  f[0] = funHKernel;
  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double fx[m],s_h[n],z,h;
  double tmp;
  
  h = bw[0];

  for(i=0;i<m;i++){
    fx[i] = 0.0;
    Sn[i] = 0.0;
  }

  for(i=0;i<n;i++){
    s_h[i] = -0.5 * pow(sig[i]/h, 2.0);
  }
  
  /* nt[] and pt[] will be used as gloabl variable within this
     function. */
  int k, K=NOPZ16;  //K is changable (Integral & here)
  double nt, nt2, nt3, pt, pt2, pt3, dum, psum, nsum;
  double nsig[K][n], psig[K][n];
  double nsig2[K], psig2[K]; //to be pass co compute integrals  
  const double *pB = &B16[K - 1];  //K is changable 
  //  const double *pA = &A16[K - 1];  //K is changable

  k=K-1;
  for (; pB >= B16; pB--,k--) {
    dum = 0.5 * *pB;

    nt = 0.5 - dum;
    nt2 = pow(nt,2.0);
    nt3 = pow(1.0-nt2,3.0);

    pt = 0.5 + dum;
    pt2 = pow(pt,2.0);
    pt3 = pow(1.0-pt2,3.0);

    nsum = 0.0;psum=0.0;//sum of exp part;
    for(j=0;j<n;j++){
      nsig[k][j] = exp(s_h[j]*nt2);
      psig[k][j] = exp(s_h[j]*pt2);
      nsum += pow(nsig[k][j],2.0);
      psum += pow(psig[k][j],2.0);
    }
    for(j=0;j<n;j++){
      nsig[k][j] = nsig[k][j]/nsum*n*nt3;
      psig[k][j] = psig[k][j]/psum*n*pt3;
    }
  }

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      z = (x0[i]-x[j])/h;
      for(k=0; k<K; k++){
	psig2[k] = psig[k][j];
	nsig2[k] = nsig[k][j];
      }
      tmp = GLInt3d(0.0,1.0,f[0],z, psig2, nsig2);
      fx[i] += tmp;
      Sn[i] += tmp * tmp;
    }
  }
  for(i=0;i<m;i++){
    x0[i] = fx[i] /(n * h);
    if(x0[i] < 0.0) x0[i] = 0.0;
    Sn[i] = sqrt(Sn[i]/n)/h; 
  }
}

double lkernel(double x,double sigh){
  double out;
  out = dnorm(x,0.0,1.0,0)* (1.0 + sigh * (1.0 - x*x));
  if(out<0) out = 0.0;
  return out;
}

void pdfHLaplace(double *x0, double *Sn, int *nx0, 
		 double *x, double *sig, 
		 double *bw, int *size)
{
  Fun3p f[1];
  f[0] = funHKernel;
  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double fx[m],s_h[n],z,h;
  double tmp;
  
  h = bw[0];

  for(i=0;i<m;i++){
    fx[i] = 0.0;
    Sn[i] = 0.0;
  }

  for(i=0;i<n;i++){
    s_h[i] = pow(sig[i]/h, 2.0);
  }
  
  /* nt[] and pt[] will be used as gloabl variable within this
     function. */
  int k, K=NOPZ16;  //K is changable (Integral & here)
  double nt, nt2, nt3, pt, pt2, pt3, dum, psum, nsum;
  double nsig[K][n], psig[K][n];
  double nsig2[K], psig2[K]; //to be pass co compute integrals  
  const double *pB = &B16[K - 1];  //K is changable 
  //  const double *pA = &A16[K - 1];  //K is changable

  k=K-1;
  for (; pB >= B16; pB--,k--) {
    dum = 0.5 * *pB;

    nt = 0.5 - dum;
    nt2 = pow(nt,2.0);
    nt3 = exp(-0.5 * nt2);

    pt = 0.5 + dum;
    pt2 = pow(pt,2.0);
    pt3 = exp(-0.5 * pt2);

    nsum = 0.0;psum=0.0;//sum of exp part;
    for(j=0;j<n;j++){
      nsig[k][j] = 1.0/(1.0 + s_h[j]*nt2);
      psig[k][j] = 1.0/(1.0 + s_h[j]*pt2);
      nsum += pow(nsig[k][j],2.0);
      psum += pow(psig[k][j],2.0);
    }
    for(j=0;j<n;j++){
      nsig[k][j] = nsig[k][j]/nsum*n*nt3;
      psig[k][j] = psig[k][j]/psum*n*pt3;
    }
  }

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      z = (x0[i]-x[j])/h;
      for(k=0; k<K; k++){
	psig2[k] = psig[k][j];
	nsig2[k] = nsig[k][j];
      }
      tmp = GLInt3d(0.0,1.0,f[0],z, psig2, nsig2);
      fx[i] += tmp;
      Sn[i] += tmp * tmp;
    }
  }
  for(i=0;i<m;i++){
    x0[i] = fx[i] /(n * h);
    if(x0[i] < 0.0) x0[i] = 0.0;
    Sn[i] = sqrt(Sn[i]/n)/h; 
  }
}

void pdfSupport(double *x0, double *Sn, int *nx0, 
		 double *x, double *sig, 
		 double *bw, int *size)
{
  Fun3p f[1];
  f[0] = funSupport;
  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double fx[m],s_h,z,h, tmp;

  h = bw[0];
  s_h = sig[0]/h;

  for(i=0;i<m;i++){
    fx[i] = 0.0;
    Sn[i] = 0.0;
  }
  
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      z = (x0[i]-x[j])/h;
      tmp = GLInt3p(0.0,1.0,f[0],z,s_h);
      fx[i] += tmp;
      Sn[i] += tmp * tmp;
    }
  }
  for(i=0;i<m;i++){
    tmp = h*M_PI;
    x0[i] = fx[i] /(n * tmp);
    if(x0[i] < 0.0) x0[i] = 0.0;
    Sn[i] = sqrt(Sn[i]/n)/tmp; 
  }
}


void pdfLaplace(double *x0, double *Sn, int *nx0, 
		 double *x, double *sig, 
		 double *bw, int *size)
{
  int i,j,m,n;
  m = nx0[0];
  n = size[0];
  double fx[m],s_h,z,h, tmp;

  h = bw[0];
  s_h = pow(sig[0]/h,2.0);

  for(i=0;i<m;i++){
    fx[i] = 0.0;
    Sn[i] = 0.0;
  }
  
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      z = (x0[i]-x[j])/h;
      tmp = lkernel(z,s_h);
      fx[i] += tmp;
      Sn[i] += tmp * tmp;
    }
  }
  for(i=0;i<m;i++){
    x0[i] = fx[i] /( n*h );
    Sn[i] = sqrt(Sn[i]/n)/h; 
    if(x0[i] < 0.0){
      x0[i] = 0.0;
      Sn[i] = 999.;
    }
  }
}

/* Non-parametric regression (Nadaraya-Watson type estimator) with
heteroscedastic Laplacian measurement errors.

To reduce computational burdens, we compute the commonly used
variables first, then compute the integrals based on these
intermediate results. Only the Real part of the Eular formula is
computed in the kernel function. */


void nprHLapSAVE(double *x0,  int *nx0, 
	     double *x, double *y, double *sig, int *size,
	     double *bw, double *gcv)
{
  Fun3p f[1];
  f[0] = KernHLap;
  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double mx[m],fx[m],z,h;
  double tmp;
  
  h = bw[0];

  for(i=0;i<m;i++){
    mx[i] = 0.0; //
    fx[i] = 0.0; //
  }
  
  /* nt[] and pt[] will be used as gloabl variable within this
     function. */
  int k, K=NOPZ16;  //K is changable (Integral & here)
  double nt, pt, dum, psum, nsum;
  double nsig[K][n], psig[K][n];
  double nsig2[K], psig2[K]; //to be pass co compute integrals  
  const double *pB = &B16[K - 1];  //K is changable 
  //  const double *pA = &A16[K - 1];  //K is changable

  k=K-1;
  for (; pB >= B16; pB--,k--) {
    dum = 0.5 * *pB;

    nt = 0.5 - dum; 
    pt = 0.5 + dum;

    nsum = 0.0;psum=0.0;//sum of exp part;
    for(j=0;j<n;j++){
      //if(fabs(nt) < 1.414*h/sig[i]){
      //CF of h.lap(0,sig_j) -> t=t/h
      nsig[k][j] = 1.0/(1.0 + 0.5 * pow(sig[i]*nt/h,2.0));
      //}else{
      //nsig[k][j] = 0.0;
      //}
      //if(fabs(pt) < 1.414*h/sig[i]){
      psig[k][j] = 1.0/(1.0 + 0.5 * pow(sig[i]*pt/h,2.0));
      //}else{
      //psig[k][j] = 0.0;
      //}
      nsum += pow(nsig[k][j],2.0);//sum squares of varphi[Uj(t/h)]
      psum += pow(psig[k][j],2.0);
    }
    for(j=0;j<n;j++){
      nsig[k][j] = n * nsig[k][j]/nsum;// 1/psi[Uj(t/h)]
      psig[k][j] = n * psig[k][j]/psum;
    }
  }

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      z = (x0[i]-x[j])/h;
      for(k=0; k<K; k++){
	psig2[k] = psig[k][j];
	nsig2[k] = nsig[k][j];
      }
      tmp = GLInt3d(-10.,10.,f[0],z, psig2, nsig2); //deconvolution kernel Lj
      if(tmp > 0.0){
	fx[i] += tmp;
	mx[i] += tmp * y[i];
      }
    }
  }

  // reuse 'x0' to export the results of mH(x)
  for(i=0;i<m;i++){
    if(fx[i] <= 0.0){
      x0[i] = 0.0;
    }else{
      x0[i] = mx[i] / fx[i]; // return mH(x)
    }
    //printf("ID=%d, m(x)=%f10.2\n",i, x0[i]);
  }
}

double subhlap(double t,double z,double h,double *sig,
		int j,int n){
  double res=0.0,psi=0.0,sum,th,tmp;
  int i;

  th = t/h;
  psi = 1.0/(1.0 + 0.5 * pow(sig[j] * th, 2.0)); //symmetric: -t vs t

  sum = 0.0;
  for(i=0;i<n;i++){
    tmp = 1.0/(1.0 +  0.5 * pow(sig[i] * th,2.0));
    sum += pow(tmp, 2.0);
  }
  psi = sum /(n*psi);

  res = cos(-t * z) * exp(-0.5*t*t) / psi;
  return res;
}

double GLInt5p(double a, double b, 
	       double (*f)(double, double,double,double*,int,int),
	       double z, double h, double *sig, int j,int n)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B16[NOPZ16 - 1];
   const double *pA = &A16[NOPZ16 - 1];

   for (; pB >= B16; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*f)(d - dum,z,h,sig,j,n) + (*f)(d + dum,z,h,sig,j,n) );
   }

   return c * integral;
}

void nprHLap(double *x0,  int *nx0, 
	     double *x, double *y, double *sig, int *size,
	     double *bw, double *gcv)
{
  Fun5p f[1];
  f[0] = subhlap;

  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double mx,fx,z,h,Lj;
  h = bw[0]; 

  for(i=0;i<m;i++){
    fx = 0.0;
    mx = 0.0;
    for(j=0;j<n;j++){
      z = (x0[i] - x[j])/h;
      Lj = GLInt5p(z-10.,z+10.,f[0],z,h,sig,j,n); //deconvolution kernel Lj      
      //if(Lj > 0.0){
      fx += Lj;
      mx += y[j] * Lj;
      //}
    }
    x0[i] = mx/fx;
  }  

  gcv[0] = 0.0;
  for(i=0;i<n;i++){
    fx = 0.0;
    mx = 0.0;
    for(j=0;j<n;j++){
      if(i!=j){
	z = (x0[i] - x[j])/h;
	Lj = GLInt5p(z-10.,z+10.,f[0],z,h,sig,j,n); //deconvolution kernel Lj
	//if(Lj > 0.0){
	fx += Lj;
	mx += y[j] * Lj;
	//}
      }
    }
    gcv[0] += pow(mx/fx-y[i],2.0);
  } 
  gcv[0] /= (double)n;
}

double nwreg(double x, double Z[], double Y[],int n, double h){
  int i;
  double usum=0.,dsum=0.,tmp=0.0;
  for(i=0;i<n;i++){
    tmp = exp(-0.5*pow((x-Z[i])/h, 2.0));
    usum += Y[i] * tmp;
    dsum += tmp;
  }
  return(usum/dsum);
}

double nwreg2(int j, double Z[], double Y[],int n, double h){
  int i;
  double usum=0.,dsum=0.,tmp=0.0,x=Z[j];
  for(i=0;i<n;i++){
    if(i != j){
      tmp = exp(-0.5*pow((x-Z[i])/h, 2.0));
      usum += Y[i] * tmp;
      dsum += tmp;
    }
  }
  return(usum/dsum);
}

void NWReg(double *Z,double *Y,int *size, 
	   double *bandwidth, double *X, int *nx, 
	   double *loo, int *optim, double *opt){
  int i,j,n=size[0];
  double hopt=bandwidth[0];

  // choose optimal bandwidth using the risk as the leave-one-out
  // cross-validation score:
  double h,hstep,rh,t,Rh=1.0e9;
  if(optim[0] > 0){
    h = bandwidth[0]*0.8;
    hstep = 0.0035*bandwidth[0];
    for(j=0;j<400;j++){
      rh = 0.0;
      for(i=0;i<size[0];i++){
	if(loo[0] > 0.0){
	  t = nwreg2(i,Z,Y,n,h);
	}else{
	  t = nwreg(Z[i],Z,Y,n,h);
	}
	rh += pow(t-Y[i],2.0);
      }
      rh = rh / size[0];
      if(rh < Rh){
	Rh = rh;
	hopt = h;
      }
      h += hstep;
    }  
    bandwidth[0] = hopt;
    opt[0] = Rh;
  }else{
    rh = 0.0;
    for(i=0;i<size[0];i++){
      if(loo[0] > 0.0){
	t = nwreg2(i,Z,Y,n,hopt);
      }else{
	t = nwreg(Z[i],Z,Y,n,hopt);
      }
      rh += pow(t-Y[i],2.0);
    }
    rh = rh / size[0];
    opt[0] = rh;
  }
  
  for(i=0;i<nx[0];i++){
    X[i] = nwreg(X[i],Z,Y,n,hopt);
  }
}

void lprHLap(double *x0,  int *nx0, 
	     double *x, double *y, double *sig, int *size,
	     double *bw, double *gcv)
{
  Fun5p f[1];
  f[0] = subhlap;

  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double mx,fx,z,h,Lj;
  h = bw[0]; 

  double bi[n],d,S1,S2;

  for(i=0;i<m;i++){
    S1 = 0.0;
    S2 = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      z = -d/h;
      Lj = GLInt5p(z-10.,z+10.,f[0],z,h,sig,j,n); //deconvolution kernel Lj      
      S1 += Lj*d;
      S2 += Lj * d * d;
      bi[j] = Lj;
    }
    fx = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      bi[j] = bi[j] * (S2-S1*d);
      fx += bi[j];
    }
    mx = 0.0;
    for(j=0;j<n;j++){
      bi[j] /= fx;
      mx += bi[j]*y[j];
    }
    x0[i] = mx;
  }  

  gcv[0] = 0.0;
  for(i=0;i<n;i++){
    S1 = 0.0;
    S2 = 0.0;
    for(j=0;j<n;j++){
      if(i!=j){
	d = x[j] - x[i]; 
	z = -d/h;
	Lj = GLInt5p(z-10.,z+10.,f[0],z,h,sig,j,n); //deconvolution kernel Lj
	S1 += Lj*d;
	S2 += Lj * d * d;
	bi[j] = Lj;
      }
    }
    fx = 0.0;
    for(j=0;j<n;j++){
      if(i!=j){
	d = x[j] - x[i]; 
	bi[j] = bi[j] * (S2-S1*d);
	fx += bi[j];
      }
    }
    mx = 0.0;
    for(j=0;j<n;j++){
      if(i!=j){
	bi[j] /= fx;
	mx += bi[j]*y[j];
      }
    }
    gcv[0] += pow(mx-y[i],2.0);
  } 
  gcv[0] /= (double)n;
}


void lprLap(double *x0,  int *nx0, 
	     double *x, double *y, double *sig, int *size,
	     double *bw, double *gcv)
{
  // Fun5p f[1];
  // f[0] = subhlap;

  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double mx,fx,z,h,Lj;
  h = bw[0]; 

  double bi[n],d,S1,S2;

  for(i=0;i<m;i++){
    S1 = 0.0;
    S2 = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      z = -d/h;
      if(sig[0]==0.0){
	Lj = dnorm(z,0.0,1.0,0);
      }else{
	Lj = KernGL(z,sig[0]);
      }
      S1 += Lj*d;
      S2 += Lj * d * d;
      bi[j] = Lj;
    }
    fx = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      bi[j] = bi[j] * (S2-S1*d);
      fx += bi[j];
    }
    mx = 0.0;
    for(j=0;j<n;j++){
      bi[j] /= fx;
      mx += bi[j]*y[j];
    }
    x0[i] = mx;
  }  

  gcv[0] = 0.0;
  for(i=0;i<n;i++){
    S1 = 0.0;
    S2 = 0.0;
    for(j=0;j<n;j++){
      if(i!=j){
	d = x[j] - x[i]; 
	z = -d/h;
	if(sig[0]==0.0){
	  Lj = dnorm(z,0.0,1.0,0);
	    }else{
	  Lj = KernGL(z,sig[0]);
	}
	S1 += Lj*d;
	S2 += Lj * d * d;
	bi[j] = Lj;
      }
    }
    fx = 0.0;
    for(j=0;j<n;j++){
      if(i!=j){
	d = x[j] - x[i]; 
	bi[j] = bi[j] * (S2-S1*d);
	fx += bi[j];
      }
    }
    mx = 0.0;
    for(j=0;j<n;j++){
      if(i!=j){
	bi[j] /= fx;
	mx += bi[j]*y[j];
      }
    }
    gcv[0] += pow(mx-y[i],2.0);
  } 
  gcv[0] /= (double)n;
}

double hwlprNorm(double x[], double y[], double w[],
		 int n, double h0){
  int i,j,k;
  double mx,fx,z,h,hstep,Lj,hopt=h0;
  double bi[n],d,S1,S2,gmin=9.e+10,gcv;
  //search over [h/2,2h] with ngrid=30
  hstep = 0.25*h0; //(10.3-0.3)*h/500
  h = 0.3*h0;
  for(k=0;k<500;k++){
    h += hstep; 
    gcv = 0.0;
    for(i=0;i<n;i++){
      S1 = 0.0;
      S2 = 0.0;
      for(j=0;j<n;j++){
	if(i!=j){
	  d = x[j] - x[i]; 
	  z = -d/h;
	  Lj = dnorm(z,0.0,1.0,0);
	  S1 += Lj*d;
	  S2 += Lj * d * d;
	  bi[j] = Lj;
	}
      }
      fx = 0.0;
      for(j=0;j<n;j++){
	if(i!=j){
	  d = x[j] - x[i]; 
	  bi[j] = bi[j] * (S2-S1*d);
	  fx += bi[j]*w[j];
	}
      }
      mx = 0.0;
      for(j=0;j<n;j++){
	if(i!=j){
	  bi[j] /= fx;
	  mx += bi[j]*y[j]*w[j];
	}
      }
      gcv += pow(mx-y[i],2.0);
    } 
    if(gcv < gmin){
      gmin = gcv;
      hopt = h;
    }
    //printf("LSCV = %f10.2, bw = %f10.2, hopt=%f10.2\n", gcv , h, hopt);
  }
  return(hopt);
}


void wlprNorm(double *x0,  int *nx0, 
	      double *x, double *y, double *w, int *size,
	      double *bw, int *opt)
{
  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double mx,fx,z,h,Lj;
  h = bw[0]; 
  
  if(opt[0]>0){//find optimal bandwidth
    h = hwlprNorm(x,y,w,n,h);
    bw[0] = h; //save bw
  }

  double bi[n],d,S1,S2;

  for(i=0;i<m;i++){
    S1 = 0.0;
    S2 = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      z = -d/h;
      Lj = dnorm(z,0.0,1.0,0);
      S1 += Lj*d;
      S2 += Lj * d * d;
      bi[j] = Lj;
    }
    fx = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      bi[j] = bi[j] * (S2-S1*d);
      fx += bi[j]*w[j];
    }
    mx = 0.0;
    for(j=0;j<n;j++){
      bi[j] /= fx;
      mx += bi[j]*y[j]*w[j];
    }
    x0[i] = mx;
  }  
}

double ahwlprNorm(double x[], double y[], double w[],
		  int n, double bw[], double h0){
  int i,j,k;
  double mx,fx,z,h,hstep,Lj,hopt=h0;
  double bi[n],d,S1,S2,gmin=9.e+10,gcv;
  //search over [h/2,2h] with ngrid=30
  hstep = 0.35*h0; //(10.3-0.3)*h/500
  h = 0.3*h0;
  for(k=0;k<500;k++){
    h += hstep; 
    gcv = 0.0;
    for(i=0;i<n;i++){
      S1 = 0.0;
      S2 = 0.0;
      for(j=0;j<n;j++){
	if(i!=j){
	  d = x[j] - x[i]; 
	  z = -d/h/bw[j];
	  Lj = dnorm(z,0.0,1.0,0);
	  S1 += Lj*d;
	  S2 += Lj * d * d;
	  bi[j] = Lj;
	}
      }
      fx = 0.0;
      for(j=0;j<n;j++){
	if(i!=j){
	  d = x[j] - x[i]; 
	  bi[j] = bi[j] * (S2-S1*d);
	  fx += bi[j]*w[j];
	}
      }
      mx = 0.0;
      for(j=0;j<n;j++){
	if(i!=j){
	  bi[j] /= fx;
	  mx += bi[j]*y[j]*w[j];
	}
      }
      gcv += pow(mx-y[i],2.0);
    } 
    if(gcv < gmin){
      gmin = gcv;
      hopt = h;
    }
    //printf("LSCV = %f10.2, bw = %f10.2, hopt=%f10.2\n", gcv , h, hopt);
  }
  return(hopt);
}

void awlprNorm(double *x0,  int *nx0, 
	      double *x, double *y, double *w, int *size,
	      double *bandwidth)
{
  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double d,mx,fx,z,h,Lj,bw[n];
  h = bandwidth[0];

  for(i=0;i<n;i++){
    fx = 0.0;
    bw[i] = 0.0;
    for(j=0;j<n;j++){
      d = (x[i] - x[j])/h; 
      bw[i] += dnorm(d,0.0,1.0,0);
    }
    bw[i] = 1.0/bw[i];
    fx += bw[i];
  }

  for(i=0;i<n;i++){
    bw[i] *= h/fx; 
  }
  
  h = ahwlprNorm(x,y,w,n,bw,h);

  double bi[n],S1,S2;

  for(i=0;i<m;i++){
    S1 = 0.0;
    S2 = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      z = -d/h/bw[j];
      Lj = dnorm(z,0.0,1.0,0);
      S1 += Lj*d;
      S2 += Lj * d * d;
      bi[j] = Lj;
    }
    fx = 0.0;
    for(j=0;j<n;j++){
      d = x[j] - x0[i]; 
      bi[j] = bi[j] * (S2-S1*d);
      fx += bi[j]*w[j];
    }
    mx = 0.0;
    for(j=0;j<n;j++){
      bi[j] /= fx;
      mx += bi[j]*y[j]*w[j];
    }
    x0[i] = mx;
  }  
}


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
	       double y[], int ny, double fy[])
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
      if(x[j]<4.*h){
	fy[i] += w[j]*(dnorm(y[i]-x[j],0.,h*hs[j],0) + 
		       dnorm(y[i]+x[j],0.,h*hs[j],0));
      }else{
	fy[i] += w[j]*dnorm(y[i]-x[j],0.,h*hs[j],0);
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
  if((sp<0.)|(sp>1.)){
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
  awpdf(x,n,w,h,hs,x,n,fx);
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
  awpdf(x,n,w,h,hs,x,n,fx);
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
	   int *ysize, double *pars)
{
  int i,nx=xsize[0],ny=ysize[0],npar=2;
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
    awpdf(x,nx,w,h,hs,x,nx,fx);
    UpdateBwfactor(fx,nx,sp,hs);
    pars[0]=h; pars[1]=sp;
  }
  
  awpdf(x,nx,w,h,hs,y,ny,fy);
  awcdf(x,nx,w,h,hs,y,ny,Fy);
}

void wkde(double *x,double *w,int *xsize,
	   double *y, double *fy, double *Fy, 
	   int *ysize, double *pars)
{
  int i,nx=xsize[0],ny=ysize[0],npar=2;
  double h0,h,hs[nx],sp=0.8; //sp=sensitivity parameter;
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
  
  awpdf(x,nx,w,h,hs,y,ny,fy);
  awcdf(x,nx,w,h,hs,y,ny,Fy);
}

void ckde(double *x,double *w,int *xsize,
	   double *y, double *fy, double *Fy, 
	   int *ysize, double *pars)
{
  int i,nx=xsize[0],ny=ysize[0];
  double h,hs[nx];

  h = pars[0]; 
  for(i=0;i<nx;i++){hs[i]=1.;}
  
  awpdf(x,nx,w,h,hs,y,ny,fy);
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


void subdKDE(double y0[], double x0[], int n0, 
	  double x[], double h[], double f[], int n)
{
  int i,j;
  double xi,nsum=0.0;
  for(j=0; j<n;j++){
    nsum += f[j];
  }

  for(i=0; i<n0;i++){
    y0[i] = 0.0;
    for(j=0; j<n;j++){
      xi = (x0[i] - x[j])/h[j];
      y0[i] += dnorm(xi,0.,1.,0)/h[j]*f[j];
    }
    y0[i] /= nsum;
  }
}

void subpKDE(double y0[], double x0[], int n0, 
	  double x[], double h[], double f[], int n)
{
  int i,j;
  double xi,nsum=0.0;
  for(j=0; j<n;j++){
    nsum += f[j];
  }

  for(i=0; i<n0;i++){
    y0[i] = 0.0;
    for(j=0; j<n;j++){
      xi = (x0[i] - x[j])/h[j];
      y0[i] += pnorm(xi,0.,1.,1,0)/h[j]*f[j]; 
    }
    y0[i] /= nsum;
  }
}

void dKDE(double *x0, int *n0, 
	  double *x, double *h, double *f, int *n)
{
  int i,k=n0[0];
  double y[k];
  for(i=0; i<k;i++){
    y[i] = 0.0;
  }

  subdKDE(y,x0,k,x,h,f,n[0]);

  for(i=0; i<k;i++){
    x0[i] = y[i];
  }
}

void pKDE(double *x0, int *n0, 
	  double *x, double *h, double *f, int *n)
{
  int i,k=n0[0];
  double y[k];
  for(i=0; i<k;i++){
    y[i] = 0.0;
  }

  subpKDE(y,x0,k,x,h,f,n[0]);

  for(i=0; i<k;i++){
    x0[i] = y[i];
  }
}
