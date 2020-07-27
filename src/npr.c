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
typedef double (*Fun5p)(double,double,double,double*,int);
typedef double (*Fun6p)(double,double,double,double*,double*,int);


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

//  100pts

static const double B100[] = {
    1.56289844215430828714e-02,    4.68716824215916316162e-02,
    7.80685828134366366918e-02,    1.09189203580061115002e-01,
    1.40203137236113973212e-01,    1.71080080538603274883e-01,
    2.01789864095735997236e-01,    2.32302481844973969643e-01,
    2.62588120371503479163e-01,    2.92617188038471964730e-01,
    3.22360343900529151720e-01,    3.51788526372421720979e-01,
    3.80872981624629956772e-01,    4.09585291678301542532e-01,
    4.37897402172031513100e-01,    4.65781649773358042251e-01,
    4.93210789208190933576e-01,    5.20158019881763056670e-01,
    5.46597012065094167460e-01,    5.72501932621381191292e-01,
    5.97847470247178721259e-01,    6.22608860203707771585e-01,
    6.46761908514129279840e-01,    6.70283015603141015784e-01,
    6.93149199355801965946e-01,    7.15338117573056446485e-01,
    7.36828089802020705530e-01,    7.57598118519707176062e-01,
    7.77627909649495475605e-01,    7.96897892390314476375e-01,
    8.15389238339176254384e-01,    8.33083879888400823522e-01,
    8.49964527879591284320e-01,    8.66014688497164623416e-01,
    8.81218679385018415547e-01,    8.95561644970726986709e-01,
    9.09029570982529690453e-01,    9.21609298145333952679e-01,
    9.33288535043079545942e-01,    9.44055870136255977955e-01,
    9.53900782925491742847e-01,    9.62813654255815527284e-01,
    9.70785775763706331929e-01,    9.77809358486918288561e-01,
    9.83877540706057015509e-01,    9.88984395242991747997e-01,
    9.93124937037443459632e-01,    9.96295134733125149166e-01,
    9.98491950639595818382e-01,    9.99713726773441233703e-01
};

static const double A100[] = {
    3.12554234538633569472e-02,    3.12248842548493577326e-02,
    3.11638356962099067834e-02,    3.10723374275665165874e-02,
    3.09504788504909882337e-02,    3.07983790311525904274e-02,
    3.06161865839804484966e-02,    3.04040795264548200160e-02,
    3.01622651051691449196e-02,    2.98909795933328309169e-02,
    2.95904880599126425122e-02,    2.92610841106382766198e-02,
    2.89030896011252031353e-02,    2.85168543223950979908e-02,
    2.81027556591011733175e-02,    2.76611982207923882944e-02,
    2.71926134465768801373e-02,    2.66974591835709626611e-02,
    2.61762192395456763420e-02,    2.56294029102081160751e-02,
    2.50575444815795897034e-02,    2.44612027079570527207e-02,
    2.38409602659682059633e-02,    2.31974231852541216230e-02,
    2.25312202563362727021e-02,    2.18430024162473863146e-02,
    2.11334421125276415432e-02,    2.04032326462094327666e-02,
    1.96530874944353058650e-02,    1.88837396133749045537e-02,
    1.80959407221281166640e-02,    1.72904605683235824399e-02,
    1.64680861761452126430e-02,    1.56296210775460027242e-02,
    1.47758845274413017686e-02,    1.39077107037187726882e-02,
    1.30259478929715422855e-02,    1.21314576629794974079e-02,
    1.12251140231859771176e-02,    1.03078025748689695861e-02,
    9.38041965369445795116e-03,    8.44387146966897140266e-03,
    7.49907325546471157895e-03,    6.54694845084532276405e-03,
    5.58842800386551515727e-03,    4.62445006342211935096e-03,
    3.65596120132637518238e-03,    2.68392537155348241939e-03,
    1.70939265351810523958e-03,    7.34634490505671730396e-04
};

#define NOPZ100  sizeof(B100) / sizeof(double)
#define NOZ100   NOPZ100+NOPZ100


double GLInteg(double a, double b, double (*f)(double))
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
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
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
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
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];
   int k=100;

   for (; pB >= B100; pA--, pB--,k--) {
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
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
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
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
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
  int K=NOPZ100;  //K is changable
  double ntexp[K],ptexp[K];  //K is changable
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable
  double ppart1,npart1;

   switch(type[0]){
   case 0:
     sb2 = pow(sig[0]/bw[0],2.0)*0.5;
     k=K-1;
     for (; pB >= B100; pA--, pB--,k--) {
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
       pB = &B100[K - 1];   
       pA = &A100[K - 1];  
       for (; pB >= B100; pA--, pB--,k--) {
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
     for (; pB >= B100; pA--, pB--,k--) {
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
       pB = &B100[K - 1];   
       pA = &A100[K - 1];  
       for (; pB >= B100; pA--, pB--,k--) {
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
  int K=NOPZ100;  //K is changable
  double nt[K], pt[K], nt2[K], pt2[K],nt3[K],pt3[K];
  double sigh[ny], z[ny];
  double nsig[K][ny], psig[K][ny];
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable
  double nsum, psum;

  nsum = bw*bw;//reuse valuable
  for(j=0;j<ny;j++){
    sigh[j] = -0.5 * pow(sig[j],2.0)/nsum;
  }
  k=K-1;
  for (; pB >= B100; pA--, pB--,k--) {
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
    pB = &B100[K - 1];   
    pA = &A100[K - 1];  
    for (; pB >= B100; pA--, pB--,k--) {
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
  int K=NOPZ100;  //K is changable
  double nt[K], pt[K], nt2[K], pt2[K],nt3[K],pt3[K];
  double sigh[ny[0]], z[ny[0]];
  double nsig[K][ny[0]], psig[K][ny[0]];
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable
  double nsum, psum;

  nsum = bw[0]*bw[0];//reuse valuable
  for(j=0;j<ny[0];j++){
    sigh[j] = -0.5 * pow(sig[j],2.0)/nsum;
  }
  k=K-1;
  for (; pB >= B100; pA--, pB--,k--) {
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
       pB = &B100[K - 1];   
       pA = &A100[K - 1];  
       for (; pB >= B100; pA--, pB--,k--) {
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
       pB = &B100[K - 1];   
       pA = &A100[K - 1];  
       for (; pB >= B100; pA--, pB--,k--) {
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
    mise = 0.5/n[0]/pow(PI*hdiff,-0.5)+Rfx[0]*pow(h,4.0);
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
  int K=NOPZ100;  //K is changable
  double ntexp[K],ptexp[K];  //K is changable
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable
  double ppart1,npart1,denom,tmp0,ppart2,npart2;
  sb2 = pow(sig[0]/bw[0],2.0)*0.5;
  k=K-1;
  for (; pB >= B100; pA--, pB--,k--) {
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
    pB = &B100[K - 1];   
    pA = &A100[K - 1];  
    for (; pB >= B100; pA--, pB--,k--) {
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
  int K=NOPZ100;  //K is changable
  double nt[K], pt[K], nt2[K], pt2[K],nt3[K],pt3[K];
  double sigh[ny[0]], z[ny[0]];
  double nsig[K][ny[0]], psig[K][ny[0]];
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable
  double nsum, psum,denom,tmp0,nsum2,psum2;

  nsum = bw[0]*bw[0];//reuse valuable
  for(j=0;j<ny[0];j++){
    sigh[j] = -0.5 * pow(sig[j],2.0)/nsum;
  }
  k=K-1;
  for (; pB >= B100; pA--, pB--,k--) {
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
    pB = &B100[K - 1];   
    pA = &A100[K - 1];  
    for (; pB >= B100; pA--, pB--,k--) {
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
      mise=fint/PI;
      if(mise<mmin) {
	hopt=h; mmin=mise;
      }
    }
    break;
  case 2:
    for(i=0;i<grid[0];i++){
      h += hstep;
      fint = GLInt6p(0,4.0,f[2],h,g,sig,y,n);//approximate 
      mise=fint/PI;
      if(mise<mmin) {
	hopt=h; mmin=mise;
      }
    }
    break;
  case 3:
    for(i=0;i<grid[0];i++){
      h += hstep;
      fint = GLInt6p(-1.0,1.0,f[1],h,g,sig,y,n);
      mise=fint/PI;
      if(mise<mmin) {
	hopt=h; mmin=mise;
      }
    }
    break;
  case 4:
    for(i=0;i<grid[0];i++){
      h += hstep;
      fint = GLInt6p(0,4.0,f[3],h,g,sig,y,n);//approximate 
      mise=fint/PI;
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
  int K=NOPZ100;  //K is changable
  double ntexp[K],ptexp[K];  //K is changable
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable
  double ppart1,npart1;

  //Fun4p f[1];
  //f[0] = funSuppNorm2;

  switch(ktype[0]){
  case 0:
    sb2 = pow(sig[0]/bw[0],2.0)*0.5;
    k=K-1;
    for (; pB >= B100; pA--, pB--,k--) {
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
      pB = &B100[K - 1];   
      pA = &A100[K - 1];  
      for (; pB >= B100; pA--, pB--,k--) {
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
  int k, K=NOPZ100;  //K is changable (Integral & here)
  double nt, nt2, nt3, pt, pt2, pt3, dum, psum, nsum;
  double nsig[K][n], psig[K][n];
  double nsig2[K], psig2[K]; //to be pass co compute integrals  
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable

  k=K-1;
  for (; pB >= B100; pA--, pB--,k--) {
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
  int k, K=NOPZ100;  //K is changable (Integral & here)
  double nt, nt2, nt3, pt, pt2, pt3, dum, psum, nsum;
  double nsig[K][n], psig[K][n];
  double nsig2[K], psig2[K]; //to be pass co compute integrals  
  const double *pB = &B100[K - 1];  //K is changable 
  const double *pA = &A100[K - 1];  //K is changable

  k=K-1;
  for (; pB >= B100; pA--, pB--,k--) {
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

