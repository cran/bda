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

typedef double (*Fun2d)(double,double);
typedef double (*Fun3d)(double,double,double);

typedef double (*Fun1p)(double);
typedef double (*Fun3p)(double,double,double);
//typedef double (*Fun4p)(double,double*,double,int);
//typedef double (*Fun5p)(double,double,double,double*,int);
//typedef double (*Fun6p)(double,double,double,double*,double*,int);
typedef double (*Fvvvi)(double,double*,double*,double*,int);
typedef double (*Fvdi)(double,double*,double,int);
typedef double (*Fun7p)(double,double,double,double*,double*,double*,int);


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

double GLIntvdi(double a, double b, 
		double (*f)(double, double*,double,int),
		double *x,double h,int n)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
     dum = c * *pB;
     integral += *pA * ( (*f)(d - dum,x,h,n) + (*f)(d + dum,x,h,n) );
   }
   
   return c * integral;
}

double KGauss(double z){
  return dnorm(z,0.0,1.0,0);
}

void lpreg(double xgrid[], int m, double x[], double y[],
	   int n, double h, double rhat[], double ellx[])
{
  int i,j;
  double K[n], l[n];
  double S1, S2, t1, t2, bsum;
  
  for(i=0; i<m; i++){
    // to compute K, S1, and S2
    S1 = 0.0;
    S2 = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - xgrid[i];
      K[j] = KGauss(t1/h);
      t2 = K[j] * t1;
      S1 += t2;
      S2 += t2 * t1;
    }
    // compute li: use li[i][j] to store b[i][j]
    bsum = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - xgrid[i];
      l[j] = K[j]*(S2-S1*t1); //bi(x)
      bsum += l[j];
    }
    rhat[i] = 0.0;
    ellx[i] = 0.0;
    for(j=0; j<n; j++){
      l[j] /= bsum;
      rhat[i] += l[j]*y[j];
      ellx[i] += l[j]*l[j];
    }
    ellx[i] = sqrt(ellx[i]);
  }
}

void alpreg(double xgrid[], int m, double x[], double y[],
	    int n, double h, double alpha, double rhat[], 
	    double ellx[])
{
  int i,j;
  double K[n], l[n], fhat[n], g, ah;
  double S1, S2, t1, t2, bsum;
  
  lpreg(x,n,x,y,n,h,fhat,ellx);
  
  t1 = 0.0;  //reuse
  for(i=0;i<n;i++){
    if(fhat[i]>0) t1 += log(fhat[i]);
  }
  g = exp(t1/n);

  for(i=0; i<m; i++){
    // to compute K, S1, and S2
    S1 = 0.0;
    S2 = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - xgrid[i];
      ah = h*pow(fhat[j]/g, -alpha);
      K[j] = KGauss(t1/ah);
      t2 = K[j] * t1;
      S1 += t2;
      S2 += t2 * t1;
    }
    // compute li: use li[i][j] to store b[i][j]
    bsum = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - xgrid[i];
      l[j] = K[j]*(S2-S1*t1); //bi(x)
      bsum += l[j];
    }
    rhat[i] = 0.0;
    ellx[i] = 0.0;
    for(j=0; j<n; j++){
      l[j] /= bsum;
      rhat[i] += l[j]*y[j];
      ellx[i] += l[j]*l[j];
    }
    ellx[i] = sqrt(ellx[i]);
  }
}

double lscvscore(double x[], double y[], int n, double h){
  int i,j;
  double K[n], l[n], rhat[n];
  double nu=0.0,S1, S2, t1, t2, bsum;


  for(i=0; i<n; i++){
    // to compute S1, and S2
    S1 = 0.0;
    S2 = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - x[i];
      K[j] = KGauss(t1/h);
      t2 = K[j] * t1;
      S1 += t2;
      S2 += t2 * t1;
    }
    // compute li: use li[i][j] to store b[i][j]
    bsum = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - x[i];
      l[j] = K[j]*(S2-S1*t1); //bi(x)
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

double alscvscore(double x[], double y[], int n, double h, 
		  double alpha, double fhat[]){
  int i,j;
  double K[n], l[n], rhat[n], ah, g;
  double nu=0.0, S1, S2, t1, t2, bsum;

  t1 = 0.0;  //reuse
  for(i=0;i<n;i++){
    if(fhat[i]>0) t1 += log(fhat[i]);
  }
  g = exp(t1/n);

  for(i=0; i<n; i++){
    // to compute S1, and S2
    S1 = 0.0;
    S2 = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - x[i];
      ah = h*pow(fhat[j]/g, -alpha);
      K[j] = KGauss(t1/ah);
      t2 = K[j] * t1;
      S1 += t2;
      S2 += t2 * t1;
    }
    // compute li: use li[i][j] to store b[i][j]
    bsum = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - x[i];
      l[j] = K[j]*(S2-S1*t1); //bi(x)
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

double lprlscv(double x[], double y[], int n, double h)
{  
  int i;
  double delta=0.03*h, h0, hopt, Rhmin=1.0e7, Rhat=0.0;
  h0 = 0.3 * h; hopt = h0;
  for(i=0; i<101; i++){
    Rhat = lscvscore(x,y,n,h0);
    if(Rhat <= Rhmin && R_FINITE(Rhat)){
      hopt = h0;
      Rhmin = Rhat;
      //      Rprintf("Rhat=: %f\n",Rhat);
    }
    h0 += delta;
  }  
  return hopt;
}

void alprlscv(double x[], double y[], int n, double h, double alpha)
{  
  int i,j;
  double delta=0.1*h, h0, hopt, Rhmin=1.0e7, Rhat=0.0;
  double rhat[n],ellx[n]; //ellx[] is a pseudo variable
  double alpha0, dalpha;

  alpha = 0.0;
  for(i=0; i<n;i++){
    rhat[i] = 0.0;
    ellx[i] = 0.0;
  }
  lpreg(x,n,x,y,n,h,rhat,ellx);

  alpha0 = 0.0; dalpha = 0.1;
  h0 = 0.3 * h; hopt = h0;
  for(i=0; i<21; i++){
    for(j=0; j<11; j++){
      Rhat = alscvscore(x,y,n,h0,alpha0,rhat);
      if(Rhat <= Rhmin && R_FINITE(Rhat)){
	hopt = h0;
	alpha = alpha0;
	Rhmin = Rhat;
      }
      alpha0 += dalpha;
    }
    h0 += delta;
  }  
  h = hopt;
}

double funGauss(double t,double *x, double h, int n){
  int i;
  double Sn1,Sn2,deSn1,deSn2,dx,z,Sbi,Sdbi;
  double bi[n],dbi[n];
  double t0;
  double SLl,Sll,SLL;
  Sn1 = 0.0; Sn2 = 0.0; deSn1 = 0.0; deSn2 = 0.0;
  
  for(i=0;i<n;i++){
    dx = t-x[i];
    z = dx/h;
    bi[i] = KGauss(z); //store K to compute bi later on
    dbi[i] = z * KGauss(z)/h;//store K' to compute bi later on
    t0 = bi[i] * dx;
    Sn1 += t0;
    Sn2 += t0 * dx;
    deSn1 += dbi[i] * dx + bi[i]; 
    deSn2 += dbi[i] * dx * dx + 2.0 * t0;
  }

  Sbi = 0.0; Sdbi = 0.0;
  for(i=0;i<n;i++){
    dx = t-x[i];
    dbi[i] = dbi[i]*(Sn2-Sn1*dx)+ bi[i]*(deSn2-Sn1-deSn1*dx);
    bi[i] *= Sn2-Sn1*dx;
    Sbi += bi[i];
    Sdbi += dbi[i];
  }
  SLl=0.0; Sll=0.0;SLL=0.0;
  for(i=0;i<n;i++){
    dx = t-x[i];
    dbi[i] = dbi[i]/Sbi-bi[i]*Sdbi/Sbi;//reuse dbi for dli
    bi[i] = bi[i] / Sbi;//reuse bi for li
    SLL += bi[i] * bi[i];
    SLl += bi[i] * dbi[i];
    Sll += dbi[i] * dbi[i];
  }
  z = Sll/SLL - SLl*SLl/SLL/SLL; 
  if((z<0)|(ISNAN(z))) z=0.000001;
  return sqrt(z);
}

void lpsmooth(double *xgrid, int *ngrid, double *x, double *y, 
	      int *size, double *bw, int *lscv, double *range, 
	      int *adapt, double *ellx, double *kappa)
{
  int i;
  double rhat[ngrid[0]], rhat2[size[0]],hopt=bw[0],alpha=0.0;

  for(i=0; i<ngrid[0]; i++){
    rhat[i] = 0.0;
  }
  for(i=0; i<size[0]; i++){
    rhat2[i] = 0.0;
  }

  if(lscv[0] == 1){
    if(adapt[0]==0){
      hopt = lprlscv(x,y,size[0],hopt);
      lpreg(xgrid, ngrid[0],x,y,size[0],hopt,rhat,ellx);
      lpreg(x, size[0],x,y,size[0],hopt,rhat2,ellx);
    }else{
      alprlscv(x,y,size[0],hopt, alpha);
      alpreg(xgrid, ngrid[0],x,y,size[0],hopt,alpha,rhat,ellx);
      alpreg(x, size[0],x,y,size[0],hopt,alpha,rhat2,ellx);
    }
  }else{
    lpreg(xgrid, ngrid[0],x,y,size[0],hopt,rhat,ellx);
    lpreg(x, size[0],x,y,size[0],hopt,rhat2,ellx);
  }
  
  bw[0] = hopt;
  for(i=0; i<ngrid[0]; i++){
    xgrid[i] = rhat[i];
  }
  for(i=0; i<size[0]; i++){
    y[i] = rhat2[i];
  }

  Fvdi f[1];
  f[0] = funGauss;
  //compute kappa0
  kappa[0] = GLIntvdi(range[0],range[1],f[0],x,hopt,size[0]); 
}

void tubecv(double *kappa, double *level){
  /* compute the critical value of the simultaneous confidence band
     using the tube formula.  Input: kappa (computed in LLS) */
  int i=0;
  double dx=10.0, z=2.0;//initial value
  
  while((i<100)&(fabs(dx)>0.000001)){
    dx = (2.*(1-pnorm(z,0.,1.,1,0))+kappa[0]/M_PI*exp(-.5*z*z)
	  -1.0 + level[0])/
      (2.*dnorm(z,0.,1.,0)+z*kappa[0]/M_PI*exp(-0.5*z*z));
    z += dx;
    i++;
  }
  if(i>=100) z=-999.;
  kappa[0] = z;
}

double GLInt7p(double a, double b,  
	       double (*fn)(double,double,double,double*,double*,double*,int),
	       double h, double g,double *half, double *w,double *f,int n)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*fn)(d - dum,h,g,half,w,f,n) + 
			  (*fn)(d + dum,h,g,half,w,f,n) );
   }

   return c * integral;
}

/*
The following codes are developed to find the bootstrap MISE optimal
bandwidth.  The common factor of 1/(2PIn^2) was disregarded.
 */

double fa(double t, double h, double g, double *H, double *W, double *f,int n){
  int j;
  double h2t2 = pow(h*t,2.0),g2t2=pow(g*h,2.0),tmp,sum1,sum2,res;
  sum1=0.;sum2=0.;tmp=0.;
  tmp = (1.-1./n)*exp(-(g2t2+h2t2))-2.*exp(-(.5*h2t2+g2t2))+exp(-g2t2);
  if(fabs(t)==0.0){
    for(j=0;j<n;j++){
      sum1 += cos(t*W[j])*f[j];
      sum2 += sin(t*W[j])*f[j];
    }
  }else{
    for(j=0;j<n;j++){
      res = sin(H[j]*t)/H[j]/t;
      sum1 +=  res * cos(t*W[j])*f[j];
      sum2 +=  res * sin(t*W[j])*f[j];
    }
  }
  res= tmp * (sum1*sum1+sum2*sum2);
  return res;
}

void hbmise(double *x,double *f, double *w, int *size, double *hopt)
/* Input:
   (1) x[j]: the observations;
   (3) w[j]: = b[j]-a[j], which the support of the uniform error;
   (4) size, or n in the code: an integer shows the number of distinct values in y[.]
   (5) f: frequencies
   Output: 
   (1) hopt: the initial bandwidth, which will also be the output.
     
   Derived parameters:
   (1) h is the bandwidth, which is searched over a defined grid by
   minimizing the BMISE.  We set the grid number 100.  The initial
   bandwidth h0 was selected by using Silverman's estimate
   (h=1.06*sigma/n^.2) or the dpik or dpih in package KernSmooth,
   and we search over the range from h0/10 to k*h0.  Where
   k=1+max(bj-aj)^2/12/sd(y). We can simply select 2. If the
   bandwidth is on the boundary, expand the search.

   Controls
   (1) The integration range is set to t \in [0,U], where U>4/h=40/h0.
   (2) iter is used to control the maximum search loop.  Maximum 10.

   Remarks: check whether erf() is a build-in function.  If yes, use
   it directly.

   */
{
  int ngrid = 50, n=size[0],i,j,nsum=0;
  double h0,h,hstep,h2;
  double tmp1,tmp2,mise, fint,fb;
  double w1_2[n],w1_2sq[n];
  double mise0=9999999999.;
  Fun7p fn[1];
  fn[0] = fa; //To compute part A to be integrated.

  tmp1 = 0.0;
  tmp2 = 0.0;
  for(i=0;i<n;i++){
    nsum += f[i];
    tmp1 = tmp1 + x[i]*f[i];
    tmp2 = tmp2 + x[i] * x[i]*f[i];
    w1_2[i] = w[i]/2.0;
    w1_2sq[i]=pow(w1_2[i],2.0);
  }
  // initialize bandwidths
  h0 = hopt[0];
  h = h0/30.;
  h2 = h * h;
  hstep = 1.5*h0/ngrid;

  for(i=0;i<ngrid;i++){
    h += hstep;
    fb = 0.0;
    for(j=0;j<n;j++){
      tmp1 = exp(-w1_2sq[j]/h2);
      tmp2 = w1_2[j]/h*M_SQRT_PI*erf(w1_2[j]/h);
      fb += fabs(1./w1_2sq[j]*(tmp1+tmp2-1.))*f[j];
    }
    fint = GLInt7p(0.0,1000/h0,fn[0],h,h,w1_2,x,f,n);
    mise=M_SQRT_PI*h*fb;
    tmp1 = (mise+2*fint/nsum/nsum)/2./pow(nsum,2.0)/M_PI;
    if(tmp1<mise0){
      hopt[0] = h;
      mise0 = tmp1;
    }
  }
  //  hopt[0] = h0;
}
 
/*  End of BMISE bandwidth selector */



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


// ISNA(x): true for R's NA only
// ISNAN(x): true for R's NA and IEEE NaN
// R_FINITE(x): false for Inf,-Inf,NA,NaN
// R_IsNaN(x): true for NaN but not NA

// Rprintf:  printing from a C routine compiled into R

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


void DesignMatrix(double *xv, int *size, double *bw, double *dm){
  int i,j, n =size[0];
  double K[n][n],l[n][n];
  double S1[n], S2[n]; // two vectors S1(xj), S2(xj)
  double t1, bsum, h=bw[0];
  
  // to compute the kernel values
  for(i=0;i<n;i++){
    for(j=i; j<n;j++){
      t1 = xv[i] - xv[j];
      K[i][j] = dnorm(t1/h,0.0, 1.0, 0);
      K[j][i] = K[i][j];
    }
  }
  // compute S1, S2.
  for(j=0; j<n;j++){
    S1[j] = 0.0; S2[j] = 0.0;
    for(i=0;i<n;i++){
      t1 = xv[i] - xv[j];
      S1[j] += K[i][j] * t1;
      S2[j] += K[i][j] * t1 * t1;
    }
  }
  // compute B and Lii: store sum(bi) to B, and bi to Lii
  for(j=0; j<n;j++){
    bsum = 0.0;
    for(i=0;i<n;i++){
      t1 = xv[i] - xv[j];
      l[i][j] = K[i][j] * (S2[j]- t1 * S1[j]);
      bsum += l[i][j];
    }
    for(i=0;i<n;i++){
      dm[j*n+i] = l[i][j]/bsum;      
    }
  }
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
