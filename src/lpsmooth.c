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
typedef double (*Fvdime)(double,double*,double,int,double);
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

double GLIntvdime(double a, double b, 
		  double (*f)(double, double*,double,int, double),
		  double *x,double h,int n, double s)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
     dum = c * *pB;
     integral += *pA * ( (*f)(d - dum,x,h,n,s) + (*f)(d + dum,x,h,n,s) );
   }
   
   return c * integral;
}

double KGauss(double z){
  return dnorm(z,0.0,1.0,0);
}

double KLaplace(double x, double h, double s){
  double a=1.0,z;
  z = x/h;
  a = 1.0+s*s*(1.0-z*z)/h/h;
  if(a<0.0) a = 0.0;
  return a * dnorm(z,0.0,1.0,0);
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

/* local polynomial regression with errors-in-variables. Laplacian
   errors are assumed. For more measurement error types, please see
   Decon */


double mescore(double x[], double y[], int n, double h, double s){
  int i,j;
  double K[n], l[n], rhat[n];
  double nu=0.0,S1, S2, t1, t2, bsum;


  for(i=0; i<n; i++){
    // to compute S1, and S2
    S1 = 0.0;
    S2 = 0.0;
    for(j=0; j<n; j++){
      t1 = x[j] - x[i];
      K[j] = KLaplace(t1,h,s);
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

double melscv(double x[], double y[], int n, double h, double s)
{  
  int i;
  double delta=0.03*h, h0, hopt, Rhmin=1.0e7, Rhat=0.0;
  h0 = 0.3 * h; hopt = h0;
  for(i=0; i<101; i++){
    Rhat = lscvscore(x,y,n,h0);
    if(Rhat <= Rhmin && R_FINITE(Rhat)){
      hopt = h0;
      Rhmin = Rhat;
    }
    h0 += delta;
  }  
  return hopt;
}

void mereg(double xgrid[], int m, double x[], double y[],
	   int n, double h, double rhat[], double ellx[],
	   double s)
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
      K[j] = KLaplace(t1,h,s);
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

 void amereg(double xgrid[], int m, double x[], double y[],
	    int n, double h, double alpha, double rhat[], 
	     double ellx[], double s)
{
  int i,j;
  double K[n], l[n], fhat[n], g, ah;
  double S1, S2, t1, t2, bsum;
  
  mereg(x,n,x,y,n,h,fhat,ellx,s);
  
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
      K[j] = KLaplace(t1,ah,s);
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

 double amescore(double x[], double y[], int n, double h, 
		 double alpha, double fhat[], double s){
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
      K[j] = KLaplace(t1,ah,s);
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

 void amelscv(double x[], double y[], int n, double h,
	      double alpha, double s)
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
  mereg(x,n,x,y,n,h,rhat,ellx,s);

  alpha0 = 0.0; dalpha = 0.1;
  h0 = 0.3 * h; hopt = h0;
  for(i=0; i<21; i++){
    for(j=0; j<11; j++){
      Rhat = amescore(x,y,n,h0,alpha0,rhat,s);
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

double funLaplace(double t,double *x, double h, int n, double s){
  int i;
  double Sn1,Sn2,deSn1,deSn2,dx,z,Sbi,Sdbi;
  double bi[n],dbi[n];
  double t0;
  double SLl,Sll,SLL;
  Sn1 = 0.0; Sn2 = 0.0; deSn1 = 0.0; deSn2 = 0.0;
  
  for(i=0;i<n;i++){
    dx = t-x[i];
    z = dx/h;
    bi[i] = KLaplace(dx,h,s); //store K to compute bi later on
    dbi[i] = z * bi[i]/h;//store K' to compute bi later on
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

void llrme(double *xgrid, int *ngrid, double *x, double *y, 
	   int *size, double *bw, int *lscv, double *range, 
	   int *adapt, double *ellx, double *kappa)
{
  int i;
  double rhat[ngrid[0]], rhat2[size[0]],hopt=bw[0],alpha=0.0,sdme=kappa[0];
  
  for(i=0; i<ngrid[0]; i++){
    rhat[i] = 0.0;
  }
  for(i=0; i<size[0]; i++){
    rhat2[i] = 0.0;
  }
  
  if(lscv[0] == 1){
    if(adapt[0]==0){
      hopt = melscv(x,y,size[0],hopt,sdme); //done
      mereg(xgrid, ngrid[0],x,y,size[0],hopt,rhat,ellx,sdme); //done
      mereg(x, size[0],x,y,size[0],hopt,rhat2,ellx,sdme);//done
    }else{
      amelscv(x,y,size[0],hopt, alpha,sdme);//done
      amereg(xgrid, ngrid[0],x,y,size[0],hopt,alpha,rhat,ellx,sdme);//done
      amereg(x, size[0],x,y,size[0],hopt,alpha,rhat2,ellx,sdme);//done
    }
  }else{
    mereg(xgrid, ngrid[0],x,y,size[0],hopt,rhat,ellx,sdme);//done
    mereg(x, size[0],x,y,size[0],hopt,rhat2,ellx,sdme);//done
  }
  
  bw[0] = hopt;
  for(i=0; i<ngrid[0]; i++){
    xgrid[i] = rhat[i];
  }
  for(i=0; i<size[0]; i++){
    y[i] = rhat2[i];
  }

  /*Fvdime f[1];
  f[0] = funLaplace;
  //compute kappa0
  kappa[0] = GLIntvdime(range[0],range[1],f[0],x,hopt,size[0],sdme); */
  Fvdi f[1];
  f[0] = funGauss;
  //compute kappa0
  kappa[0] = GLIntvdi(range[0],range[1],f[0],x,hopt,size[0]); 

}

