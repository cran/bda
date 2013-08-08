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
  double rhat[ngrid[0]], hopt=bw[0],alpha=0.0;

  for(i=0; i<ngrid[0]; i++){
    rhat[i] = 0.0;
  }

  if(lscv[0] == 1){
    if(adapt[0]==0){
      hopt = lprlscv(x,y,size[0],hopt);
      lpreg(xgrid, ngrid[0],x,y,size[0],hopt,rhat,ellx);
    }else{
      alprlscv(x,y,size[0],hopt, alpha);
      alpreg(xgrid, ngrid[0],x,y,size[0],hopt,alpha,rhat,ellx);
    }
  }else{
    lpreg(xgrid, ngrid[0],x,y,size[0],hopt,rhat,ellx);
  }
  
  bw[0] = hopt;
  for(i=0; i<ngrid[0]; i++){
    xgrid[i] = rhat[i];
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


/*  MLE for Weibull distribution */


double NR2p(double x0, 
	       double (*g)(double,double,double), 
	       double (*pg)(double,double),
	       double kappa, double alpha)
{
  int i=0,imax=1000;
  double dx,tol=0.00000001;
  dx=1.0;
  while(((dx>tol)|(dx/fabs(x0)>tol))&(i<imax)){
    dx =  (*g)(x0,kappa,alpha)/(*pg)(x0,kappa);
    x0 -= dx; dx = fabs(dx);
    i++;
  };
  if(i>=imax) x0=-999.;
  return x0;
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


/*  
 * product limit estimate for data with right censoring
 * Call: rcple(x,y,n[0],...);
 *
 */

void rcple(double x[], double w[], int n, double y[], double h[], int m) 
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



void weibullmae(double *x,double *w,int *size,
		double *pars, double *y,int *ny)
{
  int i, j, n=size[0], m=ny[0];
  double kappa=pars[0], lambda=pars[1];
  double HX[n], hx[n], x0, t3=0.0;
  //  rcple(x,w,n,y,Hx,m);
  rcple(x,w,n,x,HX,n);
  double t1,t2;
  t1 = 0.7644174 * pow(n,-.2);

  for(i=0; i<n;i++){
    // t2 = [f0''(x)]^2/f0(x)
    //    t3 = (kappa*kappa - 3.0*kappa +2.0)/(y[i]*y[i]) +
    //      kappa*kappa/pow(lambda, 2.0*kappa+1.0)*pow(y[i], 2.0*kappa-1.0) -
    //      3.0*kappa*(kappa-1.0)/pow(lambda, kappa)*pow(y[i], kappa-2.0);
    t3 = pow((kappa-1.0)/x[i]-kappa/lambda*pow(x[i]/lambda,kappa-1.),2.0) -
      (kappa-1.0)/(x[i]*x[i]) -
      kappa/(lambda*lambda)*(kappa-1.)*pow(x[i]/lambda,kappa-2.);
    t2 = dweibull(x[i],kappa,lambda,0)*t3*t3;
    hx[i] = t1 * pow(t2*HX[i],-.2);
  }

  for(i=0; i<m;i++){
    x0 = y[i]; y[i] = 0.0; // reuse y[]
    for(j=0; j<n; j++){
      t1 = (x0 - x[j])/hx[j];
      y[i] += w[j]/(HX[j]*hx[j])*dnorm(t1,0.,1.0,0);
    }
  }
  
  for(i=0; i<m;i++){
    y[i] /= n;
  }

}


void expmae(double *x,double *w,int *size,double *y,int *ny)
{
  int i, j, n=size[0], m=ny[0];
  double lambda=0.0,delta=0.0;
  double HX[n], hx[n], x0;
  for(i=0; i<n;i++){
    lambda += x[i];
    delta  += w[i];
  }
  lambda /= delta; //mle of lambda
  //  rcple(x,w,n,y,Hx,m);
  rcple(x,w,n,x,HX,n);
  double t1,t2;
  t1 = 0.7644174 * lambda * pow(n,-.2);
  t2 = 0.2/lambda;

  for(i=0; i<n;i++){
    hx[i] = t1 * exp(x[i]*t2) * pow(HX[i],-.2);
  }

  for(i=0; i<m;i++){
    x0 = y[i]; y[i] = 0.0; // reuse y[]
    for(j=0; j<n; j++){
      t1 = (x0 - x[j])/hx[j];
      y[i] += w[j]/(HX[j]*hx[j])*dnorm(t1,0.,1.0,0);
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


/*  Computer the MLE for Weibull distribution (Raw data) */
static double WeibullLlk(int npar, double *pars, void *ex)
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
    res += w[i]*(log(kappa) + (kappa-1.0)*log(x[i]) - kappa*log(lambda)
		 -  pow(x[i]/lambda, kappa));
  }
  
  return(-res);
}

void WeibullMle(double *x,double *w, int *nx, double *pars)
{
  int i,n=nx[0],npar=2; //2-parameter distribution only
  double dpar[npar],opar[npar]; 
  dpar[0] = pars[0]; dpar[1] = pars[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[2*n+1];
  yaux[0] = n; //sample size
  for(i=0;i<n;i++){
    yaux[i+1] = x[i];
    yaux[i+n+1] = w[i];
  }

  nmmin(npar,dpar,opar,&val,WeibullLlk,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  pars[0] = opar[0]; pars[1] = opar[1];
}

void WeibullMleNMMIN(double *x,double *w, int *nx, double *pars)
{
  int i,n=nx[0],npar=2; //2-parameter distribution only
  double dpar[npar],opar[npar]; 
  dpar[0] = pars[0]; dpar[1] = pars[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[2*n+1];
  yaux[0] = n; //sample size
  for(i=0;i<n;i++){
    yaux[i+1] = x[i];
    yaux[i+n+1] = w[i];
  }

  nmmin(npar,dpar,opar,&val,WeibullLlk,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  pars[0] = opar[0]; pars[1] = opar[1];
}

/*  Computer the MLE for Normal distribution (Censored data) */
static double CNormalLlk(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  double *tmp= (double*)ex, res=0.0;
  int i,n = (int)tmp[0]; //first element is the length of x;
  double kappa = pars[0], lambda= pars[1]; 
  double w[n], a[n], b[n];
  
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    a[i] = tmp[i+1]; 
    b[i] = tmp[i+n+1]; 
    w[i] = tmp[i+2*n+1]; 
  }
  
  for(i=0;i<n;i++) {
    if(b[i] > 0.0) {
      if(b[i] > a[i]){
	res += w[i]*(pnorm(b[i],kappa,lambda,1,0)
		     - pnorm(a[i],kappa,lambda,1,0));
      }else{
	res += w[i]*dnorm(a[i],kappa,lambda,0);
      }
    }else{
      res += w[i]*(1.0 - pnorm(a[i],kappa,lambda,1,0));
    }
  }
  
  return(-res);
}

void CNormalMle(double *w, double *a, double *b,
		int *nx, double *pars)
{
  int i,n=nx[0],npar=2; //2-parameter distribution only
  double dpar[npar],opar[npar]; 
  dpar[0] = pars[0]; dpar[1] = pars[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[3*n+1];
  yaux[0] = n; //sample size
  for(i=0;i<n;i++){
    yaux[i+1] = a[i];
    yaux[i+n+1] = b[i];
    yaux[i+2*n+1] = w[i];
  }

  nmmin(npar,dpar,opar,&val,CNormalLlk,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  pars[0] = opar[0]; pars[1] = opar[1];
}


/*  Computer the MLE for Weibull distribution (Censored data) */
static double CWeibullLlk(int npar, double *pars, void *ex)
// to be called by the Nelder-Mead simplex method
{
  double *tmp= (double*)ex, res=0.0;
  int i,n = (int)tmp[0]; //first element is the length of x;
  double kappa = pars[0], lambda= pars[1]; 
  double w[n], a[n], b[n];
  
  for(i=0;i<n;i++) {//restore auxiliary information from ex;
    a[i] = tmp[i+1]; 
    b[i] = tmp[i+n+1]; 
    w[i] = tmp[i+2*n+1]; 
  }
  
  for(i=0;i<n;i++) {
    if(b[i] > 0.0) {
      if(b[i] > a[i]){
	res += w[i]*(pweibull(b[i],kappa,lambda,1,0)
		     - pweibull(a[i],kappa,lambda,1,0));
      }else{
	res += w[i]*dweibull(a[i],kappa,lambda,0);
      }
    }else{
      res += w[i]*(1.0 - pweibull(a[i],kappa,lambda,1,0));
    }
  }
  
  return(-res);
}

void CWeibullMle(double *w, double *a, double *b,
		int *nx, double *pars)
{
  int i,n=nx[0],npar=2; //2-parameter distribution only
  double dpar[npar],opar[npar]; 
  dpar[0] = pars[0]; dpar[1] = pars[1]; //initial values
  double abstol=0.00000000001,reltol=0.0000000000001,val;
  int ifail=0,trace=0, maxit=1000, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double yaux[3*n+1];
  yaux[0] = n; //sample size
  for(i=0;i<n;i++){
    yaux[i+1] = a[i];
    yaux[i+n+1] = b[i];
    yaux[i+2*n+1] = w[i];
  }

  nmmin(npar,dpar,opar,&val,CWeibullLlk,&ifail,abstol,reltol, 
	(void *)yaux,alpha,beta,gamma,trace,&fncount,maxit);
  pars[0] = opar[0]; pars[1] = opar[1];
}

double kbiweight(double x)
{
  double t2, res=0.0;
  if(fabs(x) <= 1.0){
    t2 = 1.0 - x * x;
    res = 0.9375 * t2 * t2;
  }
  return res;
}
double akbiweight(double x, double q)
{
  double t2,t3,a,b,res=0.0;
  if(x >=-1 && x<=q){
    t3 = pow(1.0+q,5.0) *
      (81.-168*q+126*q*q-40.*pow(q,3.0) + 5.0*pow(q,4.0));
    a = 64.*(8.0-24.*q+48.*q*q-45*q*q*q+15*pow(q,4.))/t3;
    b = 1120.*pow(1.-q,3.0)/t3;
    t2 = 1.0 - x * x;
    res = 0.9375 * t2 * t2*(a+b*x);
  }
  return res;
}

double kunif(double x)
{
  double res=0.0;
  if(fabs(x) <= 1.0){
    res = 0.5;
  }
  return res;
}
double akunif(double x, double q)
{
  double res=0.0;
  if(x >=-1 && x<=q){
    res = 4.*(1+pow(q,3.0))/pow(1+q,-4.) +
      6.*(1-q)*pow(1+q, -3.0) * x;
  }
  return res;
}

double kepan(double x)
{
  double res=0.0;
  if(fabs(x) <= 1.0){
    res = 0.75 * (1.-x * x);
  }
  return res;
}
double akepan(double x, double q)
{
  double t3,a,b,res=0.0;
  if(x >=-1 && x<=q){
    t3 = pow(1.+q, 4.0) * (19.0-18.0*q+3.0*q*q);
    a = 64.*(2.-4.*q+6.*q*q-3.*pow(q,3.))/t3;
    b = 240.*pow(1.-q,2.0)/t3;
    res = 0.75 * (1.-x*x)*(a+b*x);
  }
  return res;
}

void smhazard(double xgrid[], int m, double x[], double y[],
	      int n, double h, int ikernel, double ht[]) 
{
  int i,j;
  double q, z;
  switch ( ikernel ) {
  case 1: //biweight
    for(i=0; i<m; i++){
      ht[i] = 0.0;
      for(j=0; j<n; j++){
	z = (xgrid[i]-x[j])/h;
	if(xgrid[i] >= h && xgrid[i] <= x[n-1]-h){
	  ht[i] +=  kbiweight(z) * y[j];
	}else if(xgrid[i] < h){
	  q = xgrid[i]/h;
	  ht[i] +=  akbiweight(z,q) * y[j];
	}else{
	  q = (x[n-1]-xgrid[i])/h;
	  ht[i] +=  akbiweight(-z,q) * y[j];
	}
      }
      ht[i] /= h;
    }
    break;
  case 2: // uniform
    for(i=0; i<m; i++){
      ht[i] = 0.0;
      for(j=0; j<n; j++){
	z = (xgrid[i]-x[j])/h;
	if(xgrid[i] >= h && xgrid[i] <= x[n-1]-h){
	  ht[i] +=  kunif(z) * y[j];
	}else if(xgrid[i] < h){
	  q = xgrid[i]/h;
	  ht[i] +=  akunif(z,q) * y[j];
	}else{
	  q = (x[n-1]-xgrid[i])/h;
	  ht[i] +=  akunif(-z,q) * y[j];
	}
      }
      ht[i] /= h;
    }
    break;
  default: //epanechinikov
    for(i=0; i<m; i++){
      ht[i] = 0.0;
      for(j=0; j<n; j++){
	z = (xgrid[i]-x[j])/h;
	if(xgrid[i] >= h && xgrid[i] <= x[n-1]-h){
	  ht[i] +=  kepan(z) * y[j];
	}else if(xgrid[i] < h){
	  q = xgrid[i]/h;
	  ht[i] +=  akepan(z,q) * y[j];
	}else{
	  q = (x[n-1]-xgrid[i])/h;
	  ht[i] +=  akepan(-z,q) * y[j];
	}
      }
      ht[i] /= h;
    }
    break;
  }
}

void smvhazard(double xgrid[], int m, double x[], double y[],
	      int n, double h, int ikernel, double ht[]) 
{
  int i,j;
  double q, z, t;
  switch ( ikernel ) {
  case 1: //biweight
    for(i=0; i<m; i++){
      ht[i] = 0.0;
      for(j=0; j<n; j++){
	z = (xgrid[i]-x[j])/h;
	if(xgrid[i] >= h && xgrid[i] <= x[n-1]-h){
	  t = kbiweight(z);
	  ht[i] +=  t * t * y[j];
	}else if(xgrid[i] < h){
	  q = xgrid[i]/h;
	  t = akbiweight(z,q);
	  ht[i] +=  t * t * y[j];
	}else{
	  q = (x[n-1]-xgrid[i])/h;
	  t = akbiweight(-z,q);
	  ht[i] +=  t * t * y[j];
	}
      }
      ht[i] /= h*h;
    }
    break;
  case 2: // uniform
    for(i=0; i<m; i++){
      ht[i] = 0.0;
      for(j=0; j<n; j++){
	z = (xgrid[i]-x[j])/h;
	if(xgrid[i] >= h && xgrid[i] <= x[n-1]-h){
	  t = kunif(z);
	  ht[i] +=  t * t * y[j];
	}else if(xgrid[i] < h){
	  q = xgrid[i]/h;
	  t = akunif(z,q);
	  ht[i] +=  t * t * y[j];
	}else{
	  q = (x[n-1]-xgrid[i])/h;
	  t = akunif(-z,q);
	  ht[i] +=  t * t * y[j];
	}
      }
      ht[i] /= h*h;
    }
    break;
  default: //epanechinikov
    for(i=0; i<m; i++){
      ht[i] = 0.0;
      for(j=0; j<n; j++){
	z = (xgrid[i]-x[j])/h;
	if(xgrid[i] >= h && xgrid[i] <= x[n-1]-h){
	  t = kepan(z);
	  ht[i] +=  t * t * y[j];
	}else if(xgrid[i] < h){
	  q = xgrid[i]/h;
	  t = akepan(z,q);
	  ht[i] +=  t * t * y[j];
	}else{
	  q = (x[n-1]-xgrid[i])/h;
	  t = akepan(-z,q);
	  ht[i] +=  t * t * y[j];
	}
      }
      ht[i] /= h*h;
    }
    break;
  }
}

double hazardlscv(double x[], double y[],
		int n, double h, int ikernel) 
{
  int i,j;
  double q, z, t, gh=0.0;
  switch ( ikernel ) {
  case 1: //biweight
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
	if(i != j){
	  z = (x[i]-x[j])/h;
	  if(x[i] >= h && x[i] <= x[n-1] - h){
	    t = kbiweight(z);
	    gh += t * t * y[j] * y[i];
	  }else if(x[i] < h){
	    q = x[i]/h;
	    t = akbiweight(z,q);
	    gh +=  t * t * y[j] * y[i];
	  }else{
	    q = (x[n-1]-x[i])/h;
	    t = akbiweight(-z,q);
	    gh +=  t * t * y[j] * y[i];
	  }
	}
      }
    }
    break;
  case 2: // uniform
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
	if(i != j){
	  z = (x[i]-x[j])/h;
	  if(x[i] >= h && x[i] <= x[n-1]-h){
	    t = kunif(z);
	    gh +=  t * t * y[j] * y[i];
	  }else if(x[i] < h){
	    q = x[i]/h;
	    t = akunif(z,q);
	    gh +=  t * t * y[j] * y[i];
	  }else{
	    q = (x[n-1]-x[i])/h;
	    t = akunif(-z,q);
	    gh +=  t * t * y[j] * y[i];
	  }
	}
      }
    }
    break;
  default: //epanechinikov
    for(i=0; i<n; i++){
      for(j=0; j<n; j++){
	if(i != j){
	  z = (x[i]-x[j])/h;
	  if(x[i] >= h && x[i] <= x[n-1]-h){
	    t = kepan(z);
	    gh +=  t * t * y[j] * y[i];
	  }else if(x[i] < h){
	    q = x[i]/h;
	    t = akepan(z,q);
	    gh +=  t * t * y[j] * y[i];
	  }else{
	    q = (x[n-1]-x[i])/h;
	    t = akepan(-z,q);
	    gh +=  t * t * y[j] * y[i];
	  }
	}
      }
    }
    break;
  }
  return gh;
}

double hsmhazard(double xgrid[], int m, double x[], double y[],
	       int n, double h, int ikernel) 
{
  int i,j,k;
  double ht[m];
  double h0, dh = 0.05 * h, hopt, gh, gh2, gmin=1.0e+10;
  h0 = dh;
  for(i=0; i<100; i++){
    smhazard(xgrid, m,x,y,n,h0,ikernel,ht);
    gh = 0.0;
    for(j=1;j<m; j++){
      gh += 0.5 * (xgrid[j] - xgrid[j-1]) * 
	(ht[j]*ht[j]+ht[j-1]*ht[j-1]);
    }
    gh2 = 0.0;
    gh2 = hazardlscv(x,y,n,h0,ikernel);
    gh -= 2.0 * gh2/h0;
    if(gh<gmin){
      gmin = gh;
      hopt = h0;
    }
    h0 += dh;
  }
  return hopt;
}

void lpshazard(double *xgrid, double *xvar, int *ngrid, 
	       double *x, double *y, double *v, int *size, 
	       double *bw, int *lscv, int *ikernel)
{
  int i,m=ngrid[0],n=size[0];
  double hopt=bw[0], ht[m];

  // initial V(ht) and ht
  for(i=0; i<m; i++){
    xvar[i] = 0.0;
    ht[i] = 0.0;
  }

  if(lscv[0] == 1){// find hMISE
    hopt = hsmhazard(xgrid, m,x,y,n,hopt,ikernel[0]);
  }
  smhazard(xgrid, m,x,y,n,hopt,ikernel[0],ht);
  // compute V(ht)
  smvhazard(xgrid, m,x,v,n,hopt,ikernel[0],xvar);
  // return results to R
  bw[0] = hopt;
  for(i=0; i<m; i++){
    xgrid[i] = ht[i];
  }
}
