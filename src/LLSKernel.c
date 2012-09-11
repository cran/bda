/*
This function is migrated from an old ADA package (2012/02/12)

We add c subroutines for lls using gaussian kernel and the case with
laplace errors using gaussian kernels (2012/03/30).

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

double GLInt1p(double a, double b, double (*f)(double))
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

double GL3pInt(double a, double b, double (*f)(double,double,double),
	       double x,double sig)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*f)(d - dum,x,sig) + (*f)(d + dum,x,sig) );
   }

   return c * integral;
}

double GLIntvvvi(double a, double b, 
		 double (*f)(double, double*,double*,double*,int),
		 double *x,double *sig,double *h,int n)
{
   double integral = 0.0; 
   double c = 0.5 * (b - a);
   double d = 0.5 * (b + a);
   double dum;
   const double *pB = &B100[NOPZ100 - 1];
   const double *pA = &A100[NOPZ100 - 1];

   for (; pB >= B100; pA--, pB--) {
      dum = c * *pB;
      integral += *pA * ( (*f)(d - dum,x,sig,h,n) + (*f)(d + dum,x,sig,h,n) );
   }

   return c * integral;
}

double funSupport(double t,double x,double sig){
  double t2=t*t;
  return cos(t*x)*pow(1.0-t2,3.0)*exp(0.5*sig*sig*t2);
}

double NRsolve(double x0, 
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

double ftube(double x, double kappa,double alpha){
  return 2.-2.*pnorm(x,0.0,1.0,1,0)+kappa/M_PI*exp(-.5*x*x)-alpha;
}

double fdertube(double x, double kappa){
  return -2.*dnorm(x,0.0,1.0,0)-x*kappa/M_PI*exp(-.5*x*x);
}

void pdfHSupport(double *x0, int *nx0, double *x, double *sig, double *bw, int *size)
{
  Fun3p f[1];
  f[0] = funSupport;
  int i,j,m,n;
  n = size[0];
  m = nx0[0];
  double fx[m],s_h[n],z,h;

  h = bw[0];

  for(i=0;i<m;i++){
    fx[i] = 0.0;
  }
  for(i=0;i<n;i++){
    s_h[i] = sig[i]/h;
  }
  
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      z = (x0[i]-x[j])/h;
      fx[i] += GL3pInt(0.0,1.0,f[0],z,s_h[j]);
    }
  }
  for(i=0;i<m;i++){
    x0[i] = fx[i] /( n*h*M_PI );
  }
}

//  Laplace kernel

double lkernel(double x,double sigh){
  double out;
  out = dnorm(x,0.0,1.0,0)*(1.0+pow(sigh,2.0)*(1.0-x*x));
  if(out<0) out=0.0;
  return out;
}

double delkernel(double x,double sigh){
  double out,t0,sigh2;
  sigh2 = sigh*sigh;
  out = 1.0 + sigh2*(1.0-x*x);
  if(out<0){
    out=0.0000001;
  }else{
    t0 = dnorm(x,0.0,1.0,0);
    out = -x*t0*out-2.0*x*t0*sigh2;
  }
  return out;
}

void pdfHLaplace(double *x0, int *nx0, double *x, double *sig, double *bw, int *size)
{
  int i,j,m,n;
  m = nx0[0];
  n = size[0];
  double fx[m],s_h[n],z,h;

  h = bw[0];

  for(i=0;i<m;i++){
    fx[i] = 0.0;
  }
  for(i=0;i<n;i++){
    s_h[i] = sig[i]/h;
  }
  
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      z = (x0[i]-x[j])/h;
      fx[i] += lkernel(z,s_h[j]);
    }
  }
  for(i=0;i<m;i++){
    x0[i] = fx[i] /( n*h );
  }
}

double funKappa(double t,double *x,double *sig, double *h, int n){
  int i;
  double Sn1,Sn2,deSn1,deSn2,dx,z,Sbi,Sdbi;
  double bi[n],dbi[n];
  double t0;
  double SLl,Sll,SLL;
  Sn1 = 0.0; Sn2 = 0.0; deSn1 = 0.0; deSn2 = 0.0;

  for(i=0;i<n;i++){
    dx = t-x[i];
    z = dx/h[i];
    t0 = sig[i]/h[i];
    bi[i] = lkernel(z,t0); //store K to compute bi later on
    dbi[i] = delkernel(z,t0)/h[i];//store K' to compute bi later on
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

void llrHLaplace(double *range, double *x, double *sig, double *y,double *scb,
		 double *bw, int *bwyn, int *size,double *sensitivity,double *level)
{
  int i,j,k,l,m,n;
  n = size[0];m=n;
  double rn[m],rnhat[n],s_h[n],bi[n],z,dx,kz,h;
  double Sn1,Sn2,Sbi;

  double fx[n],g,alp,dalp,dh,lcvscore,xmax=999999999.;
  double alp0,h0=bw[0];
  h = bw[0];

  //  Fun2d gx[1];
  //  gx[0] = fdertube;
  //  Fun3d dgx[1];
  //  dgx[0] = ftube;

  Fvvvi f[1];
  f[0] = funKappa;

  /*lcv score
  search alp from 0 to 0.1 iter=100: increment 0.001
  search h from 0.5bw[0] to 3bw[0], iter=25: increment 0.1*bw[0]
   */
  if(bwyn[0]==1){
  dalp = 0.002;
  dh = 0.05*h;
  h = 0.5*h;h0=h;
  alp=0.0;alp0=alp;

  for(i=0;i<n;i++){//initialize scbs
    scb[i] =0.0;
  }

  for(i=0;i<50;i++){
    h += dh;
    for(k=0;k<n;k++){
      fx[k] = 0.0;
      s_h[k] = sig[k]/h;
    }
    
    Sn1 = 0.0;  //reuse
    for(k=0;k<n;k++){
      for(l=0;l<n;l++){
	z = (x[k]-x[l])/h;
	fx[k] += lkernel(z,s_h[l]);
      }
      if(fx[k]>0) Sn1 += log(fx[k]);
    }
    g = exp(Sn1/n);
    alp=0.0;
    for(j=0;j<50;j++){

      for(k=0;k<n;k++){
	rnhat[k] = 0.0;
      }
      for(k=0;k<n;k++){
	s_h[k] = h*pow(fx[i]/g,-alp);
      }
  
      lcvscore =0.0;
      for(k=0;k<n;k++){
	Sn1 = 0.0;
	Sn2 = 0.0;
	for(l=0;l<n;l++){
	  dx = x[k]-x[l];
	  z = dx/s_h[l];
	  kz = lkernel(z,sig[l]/s_h[l]);
	  bi[l] = kz;
	  Sn1 += kz * dx;
	  Sn2 += kz * dx * dx;
	}
	Sbi = 0.0;
	for(l=0;l<n;l++){
	  dx = x[k]-x[l];
	  bi[l] *= Sn2-Sn1*dx;
	  Sbi += bi[l];
	}
	for(l=0;l<n;l++){
	  rnhat[k] += bi[l] * y[l]/Sbi;
	}
	lcvscore += pow((y[k] - rnhat[k])/(1.- bi[k]/Sbi),2.0);
      }
      if(lcvscore<xmax){
	xmax = lcvscore;
	alp0 = alp; h0 = h;
      }
      alp += dalp;
    }



  }
  }else{
    alp0 = 0.0;
  }

  /*  output from this part: alp and h
   */
  h = h0;bw[0] = h;
  alp = alp0;sensitivity[0]=alp;

  for(k=0;k<n;k++){
    fx[k] = 0.0;
    s_h[k] = sig[k]/h;
  }
    
  Sn1 = 0.0;  //reuse
  for(k=0;k<n;k++){
    for(l=0;l<n;l++){
      z = (x[k]-x[l])/h;
      fx[k] += lkernel(z,s_h[l]);
    }
    if(fx[k]>0) Sn1 += log(fx[k]);
  }
  g = exp(Sn1/n);

  for(i=0;i<m;i++){
    rn[i] = 0.0;
  }

  for(i=0;i<n;i++){
    s_h[i] = h*pow(fx[i]/g,-alp);
  }
  
  for(i=0;i<m;i++){
    Sn1 = 0.0;
    Sn2 = 0.0;
    for(j=0;j<n;j++){
      dx = x[i]-x[j];
      z = dx/s_h[j];
      kz = lkernel(z,sig[j]/s_h[j]);
      bi[j] = kz;
      Sn1 += kz * dx;
      Sn2 += kz * dx * dx;
    }
    Sbi = 0.0;
    for(j=0;j<n;j++){
      dx = x[i]-x[j];
      bi[j] *= Sn2-Sn1*dx;
      Sbi += bi[j];
      scb[i] += bi[j]*bi[j];
    }
    if(Sbi==0.0) Sbi=0.000001;
    scb[i] = sqrt(scb[i])/Sbi; // this is really ||\ell(x)||
    for(j=0;j<n;j++){
      rn[i] += bi[j] * y[j]/Sbi;
    }
  }

  //compute kappa0
  range[0] = GLIntvvvi(range[0],range[1],f[0],x,sig,s_h,n); 
  //compute c
  i=0;z=2.0;//initial value
  dx=10.0;
  while((i<100)&(fabs(dx)>0.000001)){
    dx = (2.*(1-pnorm(z,0.,1.,1,0))+range[0]/M_PI*exp(-.5*z*z)-1+level[0])/
      (2.*dnorm(z,0.,1.,0)+z*range[0]/M_PI*exp(-0.5*z*z));
    z += dx;
    i++;
  }
  if(i>=100) z=-999.;
  range[1] = z;
  //out results for rn(x) and partial SCB = c * ||l(x)||
  for(i=0;i<m;i++){
    s_h[i] = log((y[i]-rn[i])*(y[i]-rn[i]));  // reuse s_h as z = log(Yi-rn(xi))^2
    if(!ISNAN(rn[i])) 
      y[i] = rn[i];
    scb[i] *= z;
  }
  //Compute s(x)

  dh = 0.05*h;
  h = 0.5*h;h0=h;
  for(i=0;i<50;i++){
    h += dh;
    for(k=0;k<n;k++){
      fx[k] = 0.0;
    }
    
    for(k=0;k<n;k++){
      rnhat[k] = 0.0;
    }
    
    lcvscore =0.0;
    for(k=0;k<n;k++){
      Sn1 = 0.0;
      Sn2 = 0.0;
      for(l=0;l<n;l++){
	dx = x[k]-x[l];
	z = dx/h;
	kz = dnorm(z,0.0,1.0,0);
	bi[l] = kz;
	Sn1 += kz * dx;
	Sn2 += kz * dx * dx;
      }
      Sbi = 0.0;
      for(l=0;l<n;l++){
	dx = x[k]-x[l];
	bi[l] *= Sn2-Sn1*dx;
	Sbi += bi[l];
      }
      for(l=0;l<n;l++){
	rnhat[k] += bi[l] * s_h[l]/Sbi;
      }
      lcvscore += pow((s_h[k] - rnhat[k])/(1.- bi[k]/Sbi),2.0);
    }
    if(lcvscore<xmax){
      xmax = lcvscore;
      h0 = h;
    }
  }


  /*  output from this part: alp and h

   */
  h = h0;bw[0] = h;
      
  for(i=0;i<m;i++){
    rn[i] = 0.0;
  }
  
  for(i=0;i<m;i++){
    Sn1 = 0.0;
    Sn2 = 0.0;
    for(j=0;j<n;j++){
      dx = x[i]-x[j];
      z = dx/h;
      kz = dnorm(z,0.0,1.0,0);
      bi[j] = kz;
      Sn1 += kz * dx;
      Sn2 += kz * dx * dx;
    }
    Sbi = 0.0;
    for(j=0;j<n;j++){
      dx = x[i]-x[j];
      bi[j] *= Sn2-Sn1*dx;
      Sbi += bi[j];
    }
    for(j=0;j<n;j++){
      rn[i] += bi[j] * s_h[j]/Sbi;
    }
    
  }
  for(i=0;i<m;i++){
    scb[i] *= sqrt(exp(rn[i]));
  }


}


void llrGauss(double *x, double *y, int *nx, double *x0, int *ngrid, double *bw, int *lscv)
{
  int i,j,k,l,m,n;
  n = nx[0]; m = ngrid[0];
  double rn[m],rnhat[n],bi[n],z,dx,kz,h;
  double Sn1,Sn2,Sbi;

  double dh,lcvscore,xmax=999999999.;
  double h0;
  h = bw[0];

  /*lcv score
    search h from 0.5bw[0] to 3bw[0], iter=25: increment 0.1*bw[0]
  */
  if(lscv[0] == 1){
    dh = 0.05*h;
    h = 0.5*h;h0=h;
    for(i=0;i<50;i++){
      h += dh;
      
      for(k=0;k<n;k++){
	rnhat[k] = 0.0;
      }
      
      lcvscore =0.0;
      for(k=0;k<n;k++){
	Sn1 = 0.0;
	Sn2 = 0.0;
	for(l=0;l<n;l++){
	  dx = x[k]-x[l];
	  z = dx/h;
	  kz = dnorm(z,0.0,1.0,0);
	  bi[l] = kz;
	  Sn1 += kz * dx;
	  Sn2 += kz * dx * dx;
	}
	Sbi = 0.0;
	for(l=0;l<n;l++){
	  dx = x[k]-x[l];
	  bi[l] *= Sn2-Sn1*dx;
	  Sbi += bi[l];
	}
	for(l=0;l<n;l++){
	  rnhat[k] += bi[l] * y[l]/Sbi;
	}
	lcvscore += pow((y[k] - rnhat[k])/(1.- bi[k]/Sbi),2.0);
      }
      if(lcvscore<xmax){
	xmax = lcvscore;
	h0 = h;
      }
    }
    /*  output from this part: alp and h   */
    h = h0;bw[0] = h;
  }

  for(i=0;i<m;i++){
    Sn1 = 0.0;
    Sn2 = 0.0;
    for(j=0;j<n;j++){
      dx = x0[i]-x[j];
      z = dx/h;
      kz = dnorm(z,0.0,1.0,0);
      bi[j] = kz;
      Sn1 += kz * dx;
      Sn2 += kz * dx * dx;
    }
    Sbi = 0.0;
    for(j=0;j<n;j++){
      dx = x0[i]-x[j];
      bi[j] *= Sn2-Sn1*dx;
      Sbi += bi[j];
    }
    rn[i] = 0.0;
    for(j=0;j<n;j++){
      rn[i] += bi[j] * y[j];
    }
    rn[i] /= Sbi;
  }
  for(i=0;i<m;i++){
    x0[i] = rn[i];
  }
}

void lprLaplace(double *range, double *x, double *sig, double *y,double *scb,
		 double *bw, int *size,double *sensitivity,double *level)
{
  int i,j,k,l,m,n;
  n = size[0];m=n;
  double rn[m],rnhat[n],s_h[n],bi[n],z,dx,kz,h;
  double Sn1,Sn2,Sbi, sigx[n], sigy;

  double fx[n],g,alp,dalp,dh,lcvscore,xmax=999999999.;
  double alp0,h0;
  h = bw[0];

  //  Fun2d gx[1];
  //  gx[0] = fdertube;
  //  Fun3d dgx[1];
  //  dgx[0] = ftube;

  Fvvvi f[1];
  f[0] = funKappa;

  /*lcv score
  search alp from 0 to 0.1 iter=100: increment 0.001
  search h from 0.5bw[0] to 3bw[0], iter=25: increment 0.1*bw[0]
   */
  dalp = 0.002;
  dh = 0.05*h;
  h = 0.5*h;h0=h;
  alp=0.0;alp0=alp;

  for(i=0;i<n;i++){//initialize scbs
    scb[i] =0.0;
    sigx[i] = sig[0];
  }
  sigy = sig[1];

  for(i=0;i<50;i++){
    h += dh;
    for(k=0;k<n;k++){
      fx[k] = 0.0;
      s_h[k] = sigx[k]/h;
    }
    
    Sn1 = 0.0;  //reuse
    for(k=0;k<n;k++){
      for(l=0;l<n;l++){
	z = (x[k]-x[l])/h;
	fx[k] += lkernel(z,s_h[l]);
      }
      if(fx[k]>0) Sn1 += log(fx[k]);
    }
    g = exp(Sn1/n);
    alp=0.0;
    for(j=0;j<50;j++){

      for(k=0;k<n;k++){
	rnhat[k] = 0.0;
      }
      for(k=0;k<n;k++){
	s_h[k] = h*pow(fx[i]/g,-alp);
      }
  
      lcvscore =0.0;
      for(k=0;k<n;k++){
	Sn1 = 0.0;
	Sn2 = 0.0;
	for(l=0;l<n;l++){
	  dx = x[k]-x[l];
	  z = dx/s_h[l];
	  kz = lkernel(z,sigx[l]/s_h[l]);
	  bi[l] = kz;
	  Sn1 += kz * dx;
	  Sn2 += kz * dx * dx;
	}
	Sbi = 0.0;
	for(l=0;l<n;l++){
	  dx = x[k]-x[l];
	  bi[l] *= Sn2-Sn1*dx;
	  Sbi += bi[l];
	}
	for(l=0;l<n;l++){
	  rnhat[k] += bi[l] * y[l]/Sbi;
	}
	lcvscore += pow((y[k] - rnhat[k])/(1.- bi[k]/Sbi),2.0);
      }
      if(lcvscore<xmax){
	xmax = lcvscore;
	alp0 = alp; h0 = h;
      }
      alp += dalp;
    }



  }


  /*  output from this part: alp and h

   */
  h = h0;bw[0] = h;
  alp = alp0;sensitivity[0]=alp;

  for(k=0;k<n;k++){
    fx[k] = 0.0;
    s_h[k] = sigx[k]/h;
  }
    
  Sn1 = 0.0;  //reuse
  for(k=0;k<n;k++){
    for(l=0;l<n;l++){
      z = (x[k]-x[l])/h;
      fx[k] += lkernel(z,s_h[l]);
    }
    if(fx[k]>0) Sn1 += log(fx[k]);
  }
  g = exp(Sn1/n);

  for(i=0;i<m;i++){
    rn[i] = 0.0;
  }

  for(i=0;i<n;i++){
    s_h[i] = h*pow(fx[i]/g,-alp);
  }
  
  for(i=0;i<m;i++){
    Sn1 = 0.0;
    Sn2 = 0.0;
    for(j=0;j<n;j++){
      dx = x[i]-x[j];
      z = dx/s_h[j];
      kz = lkernel(z,sigx[j]/s_h[j]);
      bi[j] = kz;
      Sn1 += kz * dx;
      Sn2 += kz * dx * dx;
    }
    Sbi = 0.0;
    for(j=0;j<n;j++){
      dx = x[i]-x[j];
      bi[j] *= Sn2-Sn1*dx;
      Sbi += bi[j];
      scb[i] += bi[j]*bi[j];
    }
    if(Sbi==0.0) Sbi=0.000001;
    scb[i] = sqrt(scb[i])/Sbi; // this is really ||\ell(x)||
    for(j=0;j<n;j++){
      rn[i] += bi[j] * y[j]/Sbi;
    }
  }

  //compute kappa0
  range[0] = GLIntvvvi(range[0],range[1],f[0],x,sig,s_h,n); 
  //compute c
  i=0;z=2.0;//initial value
  dx=10.0;
  while((i<100)&(fabs(dx)>0.000001)){
    dx = (2.*(1-pnorm(z,0.,1.,1,0))+range[0]/M_PI*exp(-.5*z*z)-1+level[0])/
      (2.*dnorm(z,0.,1.,0)+z*range[0]/M_PI*exp(-0.5*z*z));
    z += dx;
    i++;
  }
  if(i>=100) z=-999.;
  range[1] = z;
  //out results for rn(x) and partial SCB = c * ||l(x)||
  for(i=0;i<m;i++){
    s_h[i] = log((y[i]-rn[i])*(y[i]-rn[i]));  // reuse s_h as z = log(Yi-rn(xi))^2
    //if(!ISNAN(rn[i])) 
    y[i] = rn[i];
    scb[i] *= z;
  }
  //Compute s(x)

  dh = 0.05*h;
  h = 0.5*h;h0=h;
  for(i=0;i<50;i++){
    h += dh;
    for(k=0;k<n;k++){
      fx[k] = 0.0;
    }
    
    for(k=0;k<n;k++){
      rnhat[k] = 0.0;
    }
    
    lcvscore =0.0;
    for(k=0;k<n;k++){
      Sn1 = 0.0;
      Sn2 = 0.0;
      for(l=0;l<n;l++){
	dx = x[k]-x[l];
	z = dx/h;
	kz = dnorm(z,0.0,1.0,0);
	bi[l] = kz;
	Sn1 += kz * dx;
	Sn2 += kz * dx * dx;
      }
      Sbi = 0.0;
      for(l=0;l<n;l++){
	dx = x[k]-x[l];
	bi[l] *= Sn2-Sn1*dx;
	Sbi += bi[l];
      }
      for(l=0;l<n;l++){
	rnhat[k] += bi[l] * s_h[l]/Sbi;
      }
      //      Lii[k] = bi[k]/Sbi;
      lcvscore += pow((s_h[k] - rnhat[k])/(1.- bi[k]/Sbi),2.0);
    }
    if(lcvscore<xmax){
      xmax = lcvscore;
      h0 = h;
    }
  }

  for(i=0;i<m;i++){
    scb[i] *= sigy;
  }
}

/*
The following codes deal with hetero laplacian errors with gaussian
kernel.  Data were binned to reduce the computational burden.
 */

double llsmooth(double *x, double *y, double x0, double *sig, 
		double h, int n){
  int i;
  double z,z1,z2,sh2,t1,S2, S1,tsum;
  double fx,kx[n];

  S1 = 0.0; 
  S2 = 0.0;
  for(i=0;i<n;i++){
    z1 = x[i] - x0;
    z = z1 / h;
    z2 = z*z;
    sh2 = pow(sig[i]/h,2.);
    kx[i] = (1.0 + sh2 * (1.0 - z2)) * exp(-0.5*z2)*M_1_SQRT_2PI;
    if(kx[i] < 0.0) kx[i] = 0.0;
    t1 = kx[i] * z1;
    S1 = S1 + t1;
    S2 = S2 + t1 * z1;
  }
  tsum = 0.0;
  fx = 0.0;
  for(i=0;i<n;i++){
    z1 = x[i] - x0;
    kx[i] = kx[i] * (S2 - z1 * S1); //li(x)
    tsum = tsum + kx[i]; //sum(li(x))
    fx += kx[i] * y[i];
  }
  if(tsum <= 0.0){
    fx = 0.00001;
  }else{
    fx /= tsum;
  }
  return fx;
}


void llslaplace(double *y, double *x, double *sig, int *size,
		double *bw, double *lbound,double *ubound, int *ngrid)
{
  int i,m,n;
  n = size[0];  m = ngrid[0];
  double h,a,b,x0[m],fx[m],delta,t0,tsum;
  a = lbound[0]; b=ubound[0];h=bw[0];
  delta = (b-a)/(m-1.0);
  tsum = a;
  for(i=0;i<m;i++){
    x0[i] = tsum;
    tsum += delta;
    fx[i] = 0.0;
  }
  for(i=0;i<m;i++){
    fx[i] = llsmooth(x,y,x0[i],sig,h,n);
  }

  for(i=0;i<n;i++){
    t0 = (x[i]-a)/delta;
    y[i] = fx[(int)ceil(t0)] * (t0 - floor(t0)) 
      + fx[(int)floor(t0)] * (ceil(t0)-t0);
  }

}

/*   2012/03/02: Without considering the measurement error.  Perform a
     regular local linear smoothing using Gaussian kernel.
 */


/*
  To compute the smoothing matrix ell_i(x)
*/
void CompSM(double x[], double f[], int n, double h, double ell[]){
  int i,j;
  double K[n][n],l[n][n];
  double S1[n], S2[n]; // two vectors S1(xj), S2(xj)
  double t1, bsum;
  
  // to compute the kernel values
  for(i=0;i<n;i++){
    for(j=i; j<n;j++){
      t1 = x[i] - x[j];
      K[i][j] = dnorm(t1/h,0.,1.,0);
      K[j][i] = K[i][j];
    }
  }
  // compute S1, S2.
  for(j=0; j<n;j++){
    S1[j] = 0.0; S2[j] = 0.0;
    for(i=0;i<n;i++){
      t1 = x[i] - x[j];
      S1[j] += f[i] * K[i][j] * t1;
      S2[j] += f[i] * K[i][j] * t1 * t1;
    }
  }
  // compute B and Lii: store sum(bi) to B, and bi to Lii
  for(j=0; j<n;j++){
    bsum = 0.0;
    for(i=0;i<n;i++){
      t1 = x[i] - x[j];
      l[i][j] = K[i][j] * (S2[j]- t1 * S1[j]);
      bsum += f[i] * l[i][j];
    }
    for(i=0;i<n;i++){
      ell[j*n+i] = l[i][j]/bsum;      
    }
  }
}


void CompASM(double x[], double f[], int n, double h, double lambda[], double ell[]){
  int i,j;
  double K[n][n],l[n][n];
  double S1[n], S2[n]; // two vectors S1(xj), S2(xj)
  double t1, bsum;
  
  // to compute the kernel values
  for(i=0;i<n;i++){
    for(j=i; j<n;j++){
      t1 = x[i] - x[j];
      K[i][j] = dnorm(t1/h/lambda[i],0.,1.,0);
      K[j][i] = dnorm(t1/h/lambda[j],0.,1.,0);;
    }
  }
  // compute S1, S2.
  for(j=0; j<n;j++){
    S1[j] = 0.0; S2[j] = 0.0;
    for(i=0;i<n;i++){
      t1 = x[i] - x[j];
      S1[j] += f[i] * K[i][j] * t1;
      S2[j] += f[i] * K[i][j] * t1 * t1;
    }
  }
  // compute B and Lii: store sum(bi) to B, and bi to Lii
  for(j=0; j<n;j++){
    bsum = 0.0;
    for(i=0;i<n;i++){
      t1 = x[i] - x[j];
      l[i][j] = K[i][j] * (S2[j]- t1 * S1[j]);
      bsum += f[i] * l[i][j];
    }
    for(i=0;i<n;i++){
      ell[j*n+i] = l[i][j]/bsum;      
    }
  }
}

/*
  To compute \sum(ell(x)Yi
*/
void rnhat(double y[], int n, double ell[], double rhat[]){
  /* Rh returns the GCV: = mean([(Yi-rhat)/(1-nu/n)]^2); we further
     approximate it by Cp */
  int i, j;
  for(i=0;i<n;i++){
    rhat[i] = 0;
    for(j=0;j<n;j++){
      rhat[i] += ell[i*n+j] * y[j];
    }
  }
}

double CompCp(double y[], double rhat[], double f[], int n, double ell[]){
  int i;
  double gcv, nsum=0.0;//, nu;
  gcv = 0.0; //nu=0.0;
  for(i=0;i<n;i++){
    nsum += f[i];
    //    nu += ell[i*n+i] * f[i];
    gcv += pow((y[i]/f[i] - rhat[i])/(1.0 - ell[i*n+i]), 2.0) * f[i];
  }
  //  return gcv/nsum*(1.0+2.0*nu/nsum);
  return gcv/nsum;
}


double CompVar(double y[], double rhat[], double f[], int n, double ell[]){
  int i,j;
  double gcv, nsum=0.0, nu1,nu2,tmp;
  gcv = 0.0; nu1=0.0; nu2=0.0;
  for(i=0;i<n;i++){
    nsum += f[i];
    nu1 += ell[i*n+i] * f[i];
    tmp = 0.0;
    for(j=0;j<n;j++){
      tmp += pow(ell[i*n+j],2.0) * f[j];
    }
    nu2 += tmp;
    gcv += pow((y[i]/f[i] - rhat[i])/(1.0 - ell[i*n+i]), 2.0) * f[i];
  }
  return gcv/(nsum-2.0*nu1+nu2);
}

void Compellx(double f[], int n, double ell[], double ellx[]){
  int i,j;
  double tmp;
  for(i=0;i<n;i++){
    tmp = 0.0;
    for(j=0;j<n;j++){
      tmp += pow(ell[j*n+i],2.0) * f[i];
    }
    ellx[i] =sqrt(tmp);
  }
}


/*
In this function, bw is given.  For large data sets, we bin the data
(x) to reduce the computational burden.

Input (1) y/x/f/size; xi a vector of distinct values, and fi is the
corresponding frequency.  yi is the sum of y in the bin of xi. x also
serves as the fine grid on which the smoothing function will be
evlauted. (2) bw is the bandwidth;

Output: rhat (using notation in AONS).
 */
void lpgauss(double *y, double *x, double *f, int *size,
	      double *bw, double *bws, int *nbw, double *sensitivity,
	      double *lx, double *rhat)
{

  int i, j, n = size[0];
  double h = bw[0], Cpmin=9.9e10, Cp,delta;
  double lambda[n], alp=0.0, alp0=0.0, g; //factor for adaptive bandwidth h*lambda
  double ell[n*n]; //this will not be huge (to save memory and time)

  /* choose optiml bandwidth using LCV score.  Can be approximated by
     GCV or a Cp.*/
  if(bws[0]>0){
    for(i=0;i<nbw[0];i++){
      CompSM(x, f, n, bws[i], ell);  // compute ell := [matrix]
      rnhat(y,n,ell,rhat); // compute rhat := [vector]
      Cp = CompCp(y,rhat,f,n,ell);
      if(Cp < Cpmin){
	h = bws[i]; Cpmin = Cp;
      }
      bws[i] = Cp;
    }
    bw[0] = h;
  }
  CompSM(x, f, n, h, ell);  // compute ell := [matrix]
  rnhat(y,n,ell,rhat); // compute rhat := [vector]
  if(sensitivity[0]>0.0){
    g = 0.0; 
    for(i=0;i<n;i++){
      if(rhat[i] > 0.0)	g += log(rhat[i]);
    }
    g = exp(g/n);
    delta = 1./(nbw[0]+1);
    for(i=0;i<nbw[0];i++){
      alp += delta;
      for(j=0;j<n;j++){
	lambda[j] = pow(rhat[j]/g, -alp);
      }
      CompASM(x, f, n, h, lambda, ell);  // compute ell := [matrix]
      rnhat(y,n,ell,rhat); // compute rhat := [vector]
      Cp = CompCp(y,rhat,f,n,ell);
      if(Cp < Cpmin){
	alp0 = alp; Cpmin = Cp;
      }
    }
    sensitivity[0] = alp0;

    for(i=0;i<n;i++){
      lambda[i] = pow(rhat[i]/g, -alp0);
    }
    CompASM(x, f, n, h, lambda, ell);  // compute ell := [matrix]
    rnhat(y,n,ell,rhat); // compute rhat := [vector]
  }
  Compellx(f,n,ell,lx);
}


/*   2012/02/10: 

     In this program, we apply local linear smoothing in the presence
     of laplace error.  Gaussian kernel is used by default.

 */
/*
  To compute a kernel for laplace errors with standard deviation 's'.  

  Input: x, h, and s

  Memo: h must be positive
 */
double MEKLaplace(double x, double h, double s){
  return dnorm(x,0.0,1.0,0)*(1.0+pow(s/h,2.0)*(1.0-x*x));
}

/* li(x) is a function to compute a vector 

   DMLaplace is to compute the design matrix of li(xj).  

   xv[] is a vector of length n.
 
   output: l[n,n]
*/

void SmoothMat(double xv[], int n, double h, double sig[], double dm[]){
  int i,j;
  double K[n][n],l[n][n];
  double S1[n], S2[n]; // two vectors S1(xj), S2(xj)
  double t1, bsum;
  
  // to compute the kernel values
  for(i=0;i<n;i++){
    for(j=i; j<n;j++){
      t1 = xv[i] - xv[j];
      K[i][j] = MEKLaplace(t1/h,h,sig[i]);
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

void DMLaplace(double *xv, int *size, double *bw, double *sig, double *dm){
  SmoothMat(xv, size[0], bw[0], sig, dm);
}

// for any given value x, compute a vector ell_i(x)
void li(double x, double xv[], int n, double h, double sig[], double l[]){
  int i;
  double K[n];
  double S1, S2; 
  double t1,bsum;
  
  // to compute the kernel values
  for(i=0;i<n;i++){
    t1 = xv[i] - x;
    K[i] = MEKLaplace(t1/h,h,sig[i]);
  }
  // compute S1, S2.
  for(i=0;i<n;i++){
    S1 = 0.0; S2 = 0.0;
    t1 = xv[i] - x;
    S1 += K[i] * t1;
    S2 += K[i] * t1 * t1;
  }
  // compute B and Lii: store sum(bi) to B, and bi to Lii
  bsum = 0.0;
  for(i=0;i<n;i++){
    t1 = xv[i] - x;
    l[i] = K[i] * (S2 - t1 * S1);
    bsum += l[i];
  }
  for(i=0;i<n;i++){
    l[i] = l[i]/bsum;      
  }
}

void LLSLaplace(double *y, double *x, double *sig, int *size,
		double *bw, double *rhat)
{
  int i,j,n=size[0];
  double dm[n*n];

  SmoothMat(x, n, bw[0], sig, dm);
  
  for(i=0;i<n;i++){
    rhat[i] = 0;
    for(j=0;j<n;j++){
      rhat[i] += dm[j*n+i] * y[j];
    }
  }
}
