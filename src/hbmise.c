#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

typedef double (*Fun7p)(double,double,double,double*,double*,double*,int);

//////////////////////////////////////////////////////////////////////////    
// double Gauss_Legendre_Integration_2pts( double a, double b, double (*f)(double) ) 
// void   Gauss_Legendre_Zeros_2pts( double nodes[] )                  
//    void   Gauss_Legendre_Coefs_2pts( double wght[] )                   
//////////////////////////////////////////////////////////////////////////

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
  double xbar,s,tmp1,tmp2,mise, fint,fb;
  double w2[n],w1_2[n],w1_2sq[n];
  double mise0=9999999999.;
  Fun7p fn[1];
  fn[0] = fa; //To compute part A to be integrated.

  tmp1 = 0.0;
  tmp2 = 0.0;
  for(i=0;i<n;i++){
    nsum += f[i];
    tmp1 = tmp1 + x[i]*f[i];
    tmp2 = tmp2 + x[i] * x[i]*f[i];
    w2[i] = pow(w[i],2.0);
    w1_2[i] = w[i]/2.0;
    w1_2sq[i]=pow(w1_2[i],2.0);
  }
  xbar = tmp1 / nsum;
  s = sqrt((tmp2 - nsum * xbar * xbar) / (nsum -1.));
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
