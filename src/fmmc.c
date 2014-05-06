#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#define INLINE

typedef double (*Funx)(double, double, double, double, int);

void kspvalue(double *pv){
  double ksp=0.,t;
  int i;
  t = pv[0];
  for(i=1;i<100;i++){
    ksp += exp(-2.*pow(i*t,2.));
    i++;
    ksp -= exp(-2.*pow(i*t,2.));
  }
  pv[0] = 2.*ksp;
}

double g0(double a, double c){
  //first derivative of logarithm Gamma(a)
  return(-digamma(a)+log(a)+c);
}

double g1(double a){
  return(-trigamma(a)+1./a);
}

void FitGamma(double *x, int *n, double *para){
  double ex1=0.,ex2=0., elx1=0., lex1=0.;
  double c=0.;
  int i;
  for(i=0;i<n[0];i++){
    ex1 += x[i];
    elx1 += log(x[i]);
  }
  ex1 = ex1/n[0];
  elx1 = elx1/n[0];
  lex1 = log(ex1);
  for(i=0;i<n[0];i++){
    ex2 += pow(x[i]-ex1,2.0);
  }
  ex2 = ex2/(n[0]-1.);
  c = elx1 - lex1;
  
  double x0=0.,x1, eps=0.00000001;
  int Iter=20;
  i=0;
  x1 = pow(ex1,2.0)/ex2; //mme of alpha
  //    x1 = ex2/ex1; //mme of beta
  while(i<Iter && fabs(x1-x0)>eps){
    x0 = x1;
    x1 = x1 - g0(x1,c)/g1(x1);
    i=i+1;
  }
  x0 = x1/ex1;
    
  para[0] = x1;
  para[1] = x0;
}

void FitBeta(double *x, int *n, double *para){
  double ex1=0.,ex2=0., elx1=0., elx2=0.;
  int i;
  for(i=0;i<n[0];i++){
    ex1 += x[i];
    elx1 += log(x[i]);
    elx2 += log(1.-x[i]);
  }
  ex1 = ex1/n[0];
  elx1 = elx1/n[0];
  elx2 = elx2/n[0];
  for(i=0;i<n[0];i++){
    ex2 += pow(x[i]-ex1,2.0);
  }
  ex2 = ex2/(n[0]-1.);
  
  double a0=0.,a1=0., b0=0., b1=0., eps=0.00000001;
  int Iter=20;
  double F1, F2, a,b,c;
  i=0;
  a1 = ex1*(ex1*(1.-ex1)/ex2-1.);
  b1 = (1-ex1)*(ex1*(1.-ex1)/ex2-1.);
  while(i<Iter && fabs(a1-a0)>eps && fabs(b1-b0)>eps){
    a0 = a1; b0=b1;
    b = -trigamma(a1+b1);
    a = trigamma(a1) +b;
    c = trigamma(b1)+b;
    F1 = digamma(a1)-digamma(a1+b1)-elx1;
    F2 = digamma(b1)-digamma(a1+b1)-elx2;
    a1 = a1 - (c*F1-b*F2)/(a*c-pow(b,2.));
    b1 = b1 - (-b*F1+a*F2)/(a*c-pow(b,2.));
    i=i+1;
  }
  para[0] = a1;
  para[1] = b1;
}


void FitWeibull(double *x, int *n, double *para){
  int i,j;
  double x0=1.,x1=1.5, eps=0.00000001;
  double g0x, g1x, a1,a2,a3,a4,c1,c2;
  int Iter=20;
  i=0;
  while(i<Iter && fabs(x1-x0)>eps){
    x0 = x1;
    a1=0.;a2=0.;a3=0.;a4=0.;
    c1=0.;c2=0.;
    for(j=0;j<n[0];j++){
      c1 = pow(x[j],x1);
      c2 = log(x[j]);
      a1 += c2*c1;
      a2 += c1;
      a3 += c2;
      a4 += pow(c2,2.)*c1;
    }
    g0x= 1./x1-a1/a2+a3/n[0];
    g1x = -1./pow(x1,2.) - (a4*a2 - pow(a1,2.))/pow(a2,2.);
    x1 = x1 - g0x/g1x;
    i=i+1;
  }

  x0 = pow(a2/n[0],1./x1);
    
  para[0] = x1;
  para[1] = x0;
}



// ISNA(x): true for R's NA only
// ISNAN(x): true for R's NA and IEEE NaN
// R_FINITE(x): false for Inf,-Inf,NA,NaN
// R_IsNaN(x): true for NaN but not NA

// Rprintf:  printing from a C routine compiled into R

double g0weibull(double k, double ratio)
{
  return(exp(2.0*lgamma(1.+1./k)) - ratio * exp(lgamma(1+2.0/k)));
}

double g1weibull(double k, double ratio)
{
  double g1a,g1b,g2a,g2b;
  g1a = digamma(1. + 1./k);
  g1b = exp(lgamma(1. + 1./k));
  g2a = digamma(1. + 2./k);
  g2b = exp(lgamma(1. + 2./k));
  return(2.0/(k*k)*(ratio*g2a*g2b-g1a*g1b*g1b));
}

double mmeWeibull(double ratio){
  int i=0,Iter=1000;
  double x0 = 1.0, delta=1.0, tol=1.0e-10;
  while(i < Iter && fabs(delta) > tol){
    delta =  g0weibull(x0,ratio)/ g1weibull(x0,ratio);
    x0 -= delta;
    i++;
  }
  return(x0);
}


void mme(double x[], double counts[], double weights[],
	 int n, int idist, double par[])
{
  int i;
  double nsum = 0.0, xbar=0., sig=0., sig2=0.;
  for(i=0;i<n;i++){
    xbar += x[i] * counts[i] * weights[i];
    nsum += counts[i] * weights[i];
  }
  xbar /= nsum;
  for(i=0;i<n;i++){
    sig2 += pow(x[i]-xbar, 2.0) * counts[i] * weights[i];
  }
  sig2 /= nsum - 1.;
  sig = sqrt(sig2);
  

  switch(idist)
    {
    case 0 :
      par[0] = xbar; par[1] = sig;
      break;
    case 1 :
      sig = xbar*(1.-xbar)/sig2-1.0; // an estimate for alpha+beta
      par[0] = xbar * sig;  // mme for shape parameter
      par[1] = (1-xbar) * sig; // mme for scale parameter
      //     Rprintf("mean(x)=%f,var(x)=%f\n",xbar,sig);
      break;
    case 2 :
      par[0] = pow(xbar,2.0)/sig2;  // mme for shape parameter
      par[1] = sig2/xbar; // mme for scale parameter
      break;
    case 3 :
      sig = xbar * xbar;
      par[0] = mmeWeibull(sig/(sig+sig2));
      sig = exp(lgamma(1.+1./par[0]));
      par[1] = xbar/sig;
      break;
    default :
      Rprintf("Algorithm to be developed...");
    }
}

static double funllk(int npar, double *pars, void *ex)
{
  double *x = (double*)ex, tol = 6.123234e-17;
  int i,n = (int)x[0], idist = (int)x[1];

  double dFx=0.0,res = 0.0, t1, t2;
  for(i=0;i<n;i++){
    switch(idist)
      {
      case 0 :
	if(pars[1] > tol ){
	  t2 = pnorm(x[i+2] + x[i+2+2*n], pars[0], pars[1], 1, 0);
	  t1 = pnorm(x[i+2] - x[i+2+2*n], pars[0], pars[1], 1, 0);
	  dFx = t2 - t1;
	}else{
	  dFx = 2.0 * tol;
	}
	break;
      case 1 :
	if(pars[1] > 0 && pars[0]>0){
	  t2 = pbeta(x[i+2] + x[i+2+2*n], pars[0], pars[1], 1, 0);
	  t1 = pbeta(x[i+2] - x[i+2+2*n], pars[0], pars[1], 1, 0);
	  dFx = t2 - t1;
	}else{
	  dFx = 2.0 * tol;
	}
	break;
      case 2 :
	if(pars[1] > 0 && pars[0]>0){
	  t2 = pgamma(x[i+2] + x[i+2+2*n], pars[0], pars[1], 1, 0);
	  t1 = pgamma(x[i+2] - x[i+2+2*n], pars[0], pars[1], 1, 0);
	  dFx = t2 - t1;
	}else{
	  dFx = 2.0 * tol;
	}
	break;
      case 3 :
	if(pars[1] > 0 && pars[0]>0){
	  t2 = pweibull(x[i+2] + x[i+2+2*n], pars[0], pars[1], 1, 0);
	  t1 = pweibull(x[i+2] - x[i+2+2*n], pars[0], pars[1], 1, 0);
	  dFx = t2 - t1;
	}else{
	  dFx = 2.0 * tol;
	}
	break;
      default :
	Rprintf("Algorithm to be developed...");
      }
    if(dFx > tol)
      res += x[i+2+3*n] * x[i+2+n] * log(dFx);
    else
      res += x[i+2+3*n] * x[i+2+n] * log(tol);
  }
  return(-res);
}


 /* The NM-algorithm:
    
    The Nelder–Mead method or downhill simplex method or amoeba method
    is used to find the mle estimates numerically.

    Data could be weighted for EM-algorithms, or all with 1 for
    univariate EM with one component only.

    The current version is applicable for two-parameter distribution
    families only.  For higher dimensional distributions, the results
    might not be good, and very time-consuming as well.
 */
void NMMle(double x[], double counts[], double widths[], double weights[],
	   int n, int idist, double par[])
{
  double abstol=0.000000001,reltol=0.0000000001,val;
  int npar=2,ifail=0,trace=0, maxit=500, fncount;
  double alpha=1.0, beta=0.5, gamma=2;
  double dpar[npar],opar[npar];

  int i;
  double y[4*n+2];//push the data into y:1=size, 2=idist
  y[0] = 1. * n; y[1] = 1. * idist;
  for(i=0;i<n;i++){
    y[i+2] = x[i];
    y[i+2+n] = counts[i];
    y[i+2+2*n] = widths[i];
    y[i+2+3*n] = weights[i];
  }
  /* Estimate the initial values using method of moments (store in
     dpar[]) */

  mme(x,counts,weights,n,idist,dpar);

  //  GetRNGstate();
  nmmin(npar, dpar, opar, &val, funllk, &ifail, abstol, reltol, 
	(void *)y, alpha, beta, gamma, trace, &fncount, maxit);
  //  PutRNGstate();
  //  output the fitted MLE in opar
  //  Rprintf("par[0]=%f,par[1]=%f\n",opar[0],opar[1]);
  par[0]=1.0;
  for(i=0;i<2;i++) par[i+1] = opar[i];
}

/*
  Create on 11/26/2011

  EM(x, counts, widths, weights, nb, nc, idist[0], par,llk);

  For given inistial mean values, applied EM algorithm to estimate the
  proportions (mixing parameters) and the parameters of each
  component.


 */

double Prob(double x, double h, double shape, double scale, int idist)
{
  double t1, t2, dFx=0.0, tol = 6.123234e-17;
  switch(idist)
    {
    case 0 :
      if(scale > 0){
	t2 = pnorm(x + h, shape, scale, 1, 0);
	t1 = pnorm(x - h, shape, scale, 1, 0);
	dFx = t2 - t1;
      }else{
	dFx = 2.0 * tol;
      }
      break;
    case 1 :
      if(scale > 0 && shape>0){
	t2 = pbeta(x + h, shape, scale, 1, 0);
	t1 = pbeta(x - h, shape, scale, 1, 0);
	dFx = t2 - t1;
      }else{
	dFx = 2.0 * tol;
      }
      break;
    case 2 :
      if(scale > 0 && shape>0){
	t2 = pgamma(x + h, shape, scale, 1, 0);
	t1 = pgamma(x - h, shape, scale, 1, 0);
	dFx = t2 - t1;
      }else{
	dFx = 2.0 * tol;
      }
      break;
    case 3 :
      if(scale > 0 && shape>0){
	t2 = pweibull(x + h, shape, scale, 1, 0);
	t1 = pweibull(x - h, shape, scale, 1, 0);
	dFx = t2 - t1;
      }else{
	dFx = 2.0 * tol;
      }
      break;
    default :
      Rprintf("Algorithm to be developed...");
    }
  return(dFx);
}

double Dens(double x, double shape, double scale, int idist)
{
  double out=0.0, tol = 6.123234e-17;
  switch(idist)
    {
    case 0 :
      if(scale > 0)
	out = dnorm(x, shape, scale, 0);
      else out = 2.0 * tol;
      break;
    case 1 :
      if(scale > 0 && shape>0)
	out = dbeta(x, shape, scale, 0);
      else out = 2.0 * tol;
      break;
    case 2 :
      if(scale > 0 && shape>0)
	out = dgamma(x, shape, scale, 0);
      else out = 2.0 * tol;
      break;
    case 3 :
      if(scale > 0 && shape>0)
	out = dweibull(x, shape, scale, 0);
      else out = 2.0 * tol;
      break;
    default :
      Rprintf("Algorithm to be developed...");
    }
  return(out);
}

void EM(double x[], double counts[], double widths[], double weights[],
	int n, int m, int idist, double par[], double llk[])
{
  double mu[m], s[m], p[m], tau[n][m];
  int i, j, iter;
  double tol=0.00001, tol0 = 6.123234e-17;
  double dllk=1.0, llk2=100.,llk1=1., nsum=0.0;
  double t1, t2, tpar[2],tpar3[3];

  // compute total sample size.
  for(i = 0; i < n; i++){
    nsum += counts[i];
  }
  
  /*
    We first use the given initial values or randomly selected
    component centers to approximate the weights of each observation
    on each component.  We use the normal distribution for this
    purpuse.
   */
  mme(x,counts,weights,n,0,tpar);//get sample mean and std.dev.
  for(i = 0; i < m; i++){
    s[i] = tpar[1]; // use common sample std.dev.
    mu[i] = par[i*3+1]; //and centers passed from main program.
    p[i] = 1./m; //wild guess (to be updated).
    //    Rprintf("p=%f,mu=%f, sigma=%f\n",p[i],mu[i],s[i]);
  }

  for(i = 0; i < n; i++){
    t1 = 0.0;
    for(j = 0; j < m; j++){
      t2 = Dens(x[i],mu[j],s[j], 0); //normal densities
      tau[i][j] = p[j] * t2;
      t1 += tau[i][j];
    }
    if(t1 > tol0) 
      for(j=0;j<m;j++)
	tau[i][j] /= t1;//for each distinct x value
    else //in case zero proportion for xi for all components 
      for(j=0;j<m;j++) tau[i][j] = 1.0/m;//for each distinct x value
  }
  for(j=0;j<m;j++){
    t1 = 0.0; t2 = 0.0; 
    for(i=0;i<n;i++){
      weights[i] = tau[i][j];
      t1 = tau[i][j] * counts[i];
      t2 += t1; //sum of weights used as n
    }
    p[j] = t2/nsum;
    // update the jth compoent parameters using NM-alg
    NMMle(x, counts, widths, weights, n, idist, tpar3);
    mu[j] = tpar3[1]; s[j] = tpar3[2];
   
  }

  iter = 0;
  while(iter < 1000 && dllk > tol){
    iter++;
    llk1 = llk2; //save the old LLK to llk1

    for(i = 0; i < n; i++){
      t1 = 0.0;
      for(j = 0; j < m; j++){
	t2 = Dens(x[i],mu[j],s[j], idist); //using the specified component type
	tau[i][j] = p[j] * t2;
	t1 += tau[i][j];
      }
      if(t1 > tol0) 
	for(j=0;j<m;j++)
	  tau[i][j] /= t1;//for each distinct x value
      else //in case zero proportion for xi for all components 
	for(j=0;j<m;j++) tau[i][j] = 0.0;//for each distinct x value
    }
    for(j=0;j<m;j++){
      t1 = 0.0; t2 = 0.0; 
      for(i=0;i<n;i++){
	weights[i] = tau[i][j];
	t1 = tau[i][j] * counts[i];
	t2 += t1; //sum of weights used as n
      }
      p[j] = t2/nsum;
      // update the jth compoent parameters using NM-alg
      NMMle(x, counts, widths, weights, n, idist, tpar3);
      mu[j] = tpar3[1]; s[j] = tpar3[2];
    }
    //Update the lok-likrlihood
    llk2=0.0;
    for(i=0;i<n;i++){
      t1 = 0.0;
      for(j=0;j<m;j++)
	t1 += p[j] * Dens(x[i],mu[j],s[j], idist);
      if(t1 > tol0) llk2 +=  counts[i] * log(t1);
      else llk2 +=  counts[i] * log(tol0);
    }
    dllk = fabs(llk2 - llk1);
  }
  for(i=0; i<m;i++){
    par[i*3] = p[i];
    par[i*3+1] = mu[i];
    par[i*3+2] = s[i];
  }
  llk[0] = llk2;
}


/*
  2011/11/25: The main function of EM-algorithm.  

  Inputs: x:= centers of bins; counts := the counts of observations in
  the bins; and width := the bin widths.
 
  Outputs: ncomp := number of mixture components; par := estimated
  proportions and lambdas, gammas; llk := the (likelihood, AIC, BIC,
  AICc);

  More details: we specify the number of components and randomly
  choose initial values and select the best model based on AIC or BIC.

  AIC = -2*ln(likelihood) + 2*K,
  AICc = AIC + 2K(K+1)/(N-K-1)
  BIC = -2*ln(likelihood) + ln(N)*K 

  where:

  K = model degrees of freedom (the number of model parameters);
  N = number of observations.

  Generally speaking, a better model has smaller AIC/BIC.  We compares
  the likelihood for model selection when the number of components is
  fixed and N is fixed.


  */


void fitmm(double *x, double *counts, double *widths, int *bins, 
	   int *idist, int *comp, double *par, double *llk)
{
  int i,  nb = bins[0], nc = comp[0]; // nc := number of components;
  double weights[nb],llk0, t1, tol0 = 6.123234e-17;
  for(i = 0; i < nb; i++){
    widths[i] /= 2.0; //half width
    weights[i] = 1.0;
  }
  if(nc >1){ //if number of component is two or more
    EM(x, counts, widths, weights, nb, nc, idist[0], par,llk);
    //GetRNGstate();
    //PutRNGstate();
    //    llk0 = llk[0];
  }else{
    NMMle(x, counts, widths, weights, nb, idist[0], par);
    llk0 = 0.0;
    for(i = 0; i < nb; i++){
      t1 =  Dens(x[i],par[1], par[2], idist[0]);
      if(t1 > tol0) llk0 += counts[i] * log(t1);
      else llk0 += counts[i] * log(tol0);
    }
    llk[0] = llk0;
  }
}

 
/*
  To perform the EM algorithm for data with rounding error.  Normal
  components are assumed.
  */

void renormmle(double *w, double *f, double *b, int n, double mu, double sig)
{
  int i,iter0;
  double t1,t2,t3,t4,t5,mu1,sig1,g0,g1,g2,f0,f1,f2;
  double za[n],zb[n],pdfa[n],pdfb[n],cdfa[n],cdfb[n];
  double re1,re2, tol=0.00001,tol0 = 6.123234e-17;
  iter0=0;
  re1 = 1.0;
  re2 = 1.0;

  while(iter0 < 10000 && re1 > tol){
    f0 = 0.0;
    f1 = 0.0;
    f2 = 0.0;
    g0 = 0.0;
    g1 = 0.0;
    g2 = 0.0;
    for(i=0;i<n;i++){
      za[i] = (w[i]-b[i]-mu)/sig;
      zb[i] = (w[i]+b[i]-mu)/sig;
      pdfa[i]=dnorm(za[i],0.0,1.0,0);
      pdfb[i]=dnorm(zb[i],0.0,1.0,0);
      cdfa[i]=pnorm(za[i],0.0,1.0,1,0);
      cdfb[i]=pnorm(zb[i],0.0,1.0,1,0);
      t1 = pdfb[i]-pdfa[i];
      t2 = cdfb[i]-cdfa[i];
      t3 = zb[i]*pdfb[i]-za[i]*pdfa[i];
      t4 = zb[i]*zb[i]*pdfb[i]-za[i]*za[i]*pdfa[i];
      t5 = zb[i]*zb[i]*zb[i]*pdfb[i]-za[i]*za[i]*za[i]*pdfa[i];
      if(fabs(t2) > tol0){
	f0 = f0 + f[i]*t1/t2;
	g0 = g0 + f[i]*sig * t3/t2;
	f1 = f1 + f[i]*(t3*t2+t1*t1)/t2/t2;
	f2 = f2 + f[i]*(t4*t2+t3*t1)/t2/t2;
	g1 = g1 + f[i]*((t4-t1)*t2+t3*t1)/t2/t2;
	g2 = g2 + f[i]*(t5*t2+t3*t3)/t2/t2;
      }
    }
    f1 = f1/sig;
    f2 = f2/sig;
    mu1 = mu;
    sig1 = sig;
    mu = mu1 - (f0*g2-f2*g0)/(f1*g2-f2*g1);
    sig = sig1 - (f1*g0-f0*g1)/(f1*g2-f2*g1);
    t1 = fabs((mu-mu1)/fmin(mu,mu1));
    t2 = fabs((sig-sig1)/fmin(sig,sig1));
    t3 = fabs((mu-mu1));
    t4 = fabs((sig-sig1));
    re1 = fmax(t1,t3);
    re2 = fmax(t2,t4);
	  re1 = fmax(re1,re2);
    iter0++;
  }
}

void reemnorm(double *x, double *f, double *b, int *size, int *Iter, 
	      double *p, double *mu, double *s, double *llk)
{
  /*
    2011/06/23: for em with rounding errors.  
    The data are grouped: f is the frequency and b is half bin width 
   */
  //Iter is reused: pass the number of components and return the iteration number 
  int n=size[0], m = Iter[0],i,j,iter=0;
  //int k; // use k to record those with very small s so the distr reduce to discrete case
  double slim=1./12.;
  double tau[n][m],w[n],nsum=0.0;
  double t1,t2,tol=0.00001,llk2=100.,llk1=1.,dllk=1.,mllk=1.,tol0 = 6.123234e-17;

  for(i=0;i<n;i++){
    nsum += f[i]; // to compute the sample size of the raw data
  }

  while(iter < 10000 && dllk > tol && dllk/mllk > tol){
      llk1 = llk2; //save the old LLK to llk
      //E-step: given mu,sig and p, compute tau
      for(i=0;i<n;i++){
	t1 = 0.0;
	for(j=0;j<m;j++){
	  if(b[i]<=0.0){  //if no rounding error
	    tau[i][j] = p[j]*dnorm(x[i],mu[j],s[j],0);
	  }else{
	    tau[i][j] = p[j]*(pnorm(x[i]+b[i],mu[j],s[j],1,0)-pnorm(x[i]-b[i],mu[j],s[j],1,0));
	  }
	  t1 += tau[i][j];
	}
	if(t1 > tol0){
	  for(j=0;j<m;j++){
	    tau[i][j] /= t1;//for each distinct x value
	  }
	}else{
	  for(j=0;j<m;j++){
	    tau[i][j] = 0.0;//for each distinct x value
	  }
	}
      }
      
      //M-step: compute p, mu (data are weighted with f)
      for(j=0;j<m;j++){
	t1 = 0.0;
	t2 = 0.0;
	for(i=0;i<n;i++){
	  t1 += tau[i][j]*f[i]; 
	  t2 += tau[i][j] * x[i]*f[i];
	}
	p[j] = t1/nsum;
	mu[j] = t2/t1;
	s[j] = t1;
      }
      // further update sig with a current mu
      for(j=0;j<m;j++){
	t1 = 0.0;
	for(i=0;i<n;i++){
	  t1 += tau[i][j]*pow(x[i]-mu[j],2.0)*f[i]; 
	  w[i] = f[i]*tau[i][j]; // new weight w=f*tau
	}
	t2 = t1/s[j];
	if(t2<slim) t2 = slim;	  
	if(t2 > 2.*slim) t2 = t2-slim;	  
	s[j] = sqrt(t2);
	t1 = mu[j]; // store initial values
	t2 = s[j];
	renormmle(x,w,b,n,mu[j],s[j]);
	if(isnan(mu[j])||isnan(s[j])){
	  mu[j] = t1; s[j] = t2;
	}

      }
      
      //Update the lok-likrlihood
      llk2=0.0;
      for(i=0;i<n;i++){
	t1 = 0.0;
	for(j=0;j<m;j++){
	  t1 += p[j]*dnorm(x[i],mu[j],s[j],0)*f[i];
	}
	if(t1>0.0) llk2 += log(t1); 
      }

      dllk = fabs(llk2-llk1);
      mllk = fmax(fabs(llk1),fabs(llk2));
      iter++;
    }
  llk[0] = llk2;
  Iter[0] = iter;
}

