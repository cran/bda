#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"

typedef double (*Funx)(double, double, double, double, int);


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

    //    Rprintf("\np=%f,mu=%f, sigma=%f\n",tpar3[0],tpar3[1],tpar3[2]);
    
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

 
