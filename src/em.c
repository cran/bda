#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <stdio.h>
#include <math.h>
#include "R_ext/Applic.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


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

