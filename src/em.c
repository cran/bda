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

void em3(int *size, double *x, double *pars, double *tol)
{
  int i, j, k, n = size[0], iter = 50000;
  double mu[3], s[3], p[3], w[n][3];;
  double mu2[3], s2[3], p2[3];;
  double wsum[3], fsum[n], delta;
  //initialize parameters
  for(j=0; j<n; j++){
    for(i=0; i<3; i++){
      w[j][i] = 0.0;
    }
    fsum[j] = 0.0;
  }
  mu[0] = pars[0]; mu[1] = pars[1]; mu[2] = pars[2];
  s[0] = pars[3]; s[1] = pars[4]; s[2] = pars[5];
  p[0] = pars[6]; p[1] = pars[7]; p[2] = 1.0-p[0]-p[1];
  /* EM-algorithm */
  
  for(i=0; i<iter; i++){
    delta = 0.0;
    for(k=0; k<3; k++){
      mu2[k] = mu[k];
      s2[k] = s[k];
      p2[k] = p[k];
      wsum[k] = 0.0;
    }
    /*E-step */
    for(j=0; j<n; j++){
      fsum[j] = 0.0;
      for(k=0; k<3; k++){
	w[j][k] = p[k] * dnorm(x[j],mu[k],s[k],0);
	fsum[j] += w[j][k];
      }
      
      for(k=0; k<3; k++){
	w[j][k] /= fsum[j];
	wsum[k] += w[j][k];
      }
    }
    /* M-step */
    for(k=0; k<3; k++){ //update proportions
      p[k] = wsum[k]/n;
      delta += fabs(p2[k] - p[k]);
      mu[k] = 0.0;
      s[k] = 0.0;
    }
    for(j=0; j<n; j++){ //update means
      for(k=0; k<3; k++){
	mu[k] += w[j][k] * x[j]; 
      }
    }
    for(k=0; k<3; k++){
      mu[k] /= wsum[k];
      delta += fabs(mu2[k] - mu[k]);
    }
    
    for(j=0; j<n; j++){ //update SD
      for(k=0; k<3; k++){
	s[k] += w[j][k] * (x[j]-mu[k]) * (x[j]-mu[k]); 
      }
    }
    for(k=0; k<3; k++){
      s[k] = sqrt(s[k]/wsum[k]);
      delta += fabs(s2[k] - s[k]);
    }
    if(delta < tol[0]) break;
  }

  /* reuse variable size to return the iterations */
  size[0] = i;
  pars[0] = mu[0];
  pars[1] = mu[1];
  pars[2] = mu[2];
  pars[3] = s[0];
  pars[4] = s[1];
  pars[5] = s[2];
  pars[6] = p[0];
  pars[7] = p[1];
}
