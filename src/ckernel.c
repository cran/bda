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

void wlinbin(double *X, double *W, int *size,
	     double *xa, double *xb,
	     int *Msize, int *trun, double *gcnts) {
  double lxi, delta, rem,a=xa[0],b=xb[0];
    int i, li, n=size[0], M=Msize[0];

    // Initialize grid counts to zero
    for (i = 0; i < M; i++) {
        gcnts[i] = 0.0;
    }

    delta = (b - a) / (M - 1.0);
    for (i = 0; i < n; i++) {
        lxi = ((X[i] - a) / delta) + 1;

        // Find integer part of "lxi"
        li = (int) lxi;

        rem = lxi - li;
        if (li >= 1 && li < M) {
            gcnts[li - 1] += (1 - rem) * W[i];
            gcnts[li] += rem * W[i];
        }

        if (li < 1 && trun == 0) {
            gcnts[0] += W[i];
        }

        if (li >= M && trun == 0) {
            gcnts[M - 1] += W[i];
        }
    }
}

void probin(double *X, double *Y, int *size,
	    double *xa, double *xb, int *Msize,
	    double *gcnts) {
  double lxi, lyi, delta, rem, dx,n=size[0],a=xa[0],b=xb[0];
  int i, li, lj,M=Msize[0];

    // Initialize grid counts to zero
    for (i = 0; i < M; i++) {
        gcnts[i] = 0.0;
    }

    delta = (b - a) / M;
    for (i = 0; i < n; i++) {
        lxi = ((X[i] - a) / delta) + 1;
        lyi = ((Y[i] - a) / delta) + 1;
        dx = (Y[i] - X[i]) / delta;

        // Find integer part of "lxi" and "lyi"
        li = (int) lxi;
        lj = (int) lyi;

        if (li == lj) {
            gcnts[li - 1] += 1;
        } else {
            for (int j = li - 1; j < lj; j++) {
                if (j < lj - 1 && j == li - 1) {
                    rem = lxi - li;
                    gcnts[j] += (1 - rem) / dx;
                } else if (j < lj - 1 && j > li - 1) {
                    gcnts[j] += 1;
                } else {
                    rem = lyi - lj;
                    gcnts[j] += rem / dx;
                }
            }
        }
    }
}

void yldist(double *gcnts, int *Msize, double *Y) {
    double theta, sumcos, sumsin;
    int l, k, M=Msize[0];

    // Initialize Y_ell to zero
    for (l = 0; l < M/2; l++) {
        Y[l] = 0.0;
    }

    for (l = 0; l < M/2; l++) {
        sumcos = 0.0;
        sumsin = 0.0;
        for (k = 0; k < M; k++) {
            theta = 2.0 * M_PI * (k * 1.0) * l / M;
            sumcos += gcnts[k] * cos(theta);
            sumsin += gcnts[k] * sin(theta);
        }
        Y[l] = (sumcos * sumcos + sumsin * sumsin) / (M * M);
    }
}

void linbin(double *X, int *size, double *xa, double *xb,
	    int *Msize, int *trun, double *gcnts) {
  double lxi, delta, rem,a=xa[0],b=xb[0];
  int i, li,n=size[0],M=Msize[0];

    // Initialize grid counts to zero
    for (i = 0; i < M; i++) {
        gcnts[i] = 0.0;
    }

    delta = (b - a) / (M - 1);
    for (i = 0; i < n; i++) {
        lxi = ((X[i] - a) / delta) + 1;

        // Find integer part of "lxi"
        li = (int) lxi;

        rem = lxi - li;
        if (li >= 1 && li < M) {
            gcnts[li - 1] += (1 - rem);
            gcnts[li] += rem;
        }

        if (li < 1 && trun[0] == 0) {
            gcnts[0] += 1;
        }

        if (li >= M && trun[0] == 0) {
            gcnts[M - 1] += 1;
        }
    }
}

void ofcpdf(double *y, double *f, double *a, double *b,
	    int *my, double *x, int *mx, double *bw) {
  double nsum, tmp, t1, t2, sqrt2h, sqrt2,h=bw[0];
  int i, j,ny=my[0],nx=mx[0];
  
    sqrt2 = sqrt(2.0);
    sqrt2h = sqrt2 * h;

    nsum = 0.0;
    for (i = 0; i < ny; i++) {
        nsum += f[i];
    }

    for (i = 0; i < nx; i++) {
        tmp = 0.0;
        for (j = 0; j < ny; j++) {
            t1 = (b[j] + y[j] - x[i]) / sqrt2h;
            t2 = (x[i] - a[j] - y[j]) / sqrt2h;
            tmp += 0.5 * f[j] * (erf(t1 / sqrt2) + erf(t2 / sqrt2)) / (b[j] - a[j]);
        }
        x[i] = tmp / nsum;
    }
}


void remp(int *n2, double *y, double *f,
	  double *a, double *b,
	  int *n1, double *Fx, double *x,
	  double *u) {

  int i, j, k, l, icounter;
  double F0, F1, y0, dy, t1;
  int ny=n2[0],nx=n1[0];
  
    icounter = 0;
    k = 1;
    for (i = 0; i < ny; i++) {
        // Search for F0
        y0 = y[i] + a[i];
        dy = x[nx - 1] - x[0];
        for (j = 0; j < nx; j++) {
            t1 = fabs(x[j] - y0);
            if (t1 < dy) {
                dy = t1;
                k = j;
            }
        }
        F0 = Fx[k];

        // Search for F1
        y0 = y[i] + b[i];
        dy = x[nx - 1] - x[0];
        for (j = 0; j < nx; j++) {
            t1 = fabs(x[j] - y0);
            if (t1 < dy) {
                dy = t1;
                k = j;
            }
        }
        F1 = Fx[k];

        // Temporary using y0 to keep information
        for (l = 0; l < (int)f[i]; l++) {
            icounter++;
            y0 = u[icounter] * F1 + (1.0 - u[icounter]) * F0;
            dy = 1.0;
            for (j = 0; j < nx; j++) {
                t1 = fabs(Fx[j] - y0);
                if (t1 < dy) {
                    dy = t1;
                    k = j;
                }
            }
            u[icounter] = x[k];
        }
    }
}


void mlensimp(double *w, double *f, double *a, double *b,
	      int *size, double *theta) {
  int i, iter,n=size[0];
    double t2, t3, t4, t5, mu, sig, sig1, g0, g2, re2;
    double *za, *zb, *pdfa, *pdfb, *cdfa, *cdfb;

    // Allocate memory for arrays
    za = (double *)malloc(n * sizeof(double));
    zb = (double *)malloc(n * sizeof(double));
    pdfa = (double *)malloc(n * sizeof(double));
    pdfb = (double *)malloc(n * sizeof(double));
    cdfa = (double *)malloc(n * sizeof(double));
    cdfb = (double *)malloc(n * sizeof(double));

    // Set initial values of theta
    mu = theta[0];
    sig = theta[1];

    iter = 0;
    re2 = 1.0;
    g0 = 0.0;
    g2 = 0.0;

    while (iter < 10000 && re2 > 0.000001) {
        for (i = 0; i < n; i++) {
            za[i] = (w[i] + a[i] - mu) / sig;
            zb[i] = (w[i] + b[i] - mu) / sig;
            pdfa[i] = dnorm(za[i],0.0,1.0,0);
            pdfb[i] = dnorm(zb[i],0.0,1.0,0);
            cdfa[i] = pnorm(za[i],0.0,1.0,1,0);
            cdfb[i] = pnorm(zb[i],0.0,1.0,1,0);
            t2 = cdfb[i] - cdfa[i];
            t3 = zb[i] * pdfb[i] - za[i] * pdfa[i];
            t5 = zb[i] * zb[i] * zb[i] * pdfb[i] - za[i] * za[i] * za[i] * pdfa[i];
            g0 += f[i] * sig * t3 / t2;
            g2 += f[i] * (t5 * t2 + t3 * t3) / (t2 * t2);
        }

        sig1 = sig;
        sig = sig1 - g0 / g2;
        t2 = fabs((sig - sig1) / fmin(sig, sig1));
        t4 = fabs((sig - sig1));
        re2 = fmax(t2, t4);
        iter++;
    }

    theta[1] = sig;
    n = iter;

    // Free allocated memory
    free(za);
    free(zb);
    free(pdfa);
    free(pdfb);
    free(cdfa);
    free(cdfb);
}


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

void bin2d(double *x, double *y, int *size,
	   double *brk1, int *nbrk1,
	   double *brk2, int *nbrk2,
	   double *cnt)
{
  int i, j, ix=0, iy=0, k=0,nk;
  double xi,yi;
  //initialize parameters
  nk = (nbrk1[0]-1) * (nbrk2[0]-1);
  for(i=0; i < nk; i++){
    cnt[i] = 0.0;
  }
  
  for(i=0; i < size[0]; i++){
    xi = x[i];
    yi = y[i];
    ix = nbrk1[0] - 1;
    for(j=1; j < nbrk1[0]; j++){
      if(xi < brk1[j]){
	ix = j;
	break;
      }
    }
    iy = nbrk2[0] - 1;
    for(j=1; j < nbrk2[0]; j++){
      if(yi < brk2[j]){
	iy = j;
	break;
      }
    }
    k = (ix - 1)*(nbrk2[0]-1) + (iy - 1);
    if(k >= nk) k = nk - 1;
    cnt[k] += 1.0;
  }
}

void chgpt(double *x, double *y, int *size, int *k, double *ans)
{
  /* divide the data (x,y) into two parts (1:k) and (k+1):n; fitting
     two least square straight lines, respectively, and compute the
     residual sum of squares (RSS). Changing k to find a split with
     minimum RSS.

     LSE for y=a+bx, b = (n*SXY - SX*SY)/(n*SXX - SX^2), a=(SY-b*SX)/n
   */
  int i, j,n1,n2;
  double sx1,sy1,sxx1,sxy1,a1,b1;
  double sx2,sy2,sxx2,sxy2,a2,b2;
  double rss1,rss2,rssmin,rsssum;
  
  //initialize parameters
  sx1  = 0.0;
  sy1  = 0.0;
  sxx1 = 0.0;
  sxy1 = 0.0;
  n1 = k[0];
  
  sx2  = 0.0;
  sy2  = 0.0;
  sxx2 = 0.0;
  sxy2 = 0.0;
  n2 = size[0] - n1;

  for(i=0; i < n1; i++){
    sx1 += x[i];
    sy1 += y[i];
    sxx1 += x[i]*x[i];
    sxy1 += x[i]*y[i];
  }
  b1 = (n1*sxy1 - sx1*sy1)/(n1*sxx1-sx1*sx1);
  a1 = (sy1-b1*sx1)/n1;
  
  rss1 = 0.0;
  for(i=0; i < n1; i++){
    rss1 += pow(y[i] - a1 - b1*x[i], 2.0);
  }

  
  for(i=n1; i < size[0]; i++){
    sx2 += x[i];
    sy2 += y[i];
    sxx2 += x[i]*x[i];
    sxy2 += x[i]*y[i];
  }
  b2 = (n2*sxy2 - sx2*sy2)/(n2*sxx2-sx2*sx2);
  a2 = (sy2-b2*sx2)/n2;

  rss2 = 0.0;
  for(i=n1; i < size[0]; i++){
    rss2 += pow(y[i] - a2 - b2*x[i], 2.0);
  }

  rssmin = rss1 + rss2;
  ans[0] = x[n1-1];

  int nstart=n1,nend=size[0]-n1;
  for(i=nstart; i < nend; i++){
    n1++;
    n2--;
    sx1 += x[i];
    sy1 += y[i];
    sxx1 += x[i]*x[i];
    sxy1 += x[i]*y[i];
    sx2 -= x[i];
    sy2 -= y[i];
    sxx2 -= x[i]*x[i];
    sxy2 -= x[i]*y[i];
    b1 = (n1*sxy1 - sx1*sy1)/(n1*sxx1-sx1*sx1);
    a1 = (sy1-b1*sx1)/n1;
    rss1 = 0.0;
    for(j=0; j < i+1; j++){
      rss1 += pow(y[j] - a1 - b1*x[j], 2.0);
    }
    b2 = (n2*sxy2 - sx2*sy2)/(n2*sxx2-sx2*sx2);
    a2 = (sy2-b2*sx2)/n2;
    rss2 = 0.0;
    for(j=i+1; j < size[0]; j++){
      rss2 += pow(y[j] - a2 - b2*x[j], 2.0);
    }
    rsssum = rss1 + rss2;
    if(rsssum < rssmin){
      rssmin = rsssum;
      ans[0] = x[n1-1];
    }
  }
}
