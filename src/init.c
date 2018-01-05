/*
 *  Copyright (C) 2009-2010  B. Wang
 *  Unlimited use and distribution (see LICENCE).
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


void F77_SUB(linbin)(double *x, int *n, double *a, double *b,
		     int *m, int *trun, double *gcounts);

void F77_SUB(lbtwod)(double *x, int *n, double *a1,
		     double *a2, double *b1, double *b2, int *m1,
		     int *m2, double *gcounts);

void F77_SUB(bintwod)(double *x, int *n, double *g1,
		      double *g2, int *m1,
		      int *m2, double *gcounts);

void F77_SUB(iterfx)(double *fx, double x0, int *n, double *x,
		     double *f, int m,
		     double *w, double *h, double *iter);
		    
void rootGldFmklBisection(double *q, double *lambdas);
void KSPvalue(double *x0);

// lprsmooth.c
void tubecv(double *kappa, double *level);

void DesignMatrix(double *xv, int *size, double *bw, double *dm);
		    
void lpsmooth(double *xgrid, int *ngrid, double *x, double *y, 
	      int *size, double *bw, int *lscv, 
	      double *range, int *adaptive, double *ellx, 
	      double *kappa);
//AS254.f
void F77_SUB(emmix)(int *y, int *ny, int *ng,
		    double *X0, double *X,
		    double *P, double *XBAR, double *VAR,
		    double *xlogl, int *wk, 
		    int *itrunc, int *nl, int *nu);

		   //smoothkde.f
void F77_SUB(smoothkde)(double *fx, double x0, int *n, 
			double *x,  double *f, int m,
			double *w, double *h, int *iter);

//  KernelWKDE.c  
void BDMLE(double *f, double *a, double *b, int *nbin,
	   double *pars, int *npar, int *dist);
void bdregmle(double *F, double *x,double *freq, 
	      int *nu, int *n, int *dist, double *pars);

//  cbootkde.c/fbootkde.f
void hbmise(double *x, double *f, double *h, int *n, double *hopt);
void F77_SUB(ofcpdf)(double *y, double *f,double *a, double *b,
		     int *ny, double *x, int *nx, double *bw);
void F77_SUB(remp)(int *n,double *y,double *f, double *a, double *b,
		   int *m, double *Fx, double *x, double *u, int *ntotal);

		  
void F77_SUB(nrlogit)(double *x0, double *betas, double *ps, int *n);
void ppower(double *p0, int *gsize, double *esize, double *alpha,
	    int *ssize, double *pwr);

void orexactl(int *counts, double *alpha, double *out);
void orexactu(int *counts, double *alpha, double *out);


static const R_FortranMethodDef FortEntries[] = {
  {"KSPvalue", (DL_FUNC) & KSPvalue, 1},
  {"emmix", (DL_FUNC) &F77_SUB(emmix), 13},
  {"smoothkde", (DL_FUNC) &F77_SUB(smoothkde),  9},
  {"nrlogit", (DL_FUNC) &F77_SUB(nrlogit),  4},
  {"ppower", (DL_FUNC) & ppower, 6},
  {"orexactl", (DL_FUNC) & orexactl, 3},
  {"orexactu", (DL_FUNC) & orexactu, 3},

  {"hbmise", (DL_FUNC) & hbmise, 5},
  {"ofcpdf", (DL_FUNC) &F77_SUB(ofcpdf), 8},
  {"remp", (DL_FUNC) &F77_SUB(remp), 10},

  {"BDMLE", (DL_FUNC) & BDMLE, 7},
  {"bdregmle", (DL_FUNC) & bdregmle, 7},

  {"tubecv", (DL_FUNC) & tubecv, 2},
  {"lpsmooth", (DL_FUNC) & lpsmooth, 11},
  {"DesignMatrix", (DL_FUNC) & DesignMatrix, 4},

  {"rootGldFmklBisection", (DL_FUNC) & rootGldFmklBisection, 2},

  {"iterfx", (DL_FUNC) &F77_SUB(iterfx),  9},
  {"linbin", (DL_FUNC) &F77_SUB(linbin),  7},
  {"lbtwod", (DL_FUNC) &F77_SUB(lbtwod),  9},
  {"bintwod", (DL_FUNC) &F77_SUB(bintwod),  7},

  {NULL, NULL, 0}
};


void R_init_bda(DllInfo *dll)
{
  //    R_registerRoutines(dll, NULL, NULL, callMethods, NULL);
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
