/*
 *  Copyright (C) 2009-2010  B. Wang
 *  Unlimited use and distribution (see LICENCE).
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void hbmise(double *x, double *f, double *h, int *n, double *hopt);
void fitmm(double *x, double *counts, double *widths, int *nbin, 
	   int *idist, int *m, double *par, double *llk);
void FitGamma(double *x0,int *m,double *l);
void FitBeta(double *x0,int *m,double *l);
void FitWeibull(double *x0,int *m,double *l);
void kspvalue(double *x0);
void reemnorm(double *x, double *f, double *b, int *n, int *k, 
	      double *p, double *mu, double *s, double *llk);
void F77_SUB(findxyz)(double *y, double *x, int *n, double *y0,  int *m);
void F77_SUB(linbin)(double *x, int *n, double *a,
		     double *b, int *m, double *gcounts);
void F77_SUB(iterfx)(double *fx, double x0, int *n, double *x, double *f, int m,
		     double *w, double *h, double *iter);
void F77_SUB(remlenorm)(double *x, double *f, double *b,int *n, double *theta);
void F77_SUB(ofcpdf)(double *y, double *f,double *a, double *b,
		     int *ny, double *x, int *nx, double *bw);
void F77_SUB(remp)(int *n,double *y,double *f, double *a, double *b,
		   int *m, double *Fx, double *x, double *u);


static const R_FortranMethodDef FortEntries[] = {
  {"findxyz", (DL_FUNC) &F77_SUB(findxyz),  5},
  {"reemnorm", (DL_FUNC) &reemnorm, 9},
  {"linbin", (DL_FUNC) &F77_SUB(linbin),  6},
  {"iterfx", (DL_FUNC) &F77_SUB(iterfx),  9},
  {"remlenorm", (DL_FUNC) &F77_SUB(remlenorm), 5},
  {"kspvalue", (DL_FUNC) & kspvalue, 1},
  {"FitGamma", (DL_FUNC) & FitGamma, 3},
  {"FitBeta", (DL_FUNC) & FitBeta, 3},
  {"FitWeibull", (DL_FUNC) & FitWeibull, 3},
  {"fitmm", (DL_FUNC) & fitmm, 8},
  {"hbmise", (DL_FUNC) & hbmise, 5},
  {"ofcpdf", (DL_FUNC) &F77_SUB(ofcpdf), 8},
  {"remp", (DL_FUNC) &F77_SUB(remp), 9},
  {NULL, NULL, 0}
};


void R_init_bda(DllInfo *dll)
{
  //    R_registerRoutines(dll, NULL, NULL, callMethods, NULL);
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
