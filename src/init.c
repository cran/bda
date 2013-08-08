/*
 *  Copyright (C) 2009-2010  B. Wang
 *  Unlimited use and distribution (see LICENCE).
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Ckernel.c
void weibullmae(double *x,double *w,int *size,
		double *pars, double *y,int *ny);
void expmae(double *x,double *w,int *size,double *y,int *ny);

void tubecv(double *kappa, double *level);

void lpsmooth(double *xgrid, int *ngrid, double *x, double *y, 
	      int *size, double *bw, int *lscv, 
	      double *range, int *adaptive, double *ellx, 
	      double *kappa);
void WeibullMleNMMIN(double *x, double *w, int *nx, double *pars);
void lpshazard(double *xgrid, double *vx, int *ngrid, 
	       double *x, double *y, double *v, 
	       int *size, double *bw, int *lscv, int *ikernel);

// Fkernel.f

void F77_SUB(binning)(double *X, double *F, double *W, double *A, 
		      double *B, int *n, double *xa, double *xb, 
		      int *M, double *gcnts, int *type);
void F77_SUB(yldist)(double *gcounts, int *m, double *y);
void F77_SUB(smoothkde)(double *fx, double x0, int *n, 
			double *x,  double *f, int m,
			double *w, double *h, int *iter);

		    /*  codes below this line need to be double checked
		     */

//  KernelWKDE.c  
void GridBinning(double *x, double *w, int *nx, double *xlo, double *bw,
		 int *ngrid, int *trun, int *linbin, double *gcnts);
void wkdemae(double *x,double *w,int *size,double *y,int *ny);
void RcMleWeibull(double *x,double *w,int *size,double *pars);
void BDMLE(double *f, double *a, double *b, int *nbin,
	   double *pars, int *npar, int *dist);


//  cbootkde.c/fbootkde.f
void hbmise(double *x, double *f, double *h, int *n, double *hopt);
void F77_SUB(ofcpdf)(double *y, double *f,double *a, double *b,
		     int *ny, double *x, int *nx, double *bw);
void F77_SUB(remp)(int *n,double *y,double *f, double *a, double *b,
		   int *m, double *Fx, double *x, double *u, int *ntotal);

// cwkde.c or fwkde.f =======================================
void awkde(double *y, double *w, int *n, double *x,double *fx, 
	   double *Fx, int *m, double *pars);

void wkde(double *y, double *w, int *n, double *x,double *fx, 
	  double *Fx, int *m, double *pars);

void wmise(double *x, double *w, int *n,
	   double *h,double *g, int *m);

void llrGauss(double *x, double *y, int *n, double *x0, int *m, double *bw, int *lscv);

		  /*
void F77_SUB(wlbcounts)(double *x, double *w, int *n, double *a,
		      double *b, int *m, int *trun, double *gcounts);

void F77_SUB(wlinbin)(double *x, double *w, int *n, double *a,
		      double *b, int *m, int *trun, double *gcounts);
		      
void F77_SUB(wbin)(double *x, double *w, int *n, double *a,
			 double *b, int *m, int *trun, double *gcounts);

void F77_SUB(wedf)(double *x, double *w, int *n, 
		   double *xgrid, double *fhat, int *m);
		  */


// cfmm.c or ffmm.f =======================================

void fitmm(double *x, double *counts, double *widths, int *nbin, 
	   int *idist, int *m, double *par, double *llk);
void FitGamma(double *x0,int *m,double *l);
void FitBeta(double *x0,int *m,double *l);
void FitWeibull(double *x0,int *m,double *l);
void F77_SUB(iterfx)(double *fx, double x0, int *n, double *x, double *f, int m,
		     double *w, double *h, double *iter);
void F77_SUB(remlenorm)(double *x, double *f, double *b,int *n, double *theta);
void reemnorm(double *x, double *f, double *b, int *n, int *k, 
	      double *p, double *mu, double *s, double *llk);
void F77_SUB(linbin)(double *x, int *n, double *a,
		     double *b, int *m, double *gcounts);
void F77_SUB(rlbin)(double *x, double *y, int *n,
		    double *a, double *b, int *m, int *trun, double *xcounts,
		    double *ycounts);

void kspvalue(double *x0);

		 
static const R_FortranMethodDef FortEntries[] = {
  {"binning", (DL_FUNC) &F77_SUB(binning),  11},
  {"yldist", (DL_FUNC) &F77_SUB(yldist),  3},
  {"smoothkde", (DL_FUNC) &F77_SUB(smoothkde),  9},
  {"tubecv", (DL_FUNC) & tubecv, 2},
  {"lpsmooth", (DL_FUNC) & lpsmooth, 11},
  {"lpshazard", (DL_FUNC) & lpshazard, 10},
  {"WeibullMleNMMIN", (DL_FUNC) & WeibullMleNMMIN, 4},
  {"weibullmae", (DL_FUNC) & weibullmae, 6},
  {"expmae", (DL_FUNC) & expmae, 5},
  //subroutines below this line need to be double checked.
  {"GridBinning", (DL_FUNC) & GridBinning, 9},
  {"wkdemae", (DL_FUNC) & wkdemae, 5},
  {"RcMleWeibull", (DL_FUNC) & RcMleWeibull, 4},
  {"BDMLE", (DL_FUNC) & BDMLE, 7},
  //codes above this line have been double-checked
  {"kspvalue", (DL_FUNC) & kspvalue, 1},
  {"wmise", (DL_FUNC) & wmise, 6},
  {"awkde", (DL_FUNC) & awkde, 8},
  {"wkde", (DL_FUNC) & wkde, 8},
  /*  {"wlinbin", (DL_FUNC) &F77_SUB(wlinbin),  8},
  {"wlbcounts", (DL_FUNC) &F77_SUB(wlinbin),  8},
  {"wbin", (DL_FUNC) &F77_SUB(wlinbin),  8},
  {"wedf", (DL_FUNC) &F77_SUB(wedf),  6},*/
  {"llrGauss", (DL_FUNC) & llrGauss, 7},
  {"fitmm", (DL_FUNC) & fitmm, 8},
  {"FitGamma", (DL_FUNC) & FitGamma, 3},
  {"FitBeta", (DL_FUNC) & FitBeta, 3},
  {"FitWeibull", (DL_FUNC) & FitWeibull, 3},
  {"reemnorm", (DL_FUNC) &reemnorm, 9},
  {"iterfx", (DL_FUNC) &F77_SUB(iterfx),  9},
  {"remlenorm", (DL_FUNC) &F77_SUB(remlenorm), 5},
  {"linbin", (DL_FUNC) &F77_SUB(linbin),  6},
  {"rlbin",  (DL_FUNC) &F77_SUB(rlbin),   9},
  {"hbmise", (DL_FUNC) & hbmise, 5},
  {"ofcpdf", (DL_FUNC) &F77_SUB(ofcpdf), 8},
  {"remp", (DL_FUNC) &F77_SUB(remp), 10},
  {NULL, NULL, 0}
};


void R_init_bda(DllInfo *dll)
{
  //    R_registerRoutines(dll, NULL, NULL, callMethods, NULL);
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
