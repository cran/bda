/*
 *  Copyright (C) 2009-2010  B. Wang
 *  Unlimited use and distribution (see LICENCE).
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// ckernel.c
void orexactl(int *counts, double *alpha, double *out);
void orexactu(int *counts, double *alpha, double *out);

//Fkernel.f
void F77_SUB(linbin)(double *x, int *n, double *a, double *b,
		     int *m, int *trun, double *gcounts);

void F77_SUB(lbtwod)(double *x, int *n, double *a1,
		     double *a2, double *b1, double *b2, int *m1,
		     int *m2, double *gcounts);

void F77_SUB(bintwod)(double *x, int *n, double *g1,
		      double *g2, int *m1,
		      int *m2, double *gcounts);

void F77_SUB(iterfx)(double *fx, double *x0, int *n, double *x,
		     double *f, int *m,
		     double *w, double *h, int *iter);

void F77_SUB(smoothkde)(double *fx, double *x0, int *n, 
			double *x,  double *f, int *m,
			double *w, double *h, int *iter);
		    
void rootGldFmklBisection(double *q, double *lambdas);
void KSPvalue(double *x0);
void KSP2x(double *x0, int *n);
void pks2(double *x0, int *m, int *n);

// permtest.c
void permtest(double *x, int *nx,  double *y, int *ny,
	      double *a, double *b,
	      double *D, double *pv, int *iter);
void permtest2(double *D, int *M, int *nx, int *ny,
	       int *F, int *iter);
void permtest3(double *xy, int *nx, int *ny, double *pv, int *iter);

// em.c
void em3(int *size, double *x, double *pars, double *tol);

//lognormal.c

void lnormBinChisq(int *size, double *x, double *fn, 
		   double *mu, double *s);

void lnormBinMLE(int *size, double *x, double *fn, 
		 double *mu, double *s);

void mclnorm(double *x, double *fn, int *size,
	     double *mu, double *s);

void mclnorm2(double *x, double *fn, double *delta,
	      int *size, double *mu, double *s);

void mleTN(double *x, double *d, double *f, int *size,
	   double *xp, double *qp, double *s);

void mlemixTN(double *x, double *d, double *f, int *size,
	      double *xp, double *qp, double *s, double *pmix, int *k);

void qtlmlnorm(double *x, int *k, double *p, double *mu, double *s);

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


//  KernelWKDE.c  
void BDMLE(double *f, double *a, double *b, int *nbin,
	   double *pars, int *npar, int *dist);
void bdregmle(double *F, double *x,double *freq, 
	      int *nu, int *n, int *dist, double *pars);

//  Pareto.c  
void qmPareto(double *p, double *q, int *n,
	      double *xm, double *alpha);
void mle1Pareto(double *cnts, double *b, int *nclass,
		double *xm, double *alpha);
void mle2Pareto(double *cnts, double *b, int *nclass,
		double *xm, double *alpha);

//  cbootkde.c/fbootkde.f
void hbmise(double *x, double *f, double *h, int *n, double *hopt);
void F77_SUB(ofcpdf)(double *y, double *f,double *a, double *b,
		     int *ny, double *x, int *nx, double *bw);
void F77_SUB(remp)(int *n,double *y,double *f, double *a, double *b,
		   int *m, double *Fx, double *x, double *u, int *ntotal);

static const R_FortranMethodDef FortEntries[] = {
  {"KSPvalue", (DL_FUNC) & KSPvalue, 1},
  {"KSP2x", (DL_FUNC) & KSP2x, 2},
  {"pks2", (DL_FUNC) & pks2, 3},

  {"emmix", (DL_FUNC) &F77_SUB(emmix), 13},
  {"smoothkde", (DL_FUNC) &F77_SUB(smoothkde),  9},

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

  {"em3", (DL_FUNC) & em3, 4},
  {"permtest", (DL_FUNC) & permtest, 9},
  {"permtest2", (DL_FUNC) & permtest2, 6},
  {"permtest3", (DL_FUNC) & permtest3, 5},

  {"lnormBinChisq", (DL_FUNC) & lnormBinChisq, 5},
  {"lnormBinMLE", (DL_FUNC) & lnormBinMLE, 5},
  {"mclnorm", (DL_FUNC) & mclnorm, 5},
  {"mclnorm2", (DL_FUNC) & mclnorm2, 6},
  {"mleTN", (DL_FUNC) & mleTN, 7},
  {"mlemixTN", (DL_FUNC) & mlemixTN, 9},
  {"qtlmlnorm", (DL_FUNC) & qtlmlnorm, 5},

  {"qmPareto", (DL_FUNC) & qmPareto, 5},
  {"mle1Pareto", (DL_FUNC) & mle1Pareto, 5},
  {"mle2Pareto", (DL_FUNC) & mle2Pareto, 5},

  {"orexactl", (DL_FUNC) & orexactl, 3},
  {"orexactu", (DL_FUNC) & orexactu, 3},

  {NULL, NULL, 0}
};


void R_init_bda(DllInfo *dll)
{
  //    R_registerRoutines(dll, NULL, NULL, callMethods, NULL);
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
