/*
 *  Copyright (C) 2009-2010  B. Wang
 *  Unlimited use and distribution (see LICENCE).
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
//ckernel.f
void wlinbin(double *x, double *w, int *n, double *a,
	     double *b, int *m, int *trun, double *gcounts);
void probin(double *x, double *y, int *n, double *a,
	    double *b, int *m, double *gcounts);
void yldist(double *gcounts, int *m, double *y);
void ofcpdf(double *y, double *f,double *a, double *b,
	    int *ny, double *x, int *nx, double *bw);
void remp(int *n,double *y,double *f, double *a, double *b,
	  int *m, double *Fx, double *x, double *u);

void mlensimp(double *w, double *f,double *a, double *b,
	      int *n, double *theta);


// fitdist.c
void bootsd(int *size,double *x,double *y,double *sx,double *sy,
	    int *m, double *sig, double *rho);
void fitdpro1(double *ll,double *ul, int *n,
	      double *mu,double *s);
void fitdpro2(double *ll,double *ul, int *n1,
	      double *x, int *n2,
	      double *mu,double *s);
// npr.c
void awlprNorm(double *X0, int *n0,
	     double *X, double *Y, double *W, int *size, 
	     double *bw);
void wlprNorm(double *X0, int *n0,
	     double *X, double *Y, double *W, int *size, 
	     double *bw, int *opt);
void lprLap(double *X0, int *n0,
	     double *X, double *Y, double *S, int *size, 
	     double *bw, double *gcv);
void lprHLap(double *X0, int *n0,
	     double *X, double *Y, double *S, int *size, 
	     double *bw, double *gcv);
void nprHLap(double *X0, int *n0,
	     double *X, double *Y, double *S, int *size, 
	     double *bw, double *gcv);

void DkNpReg(double *Z, double *Y, double *S, int *size, 
	     double *bw, double *X, int *nx, 
	     double *loo,double *opt);

void NWReg(double *Z, double *Y, int *size, 
	   double *bw, double *X, int *nx, 
	   double *loo,int *optim, double *opt);

void NormLap2(int *n, double *Rfx,double *ss, double *h1,
	      double *grid,double *ub);
// wnpr.c
void wnpr(double *xgrid, int *ngrid, double *x, double *y, 
	  double *w, int *size, double *bw, double *sx);
void wdekde(double *x, double *w, int *n, double *y, 
	    int *ngrid, double *bw, double *sx);

//wkde.c
void awkde(double *y, double *w, int *n, double *x,double *fx, 
	   double *Fx, int *m, double *pars);

void wkde(double *y, double *w, int *n, double *x,double *fx, 
	  double *Fx, int *m, double *pars);

void wmise(double *x, double *w, int *n,
	   double *h,double *g, int *m);

// cKDE.c
void dKDE(double *x0, int *n0, double *xc, double *h, 
	  double *f, int *n);
void pKDE(double *x0, int *n0, double *xc, double *h, 
	  double *f, int *n);

		    
void rootGldFmklBisection(double *q, double *lambdas);
void KSPvalue(double *x0);
void KSP2x(double *x0, int *n);
void pks2(double *x0, int *m, int *n);
void hbmise(double *x, double *f, double *h, int *n, double *hopt);


// ckernel.c
void orexactl(int *counts, double *alpha, double *out);
void orexactu(int *counts, double *alpha, double *out);
	    
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
void bin2d(double *x, double *y, int *n, double *brk1, int *nb1,
	   double *brk2, int *nb2, double *cnt);


//lognormal.c

void lnormBinChisq(int *size, double *x, double *fn, 
		   double *mu, double *s);

void lnormBinMLE(int *size, double *x, double *fn, 
		 double *mu, double *s);

void lnormBinMLE2(double *a, double *b, double *w, 
		 int *size, double *mu, double *s);

void lnormBinMLE3(double *a, double *b, double *w, 
		 int *size, double *mu, double *s);

void lnormMix1(double *x, double *w, double *xlmt, 
	       int *size, double *pars);

void lnormMixK(double *x, double *w, int *size, double *w0,
	       double *ps, double *mus, double *sigs);

void lnormMixNM(double *x, double *F, int *size,
		int *ncomp, double *mus, double *sigs);

void lnormLSEK(double *x, double *w, int *size, 
	       double *ps, double *mus, double *sigs);

void mclnorm(double *x, double *fn, int *size,
	     double *mu, double *s);

void mclnorm2(double *x, double *fn, double *delta,
	      int *size, double *mu, double *s);

void mleTN(double *x, double *d, double *f, int *size,
	   double *xp, double *qp, double *s);

void mlemixTN(double *x, double *d, double *f, int *size,
	      double *xp, double *qp, double *s, double *pmix, int *k);

void qtlmlnorm(double *x, int *k, double *p, double *mu, double *s);
void slr(double *y, double *x, int *n,
	       double *a, double *b);

// lprsmooth.c
void tubecv(double *kappa, double *level);

void DesignMatrix(double *xv, int *size, double *bw, double *dm);
		    
void lpsmooth(double *xgrid, int *ngrid, double *x, double *y, 
	      int *size, double *bw, int *lscv, 
	      double *range, int *adaptive, double *ellx, 
	      double *kappa);

void llrme(double *xgrid, int *ngrid, double *x, double *y, 
	      int *size, double *bw, int *lscv, 
	      double *range, int *adaptive, double *ellx, 
	      double *kappa);

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


static const R_FortranMethodDef FortEntries[] = {
  {"KSPvalue", (DL_FUNC) & KSPvalue, 1},
  {"KSP2x", (DL_FUNC) & KSP2x, 2},
  {"pks2", (DL_FUNC) & pks2, 3},
  {"bin2d", (DL_FUNC) & bin2d, 8},

  {"BDMLE", (DL_FUNC) & BDMLE, 7},
  {"bdregmle", (DL_FUNC) & bdregmle, 7},

  {"tubecv", (DL_FUNC) & tubecv, 2},
  {"lpsmooth", (DL_FUNC) & lpsmooth, 11},
  {"llrme", (DL_FUNC) & llrme, 11},
  {"DesignMatrix", (DL_FUNC) & DesignMatrix, 4},

  {"rootGldFmklBisection", (DL_FUNC) & rootGldFmklBisection, 2},

  {"em3", (DL_FUNC) & em3, 4},
  {"permtest", (DL_FUNC) & permtest, 9},
  {"permtest2", (DL_FUNC) & permtest2, 6},
  {"permtest3", (DL_FUNC) & permtest3, 5},

  {"lnormBinChisq", (DL_FUNC) & lnormBinChisq, 5},
  {"lnormBinMLE", (DL_FUNC) & lnormBinMLE, 5},
  {"lnormBinMLE2", (DL_FUNC) & lnormBinMLE2, 6},
  {"lnormBinMLE3", (DL_FUNC) & lnormBinMLE3, 6},
  {"lnormMix1", (DL_FUNC) & lnormMix1, 5},
  {"lnormMixK", (DL_FUNC) & lnormMixK, 7},
  {"lnormMixNM", (DL_FUNC) & lnormMixNM, 6},
  {"lnormLSEK", (DL_FUNC) & lnormLSEK, 6},
  {"mclnorm", (DL_FUNC) & mclnorm, 5},
  {"mclnorm2", (DL_FUNC) & mclnorm2, 6},
  {"mleTN", (DL_FUNC) & mleTN, 7},
  {"mlemixTN", (DL_FUNC) & mlemixTN, 9},
  {"qtlmlnorm", (DL_FUNC) & qtlmlnorm, 5},
  {"slr", (DL_FUNC) & slr, 5},

  {"qmPareto", (DL_FUNC) & qmPareto, 5},
  {"mle1Pareto", (DL_FUNC) & mle1Pareto, 5},
  {"mle2Pareto", (DL_FUNC) & mle2Pareto, 5},

  {"orexactl", (DL_FUNC) & orexactl, 3},
  {"orexactu", (DL_FUNC) & orexactu, 3},

  {"hbmise", (DL_FUNC) & hbmise, 5},
  {"ofcpdf", (DL_FUNC) &ofcpdf, 8},
  {"remp", (DL_FUNC) &remp, 9},
  {"mlensimp", (DL_FUNC) &mlensimp, 6},

  {"pKDE", (DL_FUNC) & pKDE, 6},
  {"dKDE", (DL_FUNC) & dKDE, 6},

  {"wmise", (DL_FUNC) & wmise, 6},
  {"awkde", (DL_FUNC) & awkde, 8},
  {"wkde", (DL_FUNC) & wkde, 8},
  {"wlinbin", (DL_FUNC) &wlinbin, 8},
  {"yldist", (DL_FUNC) &yldist,  3},
  {"probin", (DL_FUNC) &probin, 7},

  {"wnpr", (DL_FUNC) & wnpr, 8},
  {"wdekde", (DL_FUNC) & wdekde, 7},
  {"bootsd", (DL_FUNC) & bootsd, 8},
  {"fitdpro1", (DL_FUNC) & fitdpro1, 5},
  {"fitdpro2", (DL_FUNC) & fitdpro2, 7},

  {"DkNpReg", (DL_FUNC) & DkNpReg, 9},
  {"NWReg", (DL_FUNC) & NWReg, 9},
  {"NormLap2", (DL_FUNC) & NormLap2, 6},
  {"nprHLap", (DL_FUNC) & nprHLap, 8},
  {"lprHLap", (DL_FUNC) & lprHLap, 8},
  {"lprLap", (DL_FUNC) & lprLap, 8},
  {"wlprNorm", (DL_FUNC) & wlprNorm, 8},
  {"awlprNorm", (DL_FUNC) & awlprNorm, 7},
  {NULL, NULL, 0}
};


void R_init_bda(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
