#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP cpe1_call(SEXP y, SEXP group, SEXP guess);
SEXP cpe2_call(SEXP y, SEXP group, SEXP guess, SEXP sig);

/* ============================================================
   Result helpers
   ============================================================ */

static inline double sd_welford(const double *x, int n) {
    if (n <= 1) return NA_REAL;
    double mean = 0.0, M2 = 0.0;
    for (int i = 0; i < n; i++) {
        double d = x[i] - mean;
        mean += d / (i + 1);
        M2 += d * (x[i] - mean);
    }
    return sqrt(M2 / (n - 1));
}

/* ============================================================
   OLS: y ~ 1 + group + guess
   ============================================================ */

static void ols3(const double *y, const int *g, const double *x,
                 int n, double coef[3]) {

    double XtX[3][3] = {{0}};
    double XtY[3] = {0};

    for (int i = 0; i < n; i++) {
        double X[3] = {1.0, (double)g[i], x[i]};
        for (int j = 0; j < 3; j++) {
            XtY[j] += X[j] * y[i];
            for (int k = 0; k < 3; k++)
                XtX[j][k] += X[j] * X[k];
        }
    }

    double A[3][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) A[i][j] = XtX[i][j];
        A[i][3] = XtY[i];
    }

    for (int i = 0; i < 3; i++) {
        double piv = A[i][i];
        for (int j = i; j < 4; j++) A[i][j] /= piv;
        for (int k = 0; k < 3; k++) {
            if (k == i) continue;
            double f = A[k][i];
            for (int j = i; j < 4; j++)
                A[k][j] -= f * A[i][j];
        }
    }

    for (int i = 0; i < 3; i++)
        coef[i] = A[i][3];
}

/* ============================================================
   OLS: y ~ 1 + group
   ============================================================ */

static void ols2(const double *y, const int *g,
                 int n, double coef[2]) {

    double s1 = 0, s2 = 0, s4 = 0;
    for (int i = 0; i < n; i++) {
        s1 += 1.0;
        s2 += g[i];
        s4 += g[i] * g[i];
    }

    double det = s1 * s4 - s2 * s2;

    double b0 = 0, b1 = 0;
    for (int i = 0; i < n; i++) {
        b0 += y[i] * s4 - y[i] * g[i] * s2;
        b1 += y[i] * g[i] * s1 - y[i] * s2;
    }

    coef[0] = b0 / det;
    coef[1] = b1 / det;
}

/* ============================================================
   cpe1 (R interface)
   ============================================================ */

SEXP cpe1_call(SEXP yS, SEXP groupS, SEXP guessS) {

    int n = LENGTH(yS);
    const double *y = REAL(yS);
    const double *guess = REAL(guessS);
    const int *group = INTEGER(groupS);

    /* Initial regression */
    double coef3[3];
    ols3(y, group, guess, n, coef3);
    double usham = coef3[2];

    /* Pre-split groups */
    int n1 = 0, n0 = 0;
    for (int i = 0; i < n; i++)
        (group[i] == 1) ? n1++ : n0++;

    int *idx1 = (int*)R_alloc(n1, sizeof(int));
    int *idx0 = (int*)R_alloc(n0, sizeof(int));

    n1 = n0 = 0;
    for (int i = 0; i < n; i++)
        (group[i] == 1) ? (idx1[n1++] = i) : (idx0[n0++] = i);

    int iter = 1000;
    SEXP deltaS = PROTECT(allocVector(REALSXP, iter));
    SEXP resS   = PROTECT(allocVector(REALSXP, iter));

    double *delta = REAL(deltaS);
    double *res   = REAL(resS);

    for (int i = 0; i < iter; i++)
        delta[i] = i * (10.0 * usham) / (iter - 1);

    /* Main loop */
    for (int i = 0; i < iter; i++) {
        double d = delta[i];

        double s1 = 0.0, s0 = 0.0;
        for (int k = 0; k < n1; k++)
            s1 += y[idx1[k]] - d * guess[idx1[k]];
        for (int k = 0; k < n0; k++)
            s0 += y[idx0[k]] - d * guess[idx0[k]];

        double u1 = s1 / n1;
        double u0 = s0 / n0;
        double tau = u1 - u0;

        double mean = 0.0, M2 = 0.0;
        for (int j = 0; j < n; j++) {
            double r = y[j] - d * guess[j] - u0 - tau * group[j];
            double dm = r - mean;
            mean += dm / (j + 1);
            M2 += dm * (r - mean);
        }
        res[i] = sqrt(M2 / (n - 1));
    }

    int imin = 0;
    for (int i = 1; i < iter; i++)
        if (res[i] < res[imin]) imin = i;

    usham = delta[imin];

    /* Final regression */
    double *y2 = (double*)R_alloc(n, sizeof(double));
    for (int i = 0; i < n; i++)
        y2[i] = y[i] - usham * guess[i];

    double coef2[2];
    ols2(y2, group, n, coef2);

    double u2hat = coef2[0];
    double tauhat = coef2[1];

    double sigma = sd_welford(y2, n);

    /* Build return list */
    SEXP out = PROTECT(allocVector(VECSXP, 7));
    SET_VECTOR_ELT(out, 0, deltaS);
    SET_VECTOR_ELT(out, 1, resS);
    SET_VECTOR_ELT(out, 2, ScalarReal(tauhat));
    SET_VECTOR_ELT(out, 3, ScalarReal(u2hat));
    SET_VECTOR_ELT(out, 4, ScalarReal(u2hat + tauhat));
    SET_VECTOR_ELT(out, 5, ScalarReal(usham));
    SET_VECTOR_ELT(out, 6, ScalarReal(sigma));

    SEXP names = PROTECT(allocVector(STRSXP, 7));
    const char *nms[] = {"x","y","tau","ctrl","active","sham","sigma"};
    for (int i = 0; i < 7; i++)
        SET_STRING_ELT(names, i, mkChar(nms[i]));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(4);
    return out;
}

/* ============================================================
   cpe2 (R interface)
   ============================================================ */

SEXP cpe2_call(SEXP yS, SEXP groupS, SEXP guessS, SEXP sigS) {

    int n = LENGTH(yS);
    const double *y = REAL(yS);
    //const double *guess = REAL(guessS);
    //const int *group = INTEGER(groupS);
    double sig = REAL(sigS)[0];

    if (sig <= 0)
        error("sig must be positive");

    int k1 = 30, k2 = 30;

    SEXP gridS = PROTECT(allocVector(REALSXP, k1));
    double *grid = REAL(gridS);

    for (int i = 0; i < k1; i++)
        grid[i] = 0.25 + i * (1.5 - 0.25) / (k1 - 1);

    SEXP outS = PROTECT(allocMatrix(REALSXP, k1, 5));
    double *out = REAL(outS);

    GetRNGstate();

    for (int i = 0; i < k1; i++) {
        for (int j = 0; j < k2; j++) {

            SEXP y1S = PROTECT(allocVector(REALSXP, n));
            double *y1 = REAL(y1S);

            for (int t = 0; t < n; t++)
                y1[t] = y[t] + norm_rand() * grid[i] * sig;

            SEXP r = cpe1_call(y1S, groupS, guessS);

            out[i + k1*0] += REAL(VECTOR_ELT(r,0))[0]; /* tau */
            out[i + k1*1] += REAL(VECTOR_ELT(r,1))[0];
            out[i + k1*2] += REAL(VECTOR_ELT(r,2))[0];
            out[i + k1*3] += REAL(VECTOR_ELT(r,3))[0];
            out[i + k1*4] += REAL(VECTOR_ELT(r,4))[0];

            UNPROTECT(1);
        }

        for (int c = 0; c < 5; c++)
            out[i + k1*c] /= k2;
    }

    PutRNGstate();

    SEXP res = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(res, 0, gridS);
    SET_VECTOR_ELT(res, 1, outS);

    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("x"));
    SET_STRING_ELT(names, 1, mkChar("y"));
    setAttrib(res, R_NamesSymbol, names);

    UNPROTECT(4);
    return res;
}

//dyn.load("cpe_fast_R.so")
//res1 <- .Call("cpe1_call", y, group, guess)
//res2 <- .Call("cpe2_call", y, group, guess, sig)
