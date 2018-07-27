#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
void mable_em(int *m, int *n, double *p, double *x, int *maxit,  double *eps, double *llik);
void mable_em_group(int *m, double *n, int *N, double *p, double *t, int *maxit,  double *eps, double *llik);
void mable_optim(int *ml, int *mu, int *n, double *p, double *x, int *maxit,  double *eps, 
                    double *lk, double *lr, int *optim, double *tini);
void mable_optim_group(int *ml, int *mu, int *N, double *p, double *t, double *n, int *maxit,  double *eps, 
                    double *Llik, double *lr, double *llik, int *optim);
void mable_approx(double *u, double *p, int *m, int *n, int *cdf);

static R_NativePrimitiveArgType em_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType em_group_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType optim_t[] = {
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType optim_group_t[] = {
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};
static R_NativePrimitiveArgType approx_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};


static const R_CMethodDef cMethods[] = {
   {"mable_em", (DL_FUNC) &mable_em, 7, em_t},
   {"mable_em_group", (DL_FUNC) &mable_em_group, 8, em_group_t},
   {"mable_optim", (DL_FUNC) &mable_optim, 11, optim_t},
   {"mable_optim_group", (DL_FUNC) &mable_optim_group, 12, optim_group_t},
   {"mable_approx", (DL_FUNC) &mable_approx, 5, approx_t},
   {NULL, NULL, 0, NULL}
};


void R_init_mable(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
