#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
/* One-sample problem */
void mable_em(int *m, int *n, double *p, double *x, int *maxit,  double *eps, 
    double *llik, int *convergence, double *delta);
void mable_em_group(int *m, int *n, int *N, double *p, double *t, int *maxit, 
    double *eps, double *llik, int *convergence, double *delta);
void mable_optim(int *M, int *n, double *p, double *x, int *maxit, double *eps,  
    double *lk, double *lr, int *optim, double *pval, double *bic, int *chpts, 
    double *tini, int *progress, int *convergence, double *delta, double *level);
void mable_optim_group(int *M, int *N, double *p, double *t, int *n, int *maxit, 
    double *eps, double *lk, double *lr, int *optim, 
    int *progress, int *convergence, double *delta, double *tini, double *bic, 
    double *pval, int *chpts, double *level);
void mable_approx(double *u, double *p, int *m, int *n, int *cdf);
void mable_mvar(int *m, int *n, int *d, int *km, double *p, double *x, 
        int *maxit,  double *eps, double *llik, int *progress, int *conv);
void mable_mvdf(int *d, int *m, int *km, int *n, double *t, double *p, 
            double *mvdf, int *density);
void optim_gcp(int *M, double *lk, double *lr, int *m, double *pval, int *chpts);

static R_NativePrimitiveArgType em_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, LGLSXP, REALSXP};
static R_NativePrimitiveArgType em_group_t[] = {
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, LGLSXP, REALSXP};
static R_NativePrimitiveArgType optim_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, 
    REALSXP, INTSXP, REALSXP, LGLSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType optim_group_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,  
    INTSXP, LGLSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType approx_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, LGLSXP};
static R_NativePrimitiveArgType mvar_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, 
    LGLSXP, INTSXP};
static R_NativePrimitiveArgType mvdf_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, LGLSXP};
static R_NativePrimitiveArgType gcp_t[] = {
    INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP};

/* PH Regression */
void mable_ph(int *M, double *gama_hat, int *dm, double *p, double *pi0, 
    double *x, double *y, double *y2, int *N, double *x0, double *lk, double *lr,
    double *ddell, double *EPS, int *MAXIT, int *progress, double *level,
    double *pval, int *chpts, int *conv);
void mable_ph_gamma(int *M, double *gama, int *dm, double *pi0, double *x, 
    double *y, double *y2, int *N, double *x0, double *lk, double *lr, double *p,
    double *ddell, double *eps, int *maxit, int *progress, double *level, 
    double *pval, int *chpts, int *conv, double *delta);
void mable_ph_m(double *gama, double *p, int *dm, double *x, double *y, 
    double *y2, int *N, double *x0, double *ell, double *ddell, double *EPS, 
    int *MAXIT, int *progress, int *conv, double *delta);
void mable_ic(int *M, double *pi0, double *y, double *y2, int *N, double *lk,
    double *lr, double *p, double *eps, int *maxit, int *progress,
    double *pval, double *bic, int *chpts, int *optim, double *level, 
    int *conv, double *delta);

static R_NativePrimitiveArgType ph_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, 
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, 
    REALSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType ph_gamma_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, 
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, 
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType ph_m_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType ic_t[] = {
    INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, 
    INTSXP, LGLSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};

/* AFT Regression */
void mable_aft(int *M, double *gama, int *dm, double *p, double *x, double *y, double *y2,   
    int *N, double *x0, double *lk, double *lr, double *ddell, double *EPS, int *MAXIT, 
    int *progress, double *pval, int *chpts, double *level, int *conv);
void mable_aft_m(double *gama, double *p, int *dm, double *x, double *y, double *y2, 
    int *N, double *x0, double *ell, double *ddell, double *EPS, int *MAXIT, 
    int *progress, int *conv, double *delta);
void mable_aft_gamma(int *M, double *gama, int *dm, double *x, double *y, double *y2, int *N, 
        double *x0, double *lk, double *lr, double *p, double *ddell, double *eps, int *maxit, 
        int *progress, double *pval, int *chpts, double *level, int *conv, double *del);        

static R_NativePrimitiveArgType aft_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,  
    REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, REALSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType aft_m_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType aft_gamma_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, 
    REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP};

/* Register .C */
static const R_CMethodDef cMethods[] = {
   {"mable_em", (DL_FUNC) &mable_em, 9, em_t},
   {"mable_em_group", (DL_FUNC) &mable_em_group, 10, em_group_t},
   {"mable_optim", (DL_FUNC) &mable_optim, 17, optim_t},
   {"mable_optim_group", (DL_FUNC) &mable_optim_group, 18, optim_group_t},
   {"mable_approx", (DL_FUNC) &mable_approx, 5, approx_t},
   {"mable_mvar", (DL_FUNC) &mable_mvar, 11, mvar_t},
   {"mable_mvdf", (DL_FUNC) &mable_mvdf, 8, mvdf_t},
   {"optim_gcp", (DL_FUNC) &optim_gcp, 6, gcp_t},
   {"mable_ph", (DL_FUNC) &mable_ph, 20, ph_t},
   {"mable_ph_gamma", (DL_FUNC) &mable_ph_gamma, 21, ph_gamma_t},
   {"mable_ph_m", (DL_FUNC) &mable_ph_m, 15, ph_m_t},
   {"mable_ic", (DL_FUNC) &mable_ic, 18, ic_t},
   {"mable_aft", (DL_FUNC) &mable_aft, 19, aft_t},
   {"mable_aft_m", (DL_FUNC) &mable_aft_m, 15, aft_m_t},
   {"mable_aft_gamma", (DL_FUNC) &mable_aft_gamma, 20, aft_gamma_t},
   {NULL, NULL, 0, NULL}
};

SEXP mable_decon(SEXP args);
//static R_NativePrimitiveArgType decon_t[] = {SEXP};
/* Register .External */
static const R_ExternalMethodDef externalMethods[] = {
   {"mable_decon", (DL_FUNC) &mable_decon, 9},
   {NULL, NULL, 0}
};

void R_init_mable(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, externalMethods);
  R_useDynamicSymbols(info, TRUE);
}
 

