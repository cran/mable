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
    double *tini, int *progress, int *convergence, double *delta, double *level, int *vb);
void mable_optim_group(int *M, int *N, double *p, double *t, int *n, int *maxit, 
    double *eps, double *lk, double *lr, int *optim, 
    int *progress, int *convergence, double *delta, double *tini, double *bic, 
    double *pval, int *chpts, double *level, int *vb);
void mable_approx(double *u, double *p, int *m, int *n, int *cdf);
void rbeta_mi(int *n, int *m, int *w, double *v);
void mable_mvar(int *M0, int *M, int *n, int *d, int *search, double *phat,  
        int *mhat, double *x, int *maxit,  double *eps, double *lk, int *progress,   
        int *conv, double *D, int *cdf, int *hd);
void mable_mvdf(int *d, int *m, int *km, int *n, double *t, double *p, 
            double *mvdf, int *density);
void optim_gcp(int *M, double *lk, double *lr, int *m, double *pval, int *chpts);

static R_NativePrimitiveArgType em_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, LGLSXP, REALSXP};
static R_NativePrimitiveArgType em_group_t[] = {
    INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, LGLSXP, REALSXP};
static R_NativePrimitiveArgType optim_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, 
    REALSXP, INTSXP, REALSXP, LGLSXP, INTSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType optim_group_t[] = {
    INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,  
    INTSXP, LGLSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType approx_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, LGLSXP};
static R_NativePrimitiveArgType rbeta_t[] = {INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType mvar_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, INTSXP,  
    REALSXP, REALSXP, LGLSXP, INTSXP, REALSXP, LGLSXP, LGLSXP};
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
    double *tau, int *N, double *x0, double *lk, double *lr, double *ddell, double *EPS,  
    int *MAXIT, int *progress, double *pval, int *chpts, double *level, int *conv, int *known_tau);
void mable_aft_m(double *gama, double *p, int *dm, double *x, double *y, double *y2, 
    double *tau, int *N, double *x0, double *ell, double *ddell, double *EPS, int *MAXIT, 
    int *progress, int *conv, double *delta, int *known_tau, int *method);
void mable_aft_gamma(int *M, double *gama, int *dm, double *x, double *y, double *y2, int *N, 
        double *x0, double *lk, double *lr, double *p, double *ddell, double *eps, int *maxit, 
        int *progress, double *pval, int *chpts, double *level, int *conv, double *del, 
        double *tau, int *known_tau);        

static R_NativePrimitiveArgType aft_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,  
    REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, REALSXP, INTSXP, REALSXP, INTSXP, LGLSXP};
static R_NativePrimitiveArgType aft_m_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP, REALSXP, LGLSXP, INTSXP};
static R_NativePrimitiveArgType aft_gamma_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,  
    REALSXP, REALSXP, INTSXP, LGLSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, LGLSXP};

/* GPO Regression */
/*
void mable_po(int *M, double *gama_hat, int *dm, double *p, double *pi0, 
    double *x, double *y, double *y2, int *N, double *x0, double *lk, double *lr,
    double *ddell, double *EPS, int *MAXIT, int *progress, double *level,
    double *pval, int *chpts, int *conv, double *eta, int *eta_known);
void mable_po_gamma(int *M, double *gama, int *dm, double *pi0, double *x, 
    double *y, double *y2, int *N, double *x0, double *lk, double *lr, double *p,
    double *ddell, double *eps, int *maxit, int *progress, double *level, 
    double *pval, int *chpts, int *conv, double *delta, double *eta);
void mable_po_m(double *gama, double *p, int *dm, double *x, double *y, 
    double *y2, int *N, double *x0, double *ell, double *ddell, double *EPS, 
    int *MAXIT, int *progress, int *conv, double *delta, double *eta, int *eta_known);
void weib_gpo(double *theta, int *d, double *x, int *n0, int *n1, double *y,     
    double *y2, double *lk, double *ddell, double *eps, int *maxit, int *prog, 
    int *conv, double *delta, int *eta_known);

static R_NativePrimitiveArgType po_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, 
    INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, 
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, LGLSXP};
static R_NativePrimitiveArgType po_gamma_t[] = {
    INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, 
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, 
    REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType po_m_t[] = {
    REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP, REALSXP, REALSXP, LGLSXP};
static R_NativePrimitiveArgType weib_gpo_t[] = {
    REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,  
    REALSXP, REALSXP, INTSXP, LGLSXP, INTSXP, REALSXP, LGLSXP};
*/

/* Register .C */
static const R_CMethodDef cMethods[] = {
   {"mable_em", (DL_FUNC) &mable_em, 9, em_t},
   {"mable_em_group", (DL_FUNC) &mable_em_group, 10, em_group_t},
   {"mable_optim", (DL_FUNC) &mable_optim, 18, optim_t},
   {"mable_optim_group", (DL_FUNC) &mable_optim_group, 19, optim_group_t},
   {"mable_approx", (DL_FUNC) &mable_approx, 5, approx_t},
   {"rbeta_mi", (DL_FUNC) &rbeta_mi, 4, rbeta_t},
   {"mable_mvar", (DL_FUNC) &mable_mvar, 16, mvar_t},
   {"mable_mvdf", (DL_FUNC) &mable_mvdf, 8, mvdf_t},
   {"optim_gcp", (DL_FUNC) &optim_gcp, 6, gcp_t},
   {"mable_ph", (DL_FUNC) &mable_ph, 20, ph_t},
   {"mable_ph_gamma", (DL_FUNC) &mable_ph_gamma, 21, ph_gamma_t},
   {"mable_ph_m", (DL_FUNC) &mable_ph_m, 15, ph_m_t},
   {"mable_ic", (DL_FUNC) &mable_ic, 18, ic_t},
   {"mable_aft", (DL_FUNC) &mable_aft, 21, aft_t},
   {"mable_aft_m", (DL_FUNC) &mable_aft_m, 18, aft_m_t},
   {"mable_aft_gamma", (DL_FUNC) &mable_aft_gamma, 22, aft_gamma_t},
//   {"mable_po", (DL_FUNC) &mable_po, 22, po_t},
//   {"mable_po_gamma", (DL_FUNC) &mable_po_gamma, 22, po_gamma_t},
//   {"mable_po_m", (DL_FUNC) &mable_po_m, 17, po_m_t},
//   {"weib_gpo", (DL_FUNC) &weib_gpo, 15, weib_gpo_t},
   {NULL, NULL, 0, NULL}
};
/* Deconvolution*/
SEXP mable_decon(SEXP args);
SEXP optim_decon(SEXP args);
SEXP optim_decon_old(SEXP args);

/* Density Ratio */
SEXP C_mable_dr(SEXP args);
SEXP maple_dr(SEXP args);
SEXP C_mable_dr_group(SEXP args);
SEXP maple_dr_group(SEXP args);
SEXP mixtbeta_cdf(SEXP args);
/* Integrate */
//SEXP C_integrate_ibp(SEXP args);

/* Register .External */
static const R_ExternalMethodDef externalMethods[] = {
   {"mable_decon", (DL_FUNC) &mable_decon, 9},
   {"optim_decon", (DL_FUNC) &optim_decon, 11},
   {"C_mable_dr", (DL_FUNC) &C_mable_dr, 19},
   {"C_mable_dr_group", (DL_FUNC) &C_mable_dr_group, 19},
   {"maple_dr", (DL_FUNC) &maple_dr, 19},
   {"maple_dr_group", (DL_FUNC) &maple_dr_group, 19},
   {"mixtbeta_cdf", (DL_FUNC) &mixtbeta_cdf, 8},
//   {"C_integrate_ibp", (DL_FUNC) &C_integrate_ibp, 6},
   {NULL, NULL, 0}
};

void R_init_mable(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, externalMethods);
  R_useDynamicSymbols(info, TRUE);
}
 

