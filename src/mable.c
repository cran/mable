/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*           C Program for R Package: MABLE               */
/*   Maximum Approximate Bernstein Likelihood Estimation  */
/*                                                        */
/*////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/*                                                               */
/*  Reference:                                                   */
/*   Zhong Guan(2016), Efficient and robust density estimation   */
/*        using Bernstein type polynomials,                      */
/*          J. Nonpar. Stat. 28(2), 250-271                      */
/*   Zhong Guan(2017), Bernstein polynomial model for grouped    */
/*        continuous data, J. Nonpar. Stat. 2017, 29(4), 831-848 */
/*                                                               */
/*///////////////////////////////////////////////////////////////*/
/*  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <string.h>
#include <R.h>
#include <Rconfig.h>
#include <Rmath.h>
#include <stddef.h>
#include <stdio.h>
#include <float.h>
#include <R_ext/Complex.h>
#include <R_ext/Boolean.h>
#include <R_ext/BLAS.h>
#include <R_ext/Error.h>
#include <R_ext/Memory.h>
#include <R_ext/libextern.h>
#include <R_ext/RS.h>
#include <R_ext/Utils.h>
#include <R_ext/Applic.h>
#include <R_ext/Arith.h>
#define TINY 1.0e-20; /*A small number.*/
/*//////////////////////////////////////////*/
/*            Progress Indicators            */
/*//////////////////////////////////////////*/
//#define PRGRSS "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PRGRSS "######################################################################"
#define PBWIDTH 70
void ProgressBar (double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    Rprintf ("\r%3d%% [%.*s%*s]", val, lpad, PRGRSS, rpad, "");
    R_FlushConsole();
}
/*
void clockProgress (double percentage)
{
    char pcl[9]={'|','/','-','\\','|','/','-','\\','\0'};
    int val = (int) (percentage * 100);
    int agl = (int) (val % 8);
    Rprintf ("\r%3d%% [%s]", val, pcl[agl], "");
    R_FlushConsole();
}
*/
/*//////////////////////////////////////////*/
/*    Using Transient storage allocation    */
/* R will reclaim the memory at the end of the call to .C.*/
int *R_ivector(int nl,int nh){
        int *v;
        v=(int *)R_alloc((unsigned) (nh-nl+1), sizeof(int));
        return v-nl;
}

double *R_dvector(int nl,int nh){
        double *v;

        v=(double *)R_alloc((unsigned) (nh-nl+1), sizeof(double));
        return v-nl;
}
/*////////////////////////////////////////////////////////////*/
/* Bernstein base polynomials/(m+1) */
/* dBeta: Returns n x (m+1) matrix, n=length(u) */
/* dbeta(u[i], j+1, m+1-j)/(m+1)= the j-th column, j=0, 1, ..., m */
void dBeta(double *u, int m, int n, double *Bta) {
    int i, j;
    for(i=0; i<n; i++) Bta[i]=(m+1)*R_pow_di(1.0-u[i], m);
    for(i=0; i<n; i++) {
        if(u[i]<1){
            j=0;
            while(j<m){
                Bta[i+n*(j+1)]=(m-j)*(u[i]/(1.0-u[i]))*Bta[i+n*j]/(double)(j+1.0);
//                Rprintf("  Bta = %g\n",  Bta[i+n*(j+1)]);
                j++;
            }
        }
        else{
            for(j=1;j<m;j++) Bta[i+n*j]=0.0;
            Bta[i+n*m]=(m+1);
        }
    }
}
/*////////////////////////////////////////////////////////////*/
/*                    CDF of beta(i+1, m-i+1)                 */
/*        pBeta: Returns n x (m+1) matrix, n=length(u)        */
/*   pbeta(u[j], i+1, m+1-i)= the j-th column, i=0, 1, ..., m */
/*////////////////////////////////////////////////////////////*/
void pBeta(double *u, int m, int n, double *cdfBta) {
    int i, j;
    for(j=0; j<n;j++){
        for(i=0; i<=m; i++){
            cdfBta[j+n*i]=pbeta(u[j], i+1, m-i+1, TRUE, FALSE);
//            Rprintf("  Bta = %g\n",  cdfBta[j+n*i]);
        }
    }
}
/*////////////////////////////////////////////////////////////*/
/* CDF of beta(i+1, m-i+1)                                    */
/* cdf_Beta: Returns (m+1) x n matrix, N=length(u)            */
/* pbeta(u[j+1], i+1, m+1-i)-pbeta(u[j], i+1, m+1-i)          */
/*                          = the j-th column, i=0, 1, ..., m */
/*////////////////////////////////////////////////////////////*/
void cdf_Beta(double *u, int m, int N, double *cdfBta) {
    int i, j;
    for(j=0; j<N;j++){
        for(i=0; i<=m; i++){
            cdfBta[i+(m+1)*j]=pbeta(u[j+1], i+1, m-i+1, TRUE, FALSE)-pbeta(u[j], i+1, m-i+1, TRUE, FALSE);
//            Rprintf("  Bta = %g\n",  cdfBta[i+(m+1)*j]);
        }
    }
}
/*////////////////////////////////////////////////////////////*/
/*   Log-Bernstein-Likelihood                                 */
/*  p=(p0, p1, ..., pm),                                      */
/*  x: n-vector, sample from distribution F with support [0,1]*/
/*  Bta: (m+1) x n matrix of Bernstein base polynomials       */
/*////////////////////////////////////////////////////////////*/
double loglik_bern(double *p, double *Bta, int m, int n){
    int i,j;
    double llik, fx;
    llik = 0.0;
    for(i=0; i<n; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*Bta[i+n*j];
//            Rprintf("  p = %g\n", p[j]);
        }
        llik += log(fx);
    }
//    Rprintf("  lik = %g\n", llik);
//    return(llik+n*log(m+1.));
    return(llik);
}
/*////////////////////////////////////////////////////////////////*/
/*   Log-Bernstein-Likelihood for Grouped Data                    */
/*   p=(p0, p1, ..., pm),                                         */
/*   u: N-vector, grouped data from F with support [0,1]          */
/*   cdfBta: (m+1) x N matrix of Beta[t(j+1)]-Beta[t(j)]          */
/*////////////////////////////////////////////////////////////////*/
double loglik_bern_group(double *p, double *cdfBta, double *n, int m, int N){
    int i,j;
    double llik, fx;
    llik = 0.0;
    for(i=0; i<N; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*cdfBta[j+(m+1)*i];
//            Rprintf("  p = %g\n", p[j]);
        }
        llik += n[i]*log(fx);
    }
//    Rprintf("  lik = %g\n", llik);
    return(llik);
}
/*//////////////////////////////////////////////////*/
/*                 EM Method                        */
/*//////////////////////////////////////////////////*/
/* EM Method for mixture of beta(i+1, m+1-i), i=0,...,m. */
void em_beta_mix(double *p, double *Bta, int m, int n, int maxit, double eps, double *llik){
    int i, j, it;
    double del, llik_nu, *pBta,  *fp, *pnu;
    pBta = Calloc((m+1)*n, double);
    fp = Calloc(n, double);
    pnu = Calloc(m+1, double);
    llik[0] = loglik_bern(p, Bta, m, n);
    del = 10.0;
    it = 1;
    while(del>eps && it<maxit){
        for(j=0; j<n; j++){
            fp[j] = 0.0;
            for(i=0; i<=m; i++) {
                pBta[j+n*i] = p[i]*Bta[j+n*i];
                fp[j] += pBta[j+n*i];
            }
        }
        for(i=0; i<=m; i++){
            pnu[i] = 0.0;
            for(j=0; j<n; j++) pnu[i] += pBta[j+n*i]/fp[j];
            pnu[i] /= (double)n;
        }
        llik_nu = loglik_bern(pnu, Bta, m, n);
        del = fabs(llik[0]-llik_nu);
        it += 1;
        for(i=0; i<=m; i++) p[i] = pnu[i];
        llik[0] = llik_nu;
    }
    Free(pBta);
    Free(fp);
    Free(pnu);
}
/* end function em_beta_mix */
void mable_em(int *m, int *n, double *p, double *x, int *maxit,  double *eps, double *llik){
    double *Bta;
//    Rprintf("\n Program 'mble_em' is runing. This may take several minutes.\n\n\n");
    Bta = Calloc((*m+1)*(*n), double);
    dBeta(x, *m, *n, Bta);
    em_beta_mix(p, Bta, *m, *n,  *maxit, *eps, llik);
    Free(Bta);
//    Rprintf("\n 'mble_em' Finished. \n\n\n");
}
/* Grouped data EM Method for mixture of beta(i+1, m+1-i), i=0,...,m. */
/* n = (n1,...,n_N), t=(t0, t1,...,t_N), */
void em_beta_mix_group(double *p,  double *cdfBta, int N, int m, double *n, int maxit, double eps, double *llik){
    int i, j, it, nn;
    double del, llik_nu, *pBta,  *pnu, *fp;
    pBta = Calloc((m+1)*N, double);
    fp = Calloc(N, double);
    pnu = Calloc(m+1, double);
    llik[0] = loglik_bern_group(p, cdfBta, n, m, N);
//Rprintf("llik=%f\n",llik[0]);
    nn = 0;
    for(j=0; j<N;j++) nn+=n[j];
    del = 10.0;
    it = 1;
    while(del>eps && it<maxit){
//Rprintf("    m=%d  Iteration=%d\n",m,it);
        for(j=0; j<N;j++){
            fp[j] = 0.0;
            for(i=0;i<=m;i++){
                pBta[i+(m+1)*j] = p[i]*cdfBta[i+(m+1)*j];
                fp[j] += pBta[i+(m+1)*j];
            }
        }
        for(i=0; i<=m; i++){
            pnu[i] = 0.0;
            for(j=0; j<N; j++) pnu[i] += n[j]*pBta[i+(m+1)*j]/fp[j];
            pnu[i] /= (double)nn;
        }
        llik_nu = loglik_bern_group(pnu, cdfBta, n, m, N);
        del = fabs(llik[0]-llik_nu);
        it += 1;
        for(i=0; i<=m; i++) p[i] = pnu[i];
        llik[0] = llik_nu;
    }
//Rprintf("llik=%f\n",llik[0]);
    Free(pBta);
    Free(fp);
    Free(pnu);
}
/* end function em_beta_mix_group */

void mable_em_group(int *m, double *n, int *N, double *p, double *t, int *maxit,  double *eps, double *llik){
    double *cdfBta;
    cdfBta = Calloc((*m+1)*(*N), double);
    cdf_Beta(t, *m, *N, cdfBta);
//    Rprintf("\n Program mble_em_group is runing. This may take several minutes.\n\n\n");
    em_beta_mix_group(p, cdfBta, *N, *m, n,  *maxit, *eps, llik);
    Free(cdfBta);
//    Rprintf("\n mable_em_group Finished. \n\n\n");
}
////////////////////////////////////////////////////////
// choosing optimal degree m by change-point method
//  using exponential chpt models
////////////////////////////////////////////////////////
void mable_optim(int *ml, int *mu, int *n, double *p, double *x, int *maxit,  double *eps,
                    double *lk, double *lr, int *optim, double *tini){
    int i, j, l, m,  k;
    double *pt, *Bta, maxLR=0, *u, *beta_mi, *llik, pct, tmp;
    llik = Calloc(1, double);
    Bta = Calloc((*mu+1)*(*n), double);
    u = Calloc((*mu+1), double);
    beta_mi = Calloc((*mu+1)*(*mu+1), double);
//    Rprintf("\n Program 'mble_optim' is runing. This may take several minutes.\n");
    k=*mu-*ml;
    tmp = (double)(k*(k-1));
    m=*mu;
    pt = Calloc(m+1, double);
    for(j=0; j<=m; j++) p[j]=1.0/(double) (m+1);
    dBeta(x, m, *n, Bta);
    em_beta_mix(p, Bta, m, *n,  *maxit, *eps, llik);
    lk[k] = *llik;
    lr[k-1] = 0.0;
    m=*ml;
    for(j=0; j<=m; j++) p[j]=1.0/(double) (m+1);
    dBeta(x, m, *n, Bta);
    em_beta_mix(p, Bta, m, *n,  *maxit, *eps, llik);
    lk[0] = *llik;
    for(i=1; i<k;i++){
        // updating p
        p[m+1] = (m+1)*p[m]/(double)(m+2);
        for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
        p[0] = (m+1.)*p[0]/(double)(m+2.);
        // updating Bta
        for(j=0;j<*n;j++)
            Bta[j+*n*(m+1)]=(m+2.)*x[j]*Bta[j+*n*m]/(m+1.);
        for(j=0;j<=m;j++){
            for(l=0;l<*n;l++) Bta[l+*n*j]=(m+2.)*(1.0-x[l])*Bta[l+*n*j]/(double)(m+1-j);
        }
        m = *ml+i;
        //make sure initial p is in the interior of the simplex
        for(j=0; j<=m; j++) p[j] =(1.0-*tini)*p[j]+ *tini/(double)(m+1.0);
        em_beta_mix(p, Bta, m, *n,  *maxit, *eps, llik);
        lk[i] = *llik;
        //   lr: exponential LR of change-point
        lr[i-1] = k*log((lk[k]-lk[0])/(double)k)-i*log((lk[i]-lk[0])/(double)i)-(k-i)*log((lk[k]-lk[i])/(double)(k-i));
        if(lr[i-1]>maxLR) {
           maxLR = lr[i-1];
           *optim = m;
           for(j=0; j<=m; j++)  pt[j] = p[j];
        }
        pct = i*(i+1)/tmp;
        ProgressBar(pct);
//        clockProgress(pct);
    }
    m=*optim;
    for(j=0;j<=m;j++) p[j]=pt[j];
    Free(pt);  Free(Bta); Free(u); Free(beta_mi); Free(llik);
    Rprintf("\n");
//    Rprintf("\n Program  'mble_optim' Finished. \n\n\n");
}
/*//////////////////////////////////////////////////////////*/
void mable_optim_group(int *ml, int *mu, int *N, double *p, double *t, double *n, int *maxit,  double *eps,
                    double *Llik, double *lr, double *llik, int *optim){
    int i, j, m,  k, ii;
    double *pt, *cdfBta, tmp1, pct, maxlik, tmp2;
    cdfBta = Calloc((*mu+1)*(*N), double);
//    Rprintf("\n Program 'mble_optim_group' is runing. This may take several minutes.\n");
    k=*mu-*ml;
    tmp2 = (double)(k*(k+1));
    m=*mu;
    pt = Calloc(m+1, double);
    for(j=0; j<=m; j++) p[j]=1.0/(double) (m+1);
    cdf_Beta(t, m, *N, cdfBta);
    em_beta_mix_group(p, cdfBta, *N, m, n,  *maxit, *eps, llik);
    Llik[k] = *llik;
    for(i=0; i<=1;i++){
        m=i+*ml;
        for(j=0; j<=m; j++) p[j]=1.0/(double) (m+1);
        cdf_Beta(t, m, *N, cdfBta);
        em_beta_mix_group(p, cdfBta, *N,  m, n,  *maxit, *eps, llik);
        Llik[i] = *llik;
    }
    ii=1;
    maxlik=Llik[ii];
    *optim=m;
    tmp1 = .0;
    for(j=0;j<=m;j++) {
        pt[j]=p[j];
    }
    while(m<*mu){
        lr[ii-1]=ii*log((Llik[ii]-Llik[0])/(double)ii)+(k-ii)*log((Llik[k]-Llik[ii])/(double)(k-ii))-k*log((Llik[k]-Llik[0])/(double)k);
        if(lr[ii-1]<tmp1){
            *optim=m;
            for(j=0;j<=m;j++) pt[j]=p[j];
            tmp1 = lr[ii-1];
            maxlik = Llik[ii];
        }
        for(j=0;j<=m+1;j++) p[j]=1.0/(double) (m+1);
        // finish calculating initial for p
        m++;
        ii++;
        cdf_Beta(t, m, *N, cdfBta);
        pct = ii*(ii+1)/tmp2;
        ProgressBar(pct);
        em_beta_mix_group(p, cdfBta, *N, m, n, *maxit, *eps, llik);
        Llik[ii] = *llik;
    }
    *llik = maxlik;
    m=*optim+1;
    for(j=0;j<=m;j++) p[j]=pt[j];
    Free(pt);
    Rprintf("\n");
//    Rprintf("\n Program  'mble_optim_group' Finished. \n\n\n");
}
/*//////////////////////////////////////////////////////////*/
/*//////////////////////////////////////////*/
/*   Bernstein Polynomial Approximation     */
/*  Returns n-dim row vector, n=length(u)   */
/*//////////////////////////////////////////*/
void mable_approx(double *u, double *p, int *m, int *n, int *cdf){
    int i, j;
    double *Bta, tmp;
    Bta = Calloc(*n*(*m+1), double);
    if(*cdf==0) dBeta(u, *m, *n, Bta);
    if(*cdf==1) pBeta(u, *m, *n, Bta);
    for(j=0;j<*n;j++) {
       tmp=.0;
       for(i=0;i<=*m;i++){
          tmp += Bta[j+*n*i]*p[i];
       }
       u[j]=tmp;
    }
    Free(Bta);
}
