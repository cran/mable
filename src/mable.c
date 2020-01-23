/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*      C Utility Program for R Package: MABLE            */
/*   Maximum Approximate Bernstein Likelihood Estimation  */
/*                                                        */
/*////////////////////////////////////////////////////////*/
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
#include <Rinternals.h>
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
void ProgressBar (double percentage, char *txt) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    Rprintf ("\r%s%3d%% [%.*s%*s]",txt, val, lpad, PRGRSS, rpad, "");
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
/*//////////////////////////////////////////////////*/
/*               LU Decomposition                   */
/*//////////////////////////////////////////////////*/

/* C code modified from Numerical Recipe in C, store matrix a as vector A*/
/*  a[i][j] = A[i+j*n], i,j=0,...,n-1. */
void ludcmp(double *A, int n, int *indx, double *d){
    int i,imax,j,k;
    double big,dum,sum,temp;
    double *vv; /*vv stores the implicit scaling of each row.*/
    vv=Calloc(n, double);
    imax=0;
    *d=1.0; /*No row interchanges yet.*/
    for (i=0;i<n;i++) { /*Loop over rows to get the implicit scaling information.*/
        big=0.0;
        for (j=0;j<n;j++)
            if ((temp=fabs(A[i+j*n])) > big) big=temp;
        if (big == 0.0) error("Singular matrix in routine ludcmp\n");
        /*No nonzero largest element.*/
        vv[i]=1.0/big; /*Save the scaling.*/
    }
    for (j=0;j<n;j++) {
        /* This is the loop over columns of Crout  method. */
        for (i=0;i<j;i++) {
        /*This is equation (2.3.12) except for i = j.*/
            sum=A[i+j*n];
            for (k=0;k<i;k++) sum -= A[i+k*n]*A[k+j*n];
            A[i+j*n]=sum;
        }
        big=0.0; /*Initialize for the search for largest pivot element.*/
        for (i=j;i<n;i++){
        /*This is i = j of equation (2.3.12) and i = j+1...N of equation (2.3.13).*/
            sum=A[i+j*n];
            for(k=0;k<j;k++)
                sum -= A[i+k*n]*A[k+j*n];
            A[i+j*n]=sum;
            if((dum=vv[i]*fabs(sum)) >= big){
            /*Is the figure of merit for the pivot better than the best so far?*/
                big=dum;
                imax=i;
            }
        }
        if (j != imax) { /*Do we need to interchange rows?*/
            for(k=0;k<n;k++) { /*Yes, do so...*/
                dum= A[imax+k*n];
                A[imax+k*n]=A[j+k*n];
                A[j+k*n] =dum;
            }
            *d = -(*d); /*...and change the parity of d.*/
            vv[imax]=vv[j];/* Also interchange the scale factor.*/
        }
        indx[j]=imax;
        if (A[j+j*n] == 0.0) A[j+j*n]=TINY;
        /* If the pivot element is zero the matrix is singular (at least to the precision of the*/
        /* algorithm). For some applications on singular matrices, it is desirable to substitute*/
        /* TINY for zero.*/
        if (j != n-1) { /*Now, finally, divide by the pivot element.*/
            dum=1.0/A[j+j*n];
            for (i=j+1;i<n;i++) A[i+j*n] *= dum;
        }
    } /*Go back for the next column in the reduction.*/
    Free(vv);
 }
/*//////////////////////////////////////////////////*/
/*          Forward and backsubstitution            */
/*//////////////////////////////////////////////////*/

void lubksb(double *A, int n, int *indx, double b[]){
    int i,ii=0,ip,j;
    double sum;
    for (i=0;i<n;i++) {  /*When ii is set to a positive value, it will become the*/
        ip=indx[i];      /*index of the first nonvanishing element of b. We now*/
        sum=b[ip];       /*do the forward substitution, equation (2.3.6). The*/
        b[ip]=b[i];      /*only new wrinkle is to unscramble the permutation*/
        if(ii)           /*as we go.*/
            for (j=ii-1;j<=i-1;j++) sum -= A[i+j*n]*b[j];
        else if(sum) ii=i+1; /*A nonzero element was encountered, so from now on we*/
        b[i]=sum;           /*will have to  do the sums in the loop above.*/
    }
    for (i=n-1;i>=0;i--) {  /*Now we do the backsubstitution, equation (2.3.7).*/
        sum=b[i];
        for (j=i+1;j<n;j++) sum -= A[i+j*n]*b[j];
        b[i]=sum/A[i+i*n]; /*Store a component of the solution vector X.*/
    } /*All done!*/
}


/*//////////////////////////////////////////////////*/
/*               Matrix Inverse                     */
/*//////////////////////////////////////////////////*/

void minverse(double *A, int N){
    double *tmp, d, *col;
    int i,j,*indx;
    d=0.0;
    indx = Calloc(N, int); 
    col = Calloc(N, double);
    tmp = Calloc(N*N, double);
    ludcmp(A,N,indx,&d);
    /* Decompose the matrix just once.*/
    for(j=0;j<N;j++) {
        /*Find inverse by columns.*/
        for(i=0;i<N;i++) col[i]=0.0;
        col[j]=1.0;
        lubksb(A,N,indx,col);
        for(i=0;i<N;i++) tmp[i+j*N]=col[i];
    }
    for(j=0;j<N;j++)
        for(i=0;i<N;i++) A[i+j*N]=tmp[i+j*N];
        /* return inverse of A using A */
    Free(tmp); Free(col); Free(indx);
}
/* Print matrix */
void Print_Matrix(double *m, int nr, int nc, char *mname){
    int i,j;
    Rprintf("%s:\n", mname);
    for(i=1; i<=nr; i++) {
        for(j=1;j<=nc;j++) Rprintf("  %s[%d][%d] = %f, \t", mname, i,j, m[i+(j-1)*nr]);
        Rprintf("\n");
    }
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
/*   pbeta(u[i], j+1, m+1-j)= the j-th column, j=0, 1, ..., m */
/*////////////////////////////////////////////////////////////*/
void pBeta(double *u, int m, int n, double *pBta) {
    int i, j;
    for(i=0; i<n;i++){
        for(j=0; j<=m; j++){
            pBta[i+n*j]=pbeta(u[i], j+1, m-j+1, TRUE, FALSE);
//            Rprintf("  Bta = %g\n",  pBta[i+n*j]);
        }
    }
}
/*////////////////////////////////////////////////////////////*/
/*          Class Probability of beta(j+1, m-j+1)             */
/*    cpBeta: Returns N x (m+1) matrix, N=length(u)-1        */
/*    pbeta(u[i+1], j+1, m+1-j)-pbeta(u[i], j+1, m+1-j)       */
/*                          = the j-th column, j=0, 1, ..., m */
/*////////////////////////////////////////////////////////////*/
void cpBeta(double *u, int m, int N, double *dBta) {
    int i, j;
    for(i=0; i<N;i++){
        for(j=0; j<=m; j++){
            dBta[i+N*j]=pbeta(u[i+1], j+1, m-j+1, TRUE, FALSE)-pbeta(u[i], j+1, m-j+1, TRUE, FALSE);
        }
    }
}
/*///////////////////////////////////////////////////////////////////*/
/* Calculating 1-Beta(y): n x (m+1), [beta(y), 1-Beta(y)]: n x (m+2) */
/* Bdata: Returns n x (m+1) matrices Bta and bBta, n=length(y)       */
/*   Bta:  dbeta(y[i], j+1, m+1-j), i=0,...,n-1                      */
/*  bBta:  1-pbeta(y[i], j+1, m+1-j), i=0,...,n-1                    */
/*                         j=0, 1, ..., m                            */
/*///////////////////////////////////////////////////////////////////*/
void Bdata(double *y, int m, int n0, int n1, double *Bta) {
    int i, j, mp1=m+1, n=n0+n1;
    for(i=0; i<n0; i++){
        for(j=0;j<=m;j++) Bta[i+n*j]=dbeta(y[i], j+1, m-j+1, FALSE);
        Bta[i+n*mp1]=0.0;}
    for(i=n0; i<n; i++){ 
        if(y[i]<=1.0){
            for(j=0;j<=m;j++) Bta[i+n*j]=1-pbeta(y[i], j+1, m-j+1, TRUE, FALSE);
            Bta[i+n*mp1]=1.0;}
        else for(j=0;j<=mp1;j++) Bta[i+n*j]=0.0; //Bta[i+n*mp1]=0.0;
    }
}
/*////////////////////////////////////////////////////*/
/*    Calculate fm(y,p) and Sm(y,p)                   */
/*////////////////////////////////////////////////////*/
void fm_Sm(double *p, int m, double *BSy, double *BSy2, int n, double *Sy, double *Sy2){
    int j, k, mp1=m+1;
    for(k=0; k<n; k++){
        Sy[k] = 0.0; 
        Sy2[k] = 0.0;
        for(j=0; j<=mp1; j++){
            Sy[k] += p[j]*BSy[k+n*j];
            Sy2[k] += p[j]*BSy2[k+n*j]; 
        }
    }
}    
/*////////////////////////////////////////////////////*/
/*  Initialize p for fm(.|x0;p) using                 */
/*            Bernstein approx fm0(.|x0;p0)           */
/*  p_i=Cm*fmt(i/m,p0), i=0:m, sum(p_i)=1-pt[m0+1]    */
/*////////////////////////////////////////////////////*/
void pm(double *p0, int m0, double *p, int m){
    int i, j;
    double tmp=0.0, pi0=1-p0[m0+1];
    for(i=0; i<=m; i++){
        p[i]=0.0;
        for(j=0; j<=m0; j++){
            p[i] += p0[j]*dbeta(i/(double)m, j+1, m0-j+1, FALSE); 
        }
        tmp+=p[i];
    }
    for(i=0; i<=m; i++) p[i]=pi0*p[i]/tmp;
    p[m+1]=p0[m0+1];
}  
/*////////////////////////////////////////////////////*/
/*   Calculate exp(g*x-tilde) for a fixed x0          */
/*           x-tilde = x-x0                           */
/*////////////////////////////////////////////////////*/
void egxmx0(double *gama, int d, double *x, int n, double *egx, double *x0){
    int i,j;
    double gx0=0.0;
    for(j=0;j<d;j++) gx0+= x0[j]*gama[j];
//    egx0 = exp(egx0);
    for(i=0;i<n;i++){
        egx[i]=0.0;
        for(j=0;j<d;j++) egx[i]+= x[i+n*j]*gama[j];
        egx[i]=exp(egx[i]-gx0);
    }
//    for(i=0;i<n;i++) egx[i] /= egx0;
}
/*////////////////////////////////////////////////////*/
/*       Calculate exp(g*x-tilde) and find            */
/*      x0=x(gama)=argmin{gamma'x_i: i=1:n}           */
/*////////////////////////////////////////////////////*/
void egx_x0(double *gama, int d, double *x, int n, double *egx, double *x0){
    int i,j;
    double egx0=0.0;
    for(j=0;j<d;j++) egx0 += x0[j]*gama[j];
    egx0=exp(egx0);
    for(i=0;i<n;i++){
        egx[i]=0.0;
        for(j=0;j<d;j++) egx[i]+= x[i+n*j]*gama[j];
        egx[i]=exp(egx[i]);
        if(egx[i]<egx0){
            egx0=egx[i];
            for(j=0;j<d;j++) x0[j] = x[i+n*j];
        }
    }
    for(i=0;i<n;i++) egx[i]/= egx0;
}

/*///////////////////////////////////////*/
/*    Exponential Change-point Method    */
/*   for choosing optimal model degree   */
/*///////////////////////////////////////*/
// chpt[0] input k, ouput change-point
void chpt_exp(double *lk, double *lr, double *pv, int *chpt){
    int i, k = chpt[0]; 
    double mLR=0.0, lr0, lnk=log(k), llnk=log(lnk);
    //chpt[0]=k;
    lr0 = k*log((lk[k]-lk[0])/(double)k);
    lr[k-1]=0.0;
    for(i=1;i<k;i++){
        lr[i-1] = lr0-i*log((lk[i]-lk[0])/(double)i)-(k-i)*log((lk[k]-lk[i])/(double)(k-i));
        if(lr[i-1]>mLR){
            chpt[0] = i; 
            mLR=lr[i-1];
        }
     //Rprintf("\n lr[%d]=%f\n",i, lr[i-1]);
    }
    // p-value of change-point
    pv[0]=1.0-exp(-2*lnk*lnk*sqrt(llnk*M_1_PI)*exp(-2*sqrt(mLR*llnk)));
    //Rprintf("\n pv[%d]=%f\n",0, pv[0]);
    //Rprintf("\n chpt=%d, p-val=%f\n", chpt[0], pv[0]);
}
/*///////////////////////////////////////*/
/*       Gamma Change-point Method       */
/*   for choosing optimal model degree   */
/*///////////////////////////////////////*/
// MLE theta[0:1]=(alpha, beta) of Gamma Model 
// and loglikelihood theta[2] based on x[i], k<= i < n.
void mle_gamma(double *x, int k, int n, double *res){
    int i, it=0, nmk=n-k;
    double s=1.0, d=0.0, xbar=0.0, tmp=0.0;
    double del=1.0, eps=1.0e-10, maxit=100;
    // initial: MME
    for(i=k;i<n; i++){
        xbar += x[i];
        s *= x[i];
        tmp += x[i]*x[i];
    }
    xbar /= (double) nmk;
    d = log(xbar) - log(s)/(double)nmk;
    tmp = (tmp-nmk*xbar*xbar)/(double)(nmk-1);
    tmp = xbar*xbar/tmp; 
    while(del>eps && it<maxit){
        del = tmp*(log(tmp)-digamma(tmp)-d)/(1.0-tmp*trigamma(tmp));
        tmp = tmp-del;
        del = fabs(del);
    }
    res[0] = tmp;
    res[1] = xbar/tmp;
    res[2] = (tmp-1)*s-nmk*(tmp*(log(res[1])+1)+lgammafn(tmp));
}    

void chpt_gamma(double *lk, double *lr, double *pv, int *chpt){
    int i, k=chpt[0]; 
    double mLR=0.0, *res, *x, lnk=log(k), llnk=log(lnk);
    res = Calloc(3, double);
    x = Calloc(k, double);
    for(i=0;i<k;i++) x[i]=lk[i+1]-lk[i];
    //chpt[0]=k;
    lr[k-1]=0.0;
    for(i=1;i<k;i++){
        mle_gamma(x, 0, i, res);
        lr[i-1] = res[2];
        mle_gamma(x, i, k, res);
        lr[i-1] += res[2];
        mle_gamma(x, 0, k, res);
        lr[i-1] -= res[2];
        if(lr[i-1]>mLR){
            chpt[0] = i; 
            mLR=lr[i-1];
        }
     //Rprintf("\n lr[%d]=%f\n",i, lr[i-1]);
    }
    // p-value of change-point
    pv[0]=1.0-exp(-2*lnk*lnk*llnk*exp(-2*sqrt(mLR*llnk)));
    //Rprintf("\n pv[%d]=%f\n",0, chpt[0]);
    //Rprintf("\n chpt=%d, p-val=%f\n", chpt[0], pv[0]);
    Free(x); Free(res);
}
//choosing optimal degree based on lk using gamma chage-point model
void optim_gcp(int *M, double *lk, double *lr, int *m, 
        double *pval, int *chpts){
    int i, *cp, k=M[1]-M[0];
    double *res; 
    cp = Calloc(1, int);
    res = Calloc(1, double);
    for(i=0;i<3;i++){            
        pval[i]=1.0;
        chpts[i]=i;
    }
    for(i=3;i<=k;i++){
        cp[0]=i;
        chpt_exp(lk, lr, res, cp);
        pval[i]=res[0];
        chpts[i]=cp[0];
    }
    m[0] = cp[0]+M[0];
    Free(cp); Free(res);
}
// End of file "mable-utility.c
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
/*////////////////////////////////////////////////////////////*/
/*   Log-Bernstein-Likelihood                                 */
/*  p=(p0, p1, ..., pm),                                      */
/*  x: n-vector, sample from distribution F with support [0,1]*/
/*  Bta: n x (m+1) matrix of Bernstein base polynomials       */
/*////////////////////////////////////////////////////////////*/
double loglik(double *p, double *Bta, int m, int n){
    int i,j;
    double llik, fx;
    llik = 0.0;
    for(i=0; i<n; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*Bta[i+n*j];
        }
        llik += log(fx);
    }
    return(llik);
}
/*//////////////////////////////////////////////////*/
/*                 EM Method                        */
/*//////////////////////////////////////////////////*/
/* EM Method for mixture of beta(i+1, m+1-i), i=0,...,m. */
void em_beta_mix(double *p, double *Bta, int m, int n, int maxit, double eps,
    double *llik, int *convergence, double *delta){
    int i, j, it;
    double del, llik_nu, *pBta,  *fp, *pnu;
    pBta = Calloc((m+1)*n, double);
    fp = Calloc(n, double);
    pnu = Calloc(m+1, double);
    llik[0] = loglik(p, Bta, m, n);
    del = 10.0;
    it = 0;
    convergence[0]=0;
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
        llik_nu = loglik(pnu, Bta, m, n);
        del = fabs(llik[0]-llik_nu);
        it += 1;
        for(i=0; i<=m; i++) p[i] = pnu[i];
        llik[0] = llik_nu;
    }
    if(it==maxit) {
        convergence[0]=1;
        delta[0]=del;}
    Free(pBta); Free(fp); Free(pnu);
}
/* end function em_beta_mix */
void mable_em(int *m, int *n, double *p, double *x, int *maxit,  double *eps,
    double *llik, int *convergence, double *delta){
    double *Bta;
    Bta = Calloc((*m+1)*(*n), double);
    dBeta(x, *m, *n, Bta);
    em_beta_mix(p, Bta, *m, *n,  *maxit, *eps, llik, convergence, delta);
    Free(Bta);
}
/*////////////////////////////////////////////////////////////////*/
/*   Log-Bernstein-Likelihood for Grouped Data                    */
/*   p=(p0, p1, ..., pm),                                         */
/*   (t,n):  grouped data from F with support [0,1]               */
/*   t: N-vector of breaks in [0,1], n: N-vector of frequencies   */
/*   dBta: (m+1) x N matrix of Beta[t(j+1)]-Beta[t(j)]            */
/*////////////////////////////////////////////////////////////////*/
double loglik_group(double *p, double *dBta, int *n, int m, int N){
    int i,j;
    double llik, fx;
    llik = 0.0;
    for(i=0; i<N; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*dBta[i+N*j];
        }
        llik += n[i]*log(fx);
    }
    return(llik);
}
/* Grouped data EM Method for mixture of beta(i+1, m+1-i), i=0,...,m. */
/* n = (n1,...,n_N), t=(t0, t1,...,t_N), */
void em_beta_mix_group(double *p, double *dBta, int N, int m, int *n, int maxit,
    double eps, double *llik, int *convergence, double *delta){
    int i, j, it;
    double del, llik_nu, *pBta,  *pnu, *fp, nn;
    pBta = Calloc((m+1)*N, double);
    fp = Calloc(N, double);
    pnu = Calloc(m+1, double);
    llik[0] = loglik_group(p, dBta, n, m, N);
    nn = 0;
    for(i=0; i<N;i++) nn+=n[i];
    del = 10.0;
    it = 0;
    convergence[0]=0;
    while(del>eps && it<maxit){
        for(i=0; i<N; i++){
            fp[i] = 0.0;
            for(j=0; j<=m; j++){
                pBta[i+N*j] = p[j]*dBta[i+N*j];
                fp[i] += pBta[i+N*j];
            }
        }
        for(j=0; j<=m; j++){
            pnu[j] = 0.0;
            for(i=0; i<N; i++) pnu[j] += n[i]*pBta[i+N*j]/fp[i];
            pnu[j] /= nn;
        }
        llik_nu = loglik_group(pnu, dBta, n, m, N);
        del = fabs(llik[0]-llik_nu);
        for(j=0; j<=m; j++) del+=fabs(p[j]-pnu[j]);
        it += 1;
        for(j=0; j<=m; j++) p[j] = pnu[j];
        llik[0] = llik_nu;
        //Rprintf("del=%f, it=%d, eps=%f\n", del, it, eps);
    }
    if(it==maxit){
        convergence[0]=1;
        delta[0]=del;
    }
    Free(pBta); Free(fp); Free(pnu);
}
/* end function em_beta_mix_group */

void mable_em_group(int *m, int *n, int *N, double *p, double *t, int *maxit,
    double *eps, double *llik, int *convergence, double *delta){
    double *dBta;
    dBta = Calloc((*m+1)*(*N), double);
    cpBeta(t, *m, *N, dBta);
    em_beta_mix_group(p, dBta, *N, *m, n,  *maxit, *eps, llik, convergence, delta);
    Free(dBta);
}
////////////////////////////////////////////////////////
// choosing optimal degree m by change-point method
//  using exponential change-point models
////////////////////////////////////////////////////////
void mable_optim(int *M, int *n, double *p, double *x, int *maxit,
        double *eps, double *lk, double *lr, int *optim, double *pval,  
        double *bic, int *chpts, double *tini, int *progress, 
        int *convergence, double *delta, double *level){
    int d, i, j, l, m, *cp, lp, tmp, *diverge, k=M[1]-M[0]; 
    double *phat, *Bta, *llik, *res;//, max_bic; 
    double pct=0.0, ttl;
    lp=M[0]*(k+1)+(k+1)*(k+2)/2;
    phat = Calloc(lp, double);
    diverge = Calloc(1, int);
    cp = Calloc(1, int);
    res = Calloc(1, double);
    llik = Calloc(1, double);
    Bta = Calloc((M[1]+1)*(*n), double);
    ttl=(double) (k+2)*(k+1);
    m=M[0];
    for(j=0; j<=m; j++) p[j]=1.0/(double) (m+1);
    dBeta(x, m, *n, Bta);
    em_beta_mix(p, Bta, m, *n,  *maxit, eps[0], llik, diverge, delta);
    convergence[0]+=diverge[0];
    tmp=m+1;//i=0; tmp=(i+1)*M[0]+(i+1)*(i+2)/2;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    lk[0] = *llik;
    pct += 2.0;
    pval[0]=1.0;
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;
    bic[0]=lk[0]-.5*d*log(*n);
//    max_bic=bic[0];
//    optim[1]=0;
//    Rprintf("\n optim[0]=%d, ni=%d, max_bic=%f, bic[%d]=%f\n", optim[0], ni, max_bic, 0, bic[0]);
    chpts[0]=0;
    if(*progress==1) ProgressBar(pct/ttl,"");
    i=1;
    //stop if the p-value for change-point is small 
    while(i<=k && pval[i-1]>*level){
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
        m=M[0]+i;
        //make sure initial p is in the interior of the simplex
        for(j=0; j<=m; j++) p[j] =(p[j]+ *tini/(double)(m+1.0))/(1.0+*tini);
        em_beta_mix(p, Bta, m, *n, *maxit, eps[0], llik, diverge, delta);
        convergence[0] += diverge[0];
        for(j=0;j<=m;j++) phat[j+tmp]=p[j];
        tmp += m+1;
        lk[i] = *llik;
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        // Calculate BIC
        d=0;
        for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
        d=d-1;
        bic[i]=lk[i]-.5*d*log(*n);
//        if(bic[i]>max_bic){
//            max_bic=bic[i];
//            optim[1]=i;
//        }
//     Rprintf("\n optim[0]=%d, max_bic=%f, bic[%d]=%f\n", cp[0], max_bic, i, bic[i]);
//     Rprintf("\n chpt=%d, p-val=%f\n", chpts[i], pval_cp[i]);
        pct += 2*(i+1);
        if(*progress==1) ProgressBar(pct/ttl,"");
//        clockProgress(pct); 
        i++; 
    }
    if(convergence[0]>0) convergence[0]=1;
    if(*progress==1){
        ProgressBar(1.0,"");
        Rprintf("\n");}
    if(m==M[1]){
        convergence[0]+=1; 
        warning("The maximum candidate degree has been reached \n with a p-value of the change-point %f.\n", res[0]);}
    M[1]=m;
    tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
    optim[0]=cp[0]+M[0];
    m=optim[0];
    for(j=0;j<=m;j++) p[j]=phat[tmp+j];
//    tmp=optim[1]*(M[0]*2+(optim[1]+1))/2;
//    optim[1]+=M[0];
//    m=optim[1];
//    for(j=0;j<=m;j++) p[optim[0]+1+j]=phat[tmp+j];
    Free(phat); Free(Bta); Free(llik); 
    Free(cp); Free(res);
}
/*//////////////////////////////////////////////////////////*/
void mable_optim_group(int *M, int *N, double *p, double *t, int *n, int *maxit, double *eps, 
        double *lk, double *lr, int *optim, int *progress, int *convergence, double *delta, 
        double *tini, double *bic, double *pval, int *chpts, double *level){
    int d, nn=0, i, j, m, k, lp, tmp, *cp, *diverge; 
    double *phat, *dBta, pct=0.0, *mlik, ttl, *res;//, max_bic;//maxlr, tmp, 
    for(i=0; i<*N;i++) nn+=n[i];
    k=M[1]-M[0];
    lp=M[0]*(k+1)+(k+1)*(k+2)/2;
    cp=Calloc(1, int);
    diverge=Calloc(1, int);
    dBta = Calloc((M[1]+1)*(*N), double);
    phat = Calloc(lp, double);
    res = Calloc(2, double);
    ttl = (double)((k+2)*(k+1));
    mlik = Calloc(1, double);
    m=M[0];
    for(j=0; j<=m; j++) p[j]=1.0/(double)(m+1);
    cpBeta(t, m, *N, dBta);
    em_beta_mix_group(p, dBta, *N, m, n, *maxit, eps[0], mlik, diverge, delta);
    convergence[0]+=diverge[0];
    tmp=m+1;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    pct += 2/ttl;
    if(*progress==1) ProgressBar(pct,""); 
    lk[0] = mlik[0];
    pval[0]=1.0;
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;
    bic[0]=lk[0]-.5*d*log(nn);
//    max_bic=bic[0];
//    optim[1]=0;
    chpts[0]=0;
    i=1;
    while(i<=k && pval[i-1]>*level){
        // updating p
        p[m+1] = (m+1)*p[m]/(double)(m+2);
        for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
        p[0] = (m+1.)*p[0]/(double)(m+2.);
        m=M[0]+i;
        cpBeta(t, m, *N, dBta);
        // make sure initial p is in the interior of the simplex
        for(j=0; j<=m; j++) p[j] =(p[j]+ *tini/(double)(m+1.0))/(1.0+*tini);
//        for(j=0;j<=m;j++) p[j]=1.0/(double) (m+1);
        em_beta_mix_group(p, dBta, *N, m, n, *maxit, eps[0], mlik, diverge, delta);
        convergence[0]+=diverge[0];
        for(j=0;j<=m;j++) phat[j+tmp]=p[j];
        tmp+=m+1;
        lk[i] = mlik[0];
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        // Calculate BIC
        d=0;
        for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
        d=d-1;
        bic[i]=lk[i]-.5*d*log(nn);
//        if(bic[i]>max_bic){
//            max_bic=bic[i];
//            optim[1]=i;
//        }
//     Rprintf("\n optim[0]=%d, max_bic=%f, bic[%d]=%f\n", cp[0], max_bic, i, bic[i]);
        pct +=2*(i+1)/ttl;
        if(*progress==1)ProgressBar(pct,"");
        i++;
    }
    if(*progress==1){
        ProgressBar(1.0,"");
        Rprintf("\n");}
    if(convergence[0]>0) convergence[0]=1;
    if(m==M[1]){
        convergence[0]+=1; 
        warning("The maximum candidate degree has been reached \n with a p-value of the change-point %f.\n", res[0]);}
    M[1]=m;
    tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
    optim[0]=cp[0]+M[0];
    m=optim[0];
    for(j=0;j<=m;j++) p[j]=phat[tmp+j];
//    tmp=optim[1]*(M[0]*2+(optim[1]+1))/2;
//    optim[1]+=M[0];
//    m=optim[1];
//    for(j=0;j<=m;j++) p[optim[0]+1+j]=phat[tmp+j];
    Free(dBta); Free(phat); Free(diverge); Free(cp); Free(res);
    Free(mlik);
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
// End of file "mable.c
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*  Maximum Approximate Bernstein likelihood Estimation   */
/*  in Accelerated Failure Time Regression model based    */
/*                  Interval Censored data                */
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/*                                                            */
/*  Reference:                                                */
/*   Zhong Guan, Semiparametric Maximum Likelihood Estimation */
/*         in Accelerated Failure Time Model for              */
/*               Interval-Censored Data                       */
/*                                                            */
/*////////////////////////////////////////////////////////////*/
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
/*/////////////////////////////////////////////////////////////////////*/
/*         MABLE for AFT Model: S(t|x)=S[t exp(-gamma'x)|0]             */
/*  Maximum Approximate Bernstein Likelihood Estimation of survival    */
/*     survival function S(t|x) and regression coefficients gamma      */
/*       based on interval censored data (y=(y1,y2), x, delta)         */
/* (y1,y2): censoring interval containing event time t                 */
/*       x: d-vector of covariate values                               */
/*   delta: censoring indicator, 0 uncensored, 1: interval censored    */
/*  The data are arranged so that the 1st n0 obs are uncensored (y1=y2)*/
/*(delta=0) and the rest n1=n-n0 obs are interval censored(delta=1)    */
/*/////////////////////////////////////////////////////////////////////*/
/*/////////////////////////////////////////////////////////////////*/
/*   Log-Likelihood ell(gamma, p) for AFT model, where             */
/*   gamma=(gamma1,...,gamma_d), p=(p0, p1, ..., pm, p_{m+1}),     */
/*      gx: gamma*x.tilde, emgx=exp(-gx)                           */
/*     BSz: 1-B(y1*emgx)                                           */
/*    BSz2: (beta(y1*emgx),1-B(y2*emgx)),                          */
/*/////////////////////////////////////////////////////////////////*/

//log_blik_aft(p, m, gx, n0, n1, betay, BSz, BSz2);
double log_blik_aft(double *p, int m, double *gx, int n0, int n1, 
            double *BSz, double *BSz2){
    int i,j, n=n0+n1;
    double llkhd, fz, dSz;
    llkhd = 0.0;
    for(i=0; i<n0; i++){
        fz = 0.0; 
        for(j=0; j<=m; j++){
            fz += p[j]*BSz2[i+n*j];
        }
        llkhd += -gx[i]+log(fz);
    }
    for(i=n0; i<n; i++){
        dSz=0.0;
        for(j=0; j <= m; j++){
            dSz += p[j]*(BSz[i+n*j]-BSz2[i+n*j]); 
        }
        llkhd += log(dSz);
    }
    return(llkhd);
}
/*/////////////////////////////////////////////////////////////////*/
/* Derivatives of loglikelihood ell(gamma, p) wrt gamma, aft model */
/*/////////////////////////////////////////////////////////////////*/
// 
void logblik_aft_derv(double *gama, double *p, int d, int m, double *y, double *y2,
            double *x, double *x0, int n0, int n1, double *ell, double *dell, double *ddell){
    int i,j,k, n=n0+n1, mp2=m+2;
    double egxt, tmp1=0.0, tmp2=0.0, A, B, C;
    double *BSz, *BSz2, *bz, *bz2, *z, *z2;   
    z = Calloc(n, double);
    z2 = Calloc(n, double);
    bz = Calloc(n*mp2, double);  
    bz2 = Calloc(n*mp2, double); 
    BSz = Calloc(n*mp2, double);  
    BSz2 = Calloc(n*mp2, double); 
    ell[0]=0.0;
    for(i=0; i<d; i++){ dell[i]=0.0; 
        for(j=0; j<d; j++) ddell[i+d*j]=0.0;}
    for(k=0; k<n; k++){
        egxt=0.0;
        for(i=0; i<d; i++) egxt -= gama[i]*(x[k+n*i]-x0[i]);
        if(k<n0) ell[0] += egxt;
        egxt= exp(egxt);
        z[k] = egxt*y[k];
        z2[k] = egxt*y2[k];
        // add error check whether z[k] and z2[k] <=1????
//Rprintf("\n y[%d]=%f\n, z[%d]=%f", k, y[k], k, z[k]);
//        if (y2[k]<=1 && z2[k]>1) error("\n Error: z2 >1 \n");
    }
    Bdata(z, m, 0, n, BSz);//1-B(z),1-B(z2) 
    Bdata(z2, m, 0, n, BSz2);
    Bdata(z, m, n, 0, bz); // beta(z), beta(z2)
    Bdata(z2, m, n, 0, bz2);
    for(k=0; k<n0; k++){
        A = 0.0;
        B = 0.0;
        C = 0.0;
        for(j=0;j<m;j++){
            A += p[j]*bz[k+n*j];
            B += p[j]*(j*bz[k+n*j]-(j+1)*bz[k+n*(j+1)]);
            C += p[j]*(j*j*bz[k+n*j]-(j+1)*(2*j+1)*bz[k+n*(j+1)]+(j+1)*(j+2)*bz[k+n*(j+2)]);
        }  
        A += p[m]*bz[k+n*m]; 
        B += p[m]*m*bz[k+n*m];
        C += p[m]*m*m*bz[k+n*m];       
        ell[0] += log(A);
        tmp1 = B/A;
        tmp2 = C/A;
        for(i=0; i<d; i++){
            dell[i] -= (1+tmp1)*(x[k+n*i]-x0[i]); 
            for(j=0; j<d; j++)
                ddell[i+d*j] -= (tmp1*tmp1-tmp2)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
        }
    }
    for(k=n0; k<n; k++){
        A = 0.0;
        B = 0.0;
        C = 0.0;
        tmp1 = 0.0;
        for(j=0;j<=m;j++){
            A += p[j]*(z[k]*bz[k+n*j]-z2[k]*bz2[k+n*j]);
            tmp1 += p[j]*(BSz[k+n*j]-BSz2[k+n*j]);
            B += p[j]*(j*bz[k+n*j]-(j+1)*bz[k+n*(j+1)]);
            C += p[j]*(j*bz2[k+n*j]-(j+1)*bz2[k+n*(j+1)]);
//Rprintf("\n p[%d]=%f", j, p[j]);
//Rprintf("\n bz[%d,%d]=%f", k, j, bz[k+n*j]);
        }
        tmp2 = z[k]*B-z2[k]*C;
        ell[0] += log(tmp1);
        for(i=0; i<d; i++){
            dell[i]+=A*(x[k+n*i]-x0[i])/tmp1;           
            for(j=0;j<d;j++){
                ddell[i+d*j]-=(A/tmp1)*(A/tmp1+1)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
                ddell[i+d*j]-=(tmp2/tmp1)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
            }
        }
    }
//    for(i=0; i<d; i++)
//       for(j=0;j<d;j++) Rprintf("\n ddell[%d,%d]=%f", i,j, ell[i+d*j]);

    Free(z); Free(z2); Free(bz); Free(bz2);
    Free(BSz); Free(BSz2);
}
/*////////////////////////////////////////////////////*/
/*  Maximizer gamma of ell(gamma, p) for a gvien p    */
/*                   for AFT model                    */
/*////////////////////////////////////////////////////*/

void gofp_aft(double *gama, int d, double *p, int m, double *y, double *y2, 
            double *x, double *x0, int n0, int n1, double *ell, double *dell, 
            double *ddell, double eps, int maxit, int prog){
    int i,j, it=0;
    double del=0.0, *tmp; 
    tmp = Calloc(d, double);
    logblik_aft_derv(gama, p, d, m, y, y2, x, x0, n0, n1, ell, dell, ddell);
    for(i=0;i<d;i++) del+=fabs(dell[i]); 
    while(it<maxit && del>eps){
        minverse(ddell, d);  
        for(i=0;i<d;i++){
            tmp[i] = 0.0;
            for(j=0;j<d;j++) tmp[i] += ddell[i+d*j]*dell[j];
        }
        del = 0.0;
        for(i=0;i<d;i++){
            gama[i] -= tmp[i];
            del += fabs(tmp[i]);
        }
        logblik_aft_derv(gama, p, d, m, y, y2, x, x0, n0, n1, ell, dell, ddell);
        for(i=0;i<d;i++) del+=fabs(dell[i]);
        it++;
        R_CheckUserInterrupt();
    }
    if(prog==0) Rprintf("NT: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del, ell[0]);
    Free(tmp); 
}


/*////////////////////////////////////////////////////*/
/*   maximizer p of ell(gamma, p) for a gvien gamma   */
/*  gx:  gama*x.tilde, where gama is the given        */
/*       regression coefficient of the AFT model      */
/*  ss: step-size epsilon:  p = (1-ss)p+ss Psi(p),    */
/*        default ss=1 so that p = Psi(p)             */
/*////////////////////////////////////////////////////*/
void pofg_aft(double *p, int m, double *gx, int n0, int n1, double *BSz, double *BSz2, 
        double *llik, double eps, int maxit, int prog, int *conv, double *delta){
    int i, j, n=n0+n1, mp1=m+1, it=0;
    double  del=1.0, dSz;
    double *Tmp, *Tmp2, *pnu, llik_nu;
    Tmp=Calloc(mp1, double); // Tmp's can be independent of i ???
    Tmp2=Calloc(mp1, double);
    pnu=Calloc(mp1, double);
    llik[0]=log_blik_aft(p, m, gx, n0, n1, BSz, BSz2); 
    conv[0]=0;
    while(del>eps && it<maxit){
        for(j=0;j<=m;j++) pnu[j]=0.0;
        // p = p *Psi(p) 
        for(i=0; i<n0;i++){
            dSz=0.0; 
            for(j=0;j<=m;j++) {
                Tmp2[j]=BSz2[i+n*j]*p[j];
                dSz+=Tmp2[j];
            }
            for(j=0;j<=m;j++) {
                pnu[j]+=Tmp2[j]/dSz;
            }
        }
        for(i=n0; i<n;i++){
            dSz=0.0;  
            for(j=0;j<=m;j++){
                Tmp[j]=BSz[i+n*j]*p[j];
                Tmp2[j]=BSz2[i+n*j]*p[j];
                dSz+=Tmp[j]-Tmp2[j];
            }
            for(j=0;j<=m;j++){
                pnu[j]+=(Tmp[j]-Tmp2[j])/dSz;
            }
        }  
        for(j=0;j<=m;j++) pnu[j] /= (double) n;
        //pnu<-(1-ss)*p+ss*pnu
        llik_nu=log_blik_aft(pnu, m, gx, n0, n1, BSz, BSz2); 
        del=fabs(llik[0]-llik_nu);
        it++;  
        llik[0]=llik_nu;
        for(j=0;j<=m;j++) p[j]=pnu[j];
        R_CheckUserInterrupt();
    }
    if(prog==0) Rprintf("EM: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del, llik[0]);
    if(it==maxit){
        conv[0]+=1;
        delta[0]=del;}
    Free(Tmp);  Free(Tmp2); Free(pnu);
}

/*//////////////////////////////////////////////////////////////*/
/* Maximum Approximate Profile Likelihood Esimation of AFT model*/
/*  Select optimal degree m with a given gamma for AFT model    */
/* M: set of positive integers as candidate degrees of          */ 
/*       Bernstein poly model                                   */
/* gama: an efficient estimate of regression coefficient        */ 
/*       gamma, for data without covariate we set gama=0        */
/*   x0: baseline covariate value, default is 0                 */
/*    x: d-dim covariate                                        */
/*//////////////////////////////////////////////////////////////*/
void mable_aft_gamma(int *M, double *gama, int *dm, double *x, double *y, double *y2, int *N, 
        double *x0, double *lk, double *lr, double *p, double *ddell, double *eps, int *maxit, 
        int *progress, double *pval, int *chpts, double *level, int *conv, double *delta){        
    int i,j, d=dm[0], k=M[1]-M[0], *cp, n0=N[0], n1=N[1], n=n0+n1, tmp, itmp=0;
    int m=M[1], mp1=m+1,  mp2=m+2, lp=(k+1)*M[0]+(k+1)*(k+2)/2;//optim=m,
    double tini=.0001, pct=0.0, ttl, *z, *z2, *res;
    double *ell, *dell, *gx, *BSz, *BSz2, *phat; 
    cp = Calloc(1, int);
    res = Calloc(1, double);
    phat=Calloc(lp, double);
    ell = Calloc(1, double);
    dell = Calloc(d, double);
    BSz = Calloc(n*mp2, double);  
    BSz2 = Calloc(n*mp2, double); 
    z = Calloc(n, double);
    z2 = Calloc(n, double);
    gx = Calloc(n, double);
    if(*progress==1) {Rprintf("\n Mable fit of AFT model with the given gamma ... \n");
        ProgressBar(0.0,""); }
    ttl = (double)((k+2)*(k+1));
    egxmx0(gama, d, x, n, gx, x0);
    for(i=0;i<n;i++) {
        z[i] = y[i]/gx[i];
        z2[i] = y2[i]/gx[i];
        gx[i] = log(gx[i]);
        // add check to see any z, z2 are bigger than 1
        if (y2[i]<=1 && z2[i]>1) {
            Rprintf("\n");
            error("Try another baseline 'x0' and/or a larger truncation time 'tau'.\n");}
    }
    m=M[0]; 
    mp1=m+1;
    mp2=m+2;
    Bdata(z, m, 0, n, BSz);
    Bdata(z2, m, n0, n1, BSz2);
    for(i=0;i<=m;i++) p[i]=1.0/(double) mp1; 
    pofg_aft(p, m, gx, n0, n1, BSz, BSz2, ell, *eps, *maxit, *progress, conv, delta);
    itmp+=conv[0];
    tmp=mp1;
    for(i=0;i<mp1;i++) phat[i]=p[i];
    lk[0]=ell[0]; 
    pval[0]=1.0;
    chpts[0]=0;
    pct += 2/ttl;
    if(*progress==1) ProgressBar(pct,""); 
    i=1;
    while(i<=k && pval[i-1]>*level){
        p[mp1] = mp1*p[m]/(double)mp2;
        for(j=m; j>0; j--) p[j] = (p[j-1]*j+p[j]*(mp1-j))/(double)mp2;
        p[0] = mp1*p[0]/(double)mp2;
        m=M[0]+i;
        mp1=m+1;
        mp2=m+2;
        Bdata(z, m, 0, n, BSz);
        Bdata(z2, m, n0, n1, BSz2);
        for(j=0;j<mp1;j++) p[j]=(p[j]+tini/(double) mp1)/(1.0+tini);
        pofg_aft(p, m, gx, n0, n1, BSz, BSz2, ell, *eps, *maxit, *progress, conv, delta);
        for(j=0;j<mp1;j++) phat[j+tmp]=p[j];
        tmp+=mp1;
        lk[i]=ell[0];
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        R_CheckUserInterrupt();
        pct +=2*(i+1)/ttl;
        if(*progress==1) ProgressBar(pct,""); 
        i++;
        itmp+=conv[0];
    }
    if(*progress==1){
        ProgressBar(1.00,"");
        Rprintf("\n");}
    // Rprintf("mable-aft done!\n");
    if(itmp>0) conv[0]=1; 
    else conv[0]=0;
    if(k>0){
        if(m==M[1]){
            conv[0]+=1; 
            warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);
            delta[0]=res[0];
        }
    }
    M[1]=m;
    tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
    dm[1]=cp[0]+M[0];
    m=dm[1];
    for(j=0;j<=m;j++) p[j]=phat[tmp+j];
    Free(phat);  Free(ell);  Free(dell);
    Free(BSz); Free(BSz2); Free(z); Free(z2); 
    Free(gx);  Free(cp); Free(res);
}
/*////////////////////////////////////////////////////*/
/*  Maximum approx. Bernstein likelihood estimate of  */
/*   (gamma, p) with a fixed degree m for AFT model   */
/*////////////////////////////////////////////////////*/
void mable_aft_m(double *gama, double *p, int *dm, double *x, double *y, double *y2, 
            int *N, double *x0, double *ell, double *ddell, double *EPS, int *MAXIT,
            int *progress, int *conv, double *delta){
    int i, n0=N[0], n1=N[1], n=n0+n1, d=dm[0], m=dm[1],  mp2 =m+2, it=0;
    int maxit=MAXIT[0], maxit_em=MAXIT[1], prog=1; 
    double eps=EPS[0], eps_em=EPS[1], pct=0.0;
    double *z, *z2, *BSz, *BSz2, *gnu, del=1.0;//*xt, tini, 
    double *ell1, *dell, *gx, *tmp;
    tmp = Calloc(d, double);
    //tini=0.00001;// tini is used to make sure p is in interior of S(m+1)
    ell1 = Calloc(1, double);  
    dell = Calloc(d, double);  
    z = Calloc(n, double);  
    z2 = Calloc(n, double);  
    gx = Calloc(n, double);  
    BSz = Calloc(n*mp2, double);  
    BSz2 = Calloc(n*mp2, double);  
    gnu = Calloc(d, double);  
    egxmx0(gama, d, x, n, gx, x0);
    for(i=0;i<n;i++){
        z[i] = y[i]/gx[i];
        z2[i] = y2[i]/gx[i];
        gx[i] = log(gx[i]);
        if (y2[i]<=1 && z2[i]>1) {
            Rprintf("\n");
            error("Try another baseline 'x0' and/or a larger truncation time 'tau'.\n");}
    }
    Bdata(z, m, 0, n, BSz);
    Bdata(z2, m, n0, n1, BSz2);
    if(*progress==1) ProgressBar(pct,""); 
    pofg_aft(p, m, gx, n0, n1, BSz, BSz2, ell, eps_em, maxit_em, prog, conv, delta);
//    Rprintf("\n ell=%f, ell0=%f\n",ell[0], ell[1]);
    while(it<maxit && (del>eps || ell[0]<ell[1])){
        gofp_aft(gama, d, p, m, y, y2, x, x0, n0, n1, ell1, dell, ddell, eps, maxit, prog);
        egxmx0(gama, d, x, n, gx, x0);
        for(i=0;i<n;i++) {
            z[i] = y[i]/gx[i];
            z2[i] = y2[i]/gx[i];
            gx[i] = log(gx[i]);
            if (y2[i]<=1 && z2[i]>1) {
                Rprintf("\n");
                error("Try another baseline 'x0' and/or a larger truncation time 'tau'.\n");}
        }
        Bdata(z, m, 0, n, BSz);
        Bdata(z2, m, n0, n1, BSz2);
        pofg_aft(p, m, gx, n0, n1, BSz, BSz2, ell1, eps_em, maxit_em, prog, conv, delta);
        del = fabs(ell1[0]-ell[0]);
        ell[0]=ell1[0];
//         Rprintf("\n ell=%f, ell0=%f\n",ell[0], ell1[0]);
//         Rprintf("\n pmp1=%d, x0=%f,  %f\n", *pmp1, x0[0], x0[1]);
//        pct=fmax(1-fabs(del-eps)/(.00001+eps),1-fabs(ell[1]-ell[0])/fabs(ell[1]));
        pct=fmin(it/(double)maxit, fmax(0.0, 1-fabs(del-eps)));
        if(*progress==1) ProgressBar(pct,""); 
        it++;
        R_CheckUserInterrupt();
        // Rprintf("         mable-m: it=%d, del=%f, ell=%f\n", it, del, ell[0]);
    }
    if(*progress==1){
        ProgressBar(1.0,""); 
        Rprintf("\n");}
    delta[0]=del;
    if(it==maxit){
        conv[0] = 1;
     //  warning("\nThe maximum iterations were reached \nwith a delta = %f.\n", del);
     }
    // Rprintf("mable-m: it=%d, del=%f\n", it, del);
    minverse(ddell, d); //Sig=-n*ddell
    Free(tmp); Free(z); Free(z2); Free(BSz); Free(BSz2); 
    Free(gnu); Free(ell1); Free(dell); Free(gx); 
}
/*///////////////////////////////////////////////////////////*/
/*  MABLE of (gamma, p) and an optimal degree m              */
/*         for AFT model                                     */
/*       M: set of positive integers as candidate degrees    */
/*            of Bernstein poly model, M=m0:m1, k=m1-m0      */
/*gama_hat: k-vector, initial values of gamma, an efficient  */
/*           estimate of regression coefficient gamma        */
/*    phat: (m1+1)-vector, phat[0:m0] is initial of p        */
/*            for degree m=m0, used to return phat           */
/*     dk: (d,k)                                             */
/*      x: d-dim covariate                                   */
/*     x0: baseline covariate value, default is 0            */
/*///////////////////////////////////////////////////////////*/
void mable_aft(int *M, double *gama, int *dm, double *p, double *x, double *y,    
        double *y2, int *N, double *x0, double *lk, double *lr, double *ddell, double *EPS, 
        int *MAXIT, int *progress, double *pval, int *chpts, double *level, int *conv){
    int i, j, d=dm[0], k=M[1]-M[0], tmp=0,*cp;//, n0=N[0], n1=N[1], n=n0+n1
    int m=M[1], mp1=m+1, lp=(k+1)*M[0]+(k+1)*(k+2)/2, prg=1-*progress; 
    double *phat, *ghat, *ell, pct=0.0, ttl, lnn=-1.0e20, *res;  //maxLR=0.0, lr0, 
    cp = Calloc(1, int);
    res = Calloc(1, double);
    phat=Calloc(lp, double);
    ghat=Calloc(d*(k+1), double);
    ell=Calloc(2, double);
//    egx=Calloc(n, double);
    if(*progress==1) {Rprintf("\n Mable fit of AFT model ... \n");
        ProgressBar(0.0,""); }
    ttl=(double) (k+2)*(k+1);
    m=M[0]; 
    mp1=m+1;
    dm[1]=m;
    for(i=0;i<mp1;i++) p[i]=1.0/(double)mp1;
    ell[1]=lnn;
    mable_aft_m(gama, p, dm, x, y, y2, N, x0, ell, ddell, EPS, MAXIT, &prg, conv, res);
    tmp=mp1;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    for(i=0;i<d;i++) ghat[i]=gama[i];
    lk[0]=ell[0];
    ell[1]=ell[0];
    pval[0]=1.0;
    chpts[0]=0;
    pct += 2;
    if(*progress==1) ProgressBar(pct/ttl,""); 
    i=1;
    while(i<=k && pval[i-1]>*level){
        p[mp1] = mp1*p[m]/(double)(m+2);
        for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(mp1-j))/(double)(m+2);
        p[0] = mp1*p[0]/(double)(m+2);
        m=M[0]+i; 
        dm[1]=m;
        mp1=m+1; 
        for(j=0;j<mp1;j++) p[j]=(p[j]+.000001/(double)mp1)/1.000001;
        mable_aft_m(gama, p, dm, x, y, y2, N, x0, ell, ddell,  EPS, MAXIT, &prg, conv, res);
        for(j=0;j<mp1;j++) phat[j+tmp]=p[j];
        tmp+=mp1;
        for(j=0;j<d; j++) ghat[j+i*d]=gama[j];
        lk[i]=ell[0]; 
        ell[1]=ell[0]; 
        // Rprintf("\n lk[%d]=%f, lpt=%d, temp=%d\n",i, lk[i], lpt, tmp);
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        R_CheckUserInterrupt();
        pct +=2*(i+1);
        if(*progress==1) ProgressBar(pct/ttl,""); 
        i++;
    }
    if(*progress==1){
        //ProgressBar(1.00,"");
        Rprintf("\n");}
    // Rprintf("mable-aft done!\n"); 
    if(m==M[1]){
        conv[0]+=1; 
        warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);
    }
//    else conv[0]=0;
    M[1]=m;
    tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
    dm[1]=cp[0]+M[0];
    m=dm[1];
    for(j=0;j<=m;j++) p[j]=phat[tmp+j];
    for(j=0; j<dm[0]; j++) gama[j]=ghat[dm[0]*cp[0]+j];
    Free(cp); Free(phat); Free(ghat);
    Free(res); Free(ell); 
    //Free(egx);
}
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*  Maximum Approximate Bernstein likelihood Estimation   */
/*     in Proportional Hazards Regression model based     */
/*                  Interval Censored data                */
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/*                                                            */
/*  Reference:                                                */
/*   Zhong Guan, Maximum Approximate Bernstein Likelihood     */
/*         Estimation in Proportional Hazard Model            */
/*               for Interval-Censored Data                   */
/*                                                            */
/*////////////////////////////////////////////////////////////*/
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
/*/////////////////////////////////////////////////////////////////////*/
/*         MABLE for PH Model: S(t|x)=S(t|0)^exp[gamma'x]              */
/*  Maximum Approximate Bernstein Likelihood Estimation of survival    */
/*     survival function S(t|x) and regression coefficients gamma      */
/*       based on interval censored data (y=(y1,y2), x, delta)         */
/* (y1,y2): censoring interval containing event time t                 */
/*       x: d-vector of covariate values                               */
/*   delta: censoring indicator, 0 uncensored, 1: interval censored    */
/*  The data are arranged so that the 1st n0 obs are uncensored (y1=y2)*/
/*(delta=0) and the rest n1=n-n0 obs are interval censored(delta=1)    */
/*/////////////////////////////////////////////////////////////////////*/
/*/////////////////////////////////////////////////////////////////*/
/*   Log-Likelihood ell(gamma, p) for ph model, where              */
/*   gamma=(gamma1,...,gamma_d), p=(p0, p1, ..., pm, p_{m+1}),     */
/*     egx: exp(gamma*x.tilde)                                          */
/*     BSy: 1-B(y1)                                                */
/*    BSy2: (beta(y1),1-B(y2)),                                    */
/*/////////////////////////////////////////////////////////////////*/
//log_blik_ph(p, m, egx, n0, n1, BSy, BSy2);
double log_blik_ph(double *p, int m, double *egx, int n0, int n1, 
            double *BSy, double *BSy2){
    int i,j, n=n0+n1, mp1=m+1;
    double llkhd, fy, Sy, Sy2;
    llkhd = 0.0;
    for(i=0; i<n0; i++){
        fy = 0.0; 
        Sy=0.0;
        for(j=0; j<=m; j++){
            fy += p[j]*BSy2[i+n*j];
            Sy += p[j]*BSy[i+n*j];
        }
        Sy+=p[mp1];
        llkhd += log(egx[i]*fy)+(egx[i]-1.0)*log(Sy);
    }
//    Rprintf("lk1: lk=%f\n", llkhd);
    for(i=n0; i<n; i++){
        Sy=0.0;
        Sy2=0.0;
        for(j=0; j<=mp1; j++){
            Sy += p[j]*BSy[i+n*j];
            Sy2 += p[j]*BSy2[i+n*j]; 
        }
//        Rprintf("Sy: Sy=%f\n", Sy);
        llkhd += log(R_pow(Sy, egx[i])-R_pow(Sy2, egx[i]));
    }
//    Rprintf("lk2: lk=%f\n", llkhd);
    return(llkhd);
}
/*/////////////////////////////////////////////////////////////////*/
/*      Derivatives of loglikelihood ell(gamma, p) wrt gamma       */
/*/////////////////////////////////////////////////////////////////*/
//logblik_ph_derv(gama, d, x, x0, n0, n1, betay, Sy, Sy2);
void logblik_ph_derv(double *gama, int d, double *x, double *x0, int n0, int n1,
            double *Sy, double *Sy2, double *ell, double *dell, double *ddell){
    int i,j,k, n=n0+n1;
    double egxt, tmp=0.0;
    double dPy, dPy2, Py12;
    ell[0]=0.0;
    for(i=0; i<d; i++){ dell[i]=0.0; 
        for(j=0; j<d; j++) ddell[i+d*j]=0.0;}
    for(k=0; k<n0; k++){
        egxt=0.0;
        for(i=0; i<d; i++) egxt+= gama[i]*(x[k+n*i]-x0[i]);
        egxt= exp(egxt);
        ell[0] += log(egxt*Sy2[k])+(egxt-1.0)*log(Sy[k]);
        for(i=0; i<d; i++){
            dell[i] += (1+egxt*log(Sy[k]))*(x[k+n*i]-x0[i]); 
            for(j=0; j<d; j++)
                ddell[i+d*j]+=egxt*log(Sy[k])*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
        }
    }
    for(k=n0; k<n; k++){
        egxt=0.0;
        for(i=0; i<d; i++) egxt+= gama[i]*(x[k+n*i]-x0[i]);
        egxt= exp(egxt);
//        tmp=0.0;
        dPy=R_pow(Sy[k], egxt); 
        dPy2=R_pow(Sy2[k], egxt);  
        Py12 =dPy-dPy2;
        ell[0] += log(Py12);
        dPy*=log(Sy[k]); 
        if(Sy2[k]>0) dPy2*=log(Sy2[k]); 
        else dPy2=0.0;
        tmp=dPy-dPy2; 
        if(Sy2[k]>0) tmp+=egxt*(dPy*log(Sy[k])-dPy2*log(Sy2[k]));
        else tmp+=egxt*dPy*log(Sy[k]);
        tmp-=egxt*(dPy-dPy2)*(dPy-dPy2)/Py12;
        for(i=0; i<d; i++){
            dell[i]+=egxt*(x[k+n*i]-x0[i])*(dPy-dPy2)/Py12;           
            for(j=0;j<d;j++) 
                ddell[i+d*j]+=egxt*tmp*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j])/Py12;
        }
    }
}

/*////////////////////////////////////////////////////*/
/*  Maximizer gamma of ell(gamma, p) for a gvien p    */
/*  with f(t|x0) being approximated by Bernstein poly */
/*////////////////////////////////////////////////////*/

void gofp_ph(double *gama, int d, double *p, int m, double *x, double *x0, int n0, int n1,  
            double *BSy, double *BSy2, double *ell, double *dell, double *ddell, 
            double eps, int maxit, int prog){
    int i,j, it=0, n=n0+n1;
    double del=0.0, *tmp, *Sy, *Sy2;
    tmp = Calloc(d, double);
    Sy = Calloc(n, double);
    Sy2 = Calloc(n, double);
    fm_Sm(p, m, BSy, BSy2, n, Sy, Sy2);
    //Rprintf("NT: gama=%f\n", gama[0]);
    logblik_ph_derv(gama, d, x, x0, n0, n1, Sy, Sy2, ell, dell, ddell);
    for(i=0;i<d;i++) del+=fabs(dell[i]); 
    while(it<maxit && del>eps){
    //Rprintf("NT: m=%d, it=%d, del=%e, gama=%f\n",m,  it, del, gama[0]);
        minverse(ddell, d);  
        for(i=0;i<d;i++){
            tmp[i] = 0.0;
            for(j=0;j<d;j++) tmp[i] += ddell[i+d*j]*dell[j];
        }
        del = 0.0;
        for(i=0;i<d;i++){
            gama[i] -= tmp[i];
            del += fabs(tmp[i]);
        }
        logblik_ph_derv(gama, d, x, x0, n0, n1, Sy, Sy2, ell, dell, ddell);
        for(i=0;i<d;i++) del+=fabs(dell[i]);
        it++;
        R_CheckUserInterrupt();
    }
    if(prog==0) Rprintf("NT: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del, ell[0]);
    Free(tmp); 
    Free(Sy); 
    Free(Sy2);
}
/*////////////////////////////////////////////////////*/
/*   Initializing p for fm(.|x1; p), x1=x0(gama1)     */
/*       using fm(.|x0; p0), where x0=x0(gama0)       */
/*            dgx0 = gama1*(x1-x0)                    */
/*////////////////////////////////////////////////////*/
void initialize_p(double *p, int m, double dgx0){
    int i, j, mp1=m+1;
    double pi0=0.0, sum_p=0.0, edgx0, *tmp, *Tmp;
    tmp=Calloc(mp1, double);
    Tmp=Calloc(mp1, double);
    edgx0=exp(dgx0);
    pi0=1.0-R_pow(p[mp1], edgx0);
    //Rprintf("Init: pi0=%f\n",pi0);
    for(i=0; i<=m;i++){ 
        tmp[i]=0.0;
        Tmp[i]=0.0;
        for(j=0;j<mp1;j++){
            tmp[i]+=p[j]*dbeta(i/(double) m, j+1,m-j+1, FALSE);
            Tmp[i]+=p[j]*(1-pbeta(i/(double) m, j+1,m-j+1, TRUE, FALSE)); 
        }
//        Tmp[i]+=p[mp1];
    }
    for(i=0; i<=m;i++){ 
        p[i]=edgx0*R_pow(Tmp[i], edgx0-1.0)*tmp[i];
        sum_p+=p[i];
    }
//    Rprintf("Init: sum_p=%f\n",sum_p);
    for(i=0; i<=m; i++){
    p[i]=pi0*p[i]/sum_p;
    //Rprintf("Init: p[i]=%f\n",p[i]);
    }
    p[mp1]=1-pi0;
    Free(tmp); 
    Free(Tmp); 
}

/*////////////////////////////////////////////////////*/
/* maximizer p of ell(gamma, p) for a gvien gamma     */
/*  egx: exp(gama*x.tilde), where gama is the given   */
/*       regression coefficient                       */
/*  ss: step-size epsilon:  p = (1-ss)p+ss Psi(p),    */
/*        default ss=1 so that p = Psi(p)             */
/*////////////////////////////////////////////////////*/
void pofg_ph(double *p, int m, double *egx, int n0, int n1, double *BSy, double *BSy2, 
        double *llik, double eps, int maxit, int prog, int *conv, double *delta){
    int i, j, n=n0+n1, mp1=m+1, mp2=m+2, it=0;
    double sum_egx=0.0,  del=1.0, Sp, Sp2;
    double *Tmp, *Tmp2, *pnu, tmp1, tmp2, llik_nu;
    Tmp=Calloc(mp2, double);
    Tmp2=Calloc(mp2, double);
    pnu=Calloc(mp2, double);
    for(i=0;i<n;i++) sum_egx+=egx[i];
    llik[0]=log_blik_ph(p, m, egx, n0, n1, BSy, BSy2); 
    while(del>eps && it<maxit){
        for(j=0;j<mp2;j++) pnu[j]=0.0;
        // p = p *w(p)/sum(exp(gamma'x))
        for(i=0; i<n0;i++){
            Sp=0.0; Sp2=0.0; 
            for(j=0;j<mp2;j++){
                Tmp[j]=BSy[i+n*j]*p[j];
                Sp+=Tmp[j];
                Tmp2[j]=BSy2[i+n*j]*p[j];
                Sp2+=Tmp2[j];
            }
            for(j=0;j<=m;j++){
                pnu[j]+=Tmp2[j]/Sp2;
                pnu[j]+=Tmp[j]*(egx[i]-1.0)/Sp;
            }
            pnu[mp1]+=Tmp[mp1]*(egx[i]-1.0)/Sp;
        }
        for(i=n0; i<n;i++){
            Sp=0.0; Sp2=0.0;
            for(j=0;j<mp2;j++){
                Tmp[j]=BSy[i+n*j]*p[j];
                Sp+=Tmp[j];
                Tmp2[j]=BSy2[i+n*j]*p[j];
                Sp2+=Tmp2[j];
            }
            tmp1=R_pow(Sp,egx[i]);
            tmp2=R_pow(Sp2,egx[i]);
            for(j=0;j<mp2;j++){
                pnu[j]+=(Tmp[j]*R_pow(Sp,egx[i]-1.0)-Tmp2[j]*R_pow(Sp2,egx[i]-1.0))*egx[i]/(tmp1-tmp2);
            }
        }  
        for(j=0;j<mp2;j++) pnu[j] /= sum_egx;
        //pnu<-(1-ss)*p+ss*pnu
        llik_nu=log_blik_ph(pnu, m, egx, n0, n1, BSy, BSy2); 
        del=fabs(llik[0]-llik_nu);
        it++;  llik[0]=llik_nu;
        for(j=0;j<mp2;j++) p[j]=pnu[j];
        R_CheckUserInterrupt();
    }
    if(prog==0) Rprintf("EM: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del, llik[0]);
    conv[0]=0;
    delta[0]=del;
    if(it==maxit){
        conv[0]+=1;
        //warning("\nThe maximum iterations were reached \nwith a delta = %f.\n", del);
    }
    Free(Tmp); 
    Free(Tmp2); 
    Free(pnu);//Free(Sp);  Free(Sp2);
}
/*////////////////////////////////////////////////////////////*/
/* Maximum Approximate Profile Likelihhod Estimation in       */
/*              PH model with a given gamma                   */
/* M: set of positive integers as candidate degrees of        */ 
/*       Bernstein poly model                                 */
/* gama: an efficient estimate of regression coefficient      */ 
/*       gamma, for data without covariate we set gama=0      */
/*    x: covariate centered at x0 satisfying                  */
/*             gama'x0=min{gama'xi, i=1,...,n}                */
/*////////////////////////////////////////////////////////////*/
void mable_ph_gamma(int *M, double *gama, int *dm, double *pi0, double *x, 
        double *y, double *y2, int *N, double *x0, double *lk, double *lr, 
        double *p, double *ddell, double *eps, int *maxit, int *progress,
        double *level, double *pval, int *chpts, int *conv, double *delta){
    int i, j, d=dm[0], k=M[1]-M[0], n0=N[0], n1=N[1], n=n0+n1;
    int *cp, lp, tmp, m=M[1],  mp1=m+1, mp2=m+2, itmp=0; 
    double tini=.000001, pct, ttl, *res;
    double *ell, *dell, *egx, *BSy, *BSy2, *phat, *Sy, *Sy2;//*xt, 
    lp=M[0]*(k+1)+(k+1)*(k+4)/2;
    phat = Calloc(lp, double);
    cp = Calloc(1, int);
    res = Calloc(1, double);
    ell = Calloc(1, double);
    dell = Calloc(d, double);
    BSy = Calloc(n*mp2, double);  
    BSy2 = Calloc(n*mp2, double); 
    Sy = Calloc(n, double);
    Sy2 = Calloc(n, double);
    egx = Calloc(n, double);
    if(*progress==1) {Rprintf("\n Mable fit of PH model with the give gamma ... \n");
        ProgressBar(0.0,""); }
    ttl = (double)((k+2)*(k+1));
//    egx_x0(gama, d, x, n, egx, x0);
    egxmx0(gama, d, x, n, egx, x0);
    // add check to see if any egx is less than 1
    for(i=0;i<n;i++) 
        if (egx[i]<1) {
            Rprintf("\n");
            error("Try another baseline 'x0'.\n");}
    m=M[0]; 
    mp1=m+1;
    mp2=m+2;
    Bdata(y, m, 0, n, BSy);
    Bdata(y2, m, n0, n1, BSy2);
    for(i=0;i<=m;i++) p[i]=*pi0/(double) mp1; 
    p[mp1]=1.0-*pi0;
    pofg_ph(p, m, egx, n0, n1, BSy, BSy2, ell, *eps, *maxit, *progress, conv, delta);
    itmp+=conv[0];
    tmp=mp2;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    lk[0]=ell[0]; 
    chpts[0]=0;
    pval[0]=1.0;
    pct = 2/ttl;
    if(*progress==1) ProgressBar(pct,"");
    i=1; 
    while(i<=k && pval[i-1]>*level){
        p[mp2]=p[mp1];
        p[mp1] = mp1*p[m]/(double)mp2;
        for(j=m; j>0; j--) p[j] = (p[j-1]*j+p[j]*(mp1-j))/(double)mp2;
        p[0] = mp1*p[0]/(double)mp2;
        m=M[0]+i;
        mp1=m+1;
        mp2=m+2;
        Bdata(y, m, 0, n, BSy);
        Bdata(y2, m, n0, n1, BSy2);
        for(j=0;j<=mp1;j++) p[j]=(p[j]+tini/(double) mp2)/(1.0+tini);
        pofg_ph(p, m, egx, n0, n1, BSy, BSy2, ell, *eps, *maxit, *progress, conv, res);
        lk[i]=ell[0];
        for(j=0;j<=mp1;j++) phat[j+tmp]=p[j];
        tmp += mp2;
//        Rprintf("lk[%d]=%f\n",i, lk[i]);
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        R_CheckUserInterrupt();
        pct += 2*(i+1)/ttl;
        if(*progress==1) ProgressBar(pct,""); 
        i++;
        itmp+=conv[0];
    }
    if(*progress==1){
        ProgressBar(1.0,"");
        Rprintf("\n");}
    if(itmp>0) conv[0]=1; 
    else conv[0]=0;
    if(k>0){
        if(m==M[1]){
            conv[0]+=1; 
            warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);
            delta[0]=res[0];}
    }
    M[1]=m;
    tmp=cp[0]*(M[0]*2+(cp[0]+3))/2;
    dm[1]=cp[0]+M[0];
    m=dm[1];
    for(j=0;j<=m+1;j++) p[j]=phat[tmp+j];
    Free(phat);  
    Free(cp);  
    Free(res);  
    Free(ell);  
    Free(dell);
    Free(BSy); 
    Free(BSy2); 
    Free(Sy); 
    Free(Sy2); 
    Free(egx);  
}

/*////////////////////////////////////////////////////*/
/*  Maximum approx. Bernstein likelihood estimate of  */
/*           (gamma, p) with a fixed degree m         */
/*////////////////////////////////////////////////////*/

void mable_ph_m(double *gama, double *p, int *dm, double *x, double *y, double *y2, 
        int *N, double *x0, double *ell, double *ddell, double *EPS, int *MAXIT,
        int *progress, int *conv, double *delta){
    int i, j, n0=N[0], n1=N[1], n=n0+n1, d=dm[0], m=dm[1], mp2=m+2, it=0;
    int maxit=MAXIT[0], maxit_em=MAXIT[1], maxit_nt=MAXIT[2], prog=1; 
    double eps=EPS[0], eps_em=EPS[1], eps_nt=EPS[2];
    double tini, *BSy, *BSy2, *gnu, del, pct=0.0; 
    double *ell1, *dell, *egx;
    tini=0.00001;// tini is used to make sure p is in interior of S(m+1)
    ell1 = Calloc(1, double);  
    dell = Calloc(d, double);  
    egx = Calloc(n, double);  
    BSy = Calloc(n*mp2, double);  
    BSy2 = Calloc(n*mp2, double);  
    gnu = Calloc(d, double);  
    Bdata(y, m, 0, n, BSy);
    Bdata(y2, m, n0, n1, BSy2);
//    egx_x0(gama, d, x, n, egx, x0);
    egxmx0(gama, d, x, n, egx, x0);
    // add check to see if any egx is less than 1
    for(i=0;i<n;i++) 
        if (egx[i]<1) {
            Rprintf("\n");
            error("Try another baseline 'x0'.\n");}
    for(j=0;j<d;j++)  gnu[j]=gama[j];
    //Rprintf("gnu=%f\n", gnu[0]);
    pofg_ph(p, m, egx, n0, n1, BSy, BSy2, ell, eps_em, maxit_em, prog, conv, delta);
    gofp_ph(gnu, d, p, m, x, x0, n0, n1, BSy, BSy2, ell, dell, ddell, eps_nt, maxit_nt, prog);
    del=0.0;
    for(i=0;i<d;i++){
        del+=fabs(gnu[i]-gama[i]);
        gama[i]=gnu[i];
    }
    if(*progress==1) ProgressBar(pct,""); 
    while(it<maxit && del>eps){
        //egx_x0(gama, d, x, n, egx, x0);
        egxmx0(gama, d, x, n, egx, x0);
        // add check to see if any egx is less than 1
        for(i=0;i<n;i++) 
            if (egx[i]<1) {
                Rprintf("\n");
                error("Try another baseline 'x0'.\n");}
        for(i=0;i<mp2; i++) p[i]=(p[i]+tini/(double) mp2)/(1.0+tini); 
        pofg_ph(p, m, egx, n0, n1, BSy, BSy2, ell1, eps_em, maxit_em, prog, conv, delta);
        gofp_ph(gnu, d, p, m, x, x0, n0, n1, BSy, BSy2, ell1, dell, ddell, eps_nt, maxit_nt, prog);
        del=0.0;
        for(i=0;i<d;i++){
            del+=fabs(gnu[i]-gama[i]);
            gama[i]=gnu[i];
        }
        del+=fabs(ell1[0]-ell[0]);
        ell[0]=ell1[0];
        it++;
        pct = fmin(it/(double) maxit, (1+eps)/(1+del));
        if(*progress==1) ProgressBar(pct,""); 
        R_CheckUserInterrupt();
//        Rprintf("         mable-m: it=%d, del=%f\n", it, del);
    }
    if(*progress==1) {
        ProgressBar(1.0,"");
        Rprintf("\n");
    }
    conv[0]=0;
    delta[0]=del;
    if(it==maxit){
        conv[0]+=1;
        //warning("\nThe maximum iterations were reached \nwith a delta = %f.\n", del);
    }
//    Rprintf("mable-m: it=%d, del=%f\n", it, del);
    minverse(ddell, d); //Sig=-n*ddell
    Free(BSy); 
    Free(BSy2); 
    Free(gnu); 
    Free(ell1); 
    Free(dell); 
    Free(egx); 
}

/*////////////////////////////////////////////////////////*/
/*  MABLE of (gamma, p) and an optimal degree m           */
/*                                                        */
/*    M: set of positive integers as candidate degrees    */
/*       of Bernstein poly model                          */
/*   gama: initial value of gamma, an efficient       */
/*     estimate  of regression coefficient gamma          */
/*   phat: (m1+2)-vector, first mt+2 are mable of p       */
/*         obtained with gamma=gama_hat                   */
/*     x: covariate                                       */
/*    x0: gama'x0=min{gama'xi, i=1,...,n}                 */
/*////////////////////////////////////////////////////////*/
void mable_ph(int *M, double *gama, int *dm, double *p, double *pi0, double *x,   
        double *y, double *y2, int *N, double *x0, double *lk, double *lr, 
        double *ddell, double *EPS, int *MAXIT, int *progress, double *level,
        double *pval, int *chpts, int *conv){
    int d=dm[0], i, j, k=M[1]-M[0], prg=1-*progress;
    int m, *cp, tmp, lp;  
    double *ghat, *phat, *res, *ell, pct, ttl; 
    lp=M[0]*(k+1)+(k+1)*(k+4)/2;
    cp = Calloc(1, int);
    res = Calloc(1, double);
    phat=Calloc(lp, double);
    ghat=Calloc(d*(k+1), double);
    ell=Calloc(1, double);
//    egx=Calloc(n, double);
    if(*progress==1) {Rprintf("\n Mable fit of Cox PH regression model ... \n");
        ProgressBar(0.0,""); }
    ttl=(double)(k+2)*(k+1);
    m=M[0]; 
    dm[1]=m;
    for(i=0;i<=m;i++) p[i]=*pi0/(double)(m+1);
    p[m+1] = 1-*pi0;
    mable_ph_m(gama, p, dm, x, y, y2, N, x0, ell, ddell, EPS, MAXIT, &prg, conv, res);
    for(i=0;i<dm[0];i++) ghat[i]=gama[i];
    tmp=m+2;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    lk[0]=ell[0];
    pval[0]=1.0;
    chpts[0]=0;
    pct = 2/ttl;
    if(*progress==1) ProgressBar(pct,""); 
    i=1;
    while(i<=k && pval[i-1]>*level){
        p[m+2] = p[m+1];
        p[m+1] = (m+1)*p[m]/(double)(m+2);
        for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m+1-j))/(double)(m+2);
        p[0] = (m+1)*p[0]/(double)(m+2);
        m=M[0]+i;
        dm[1]=m;
        for(j=0;j<=m+1;j++) p[j]=(p[j]+.000001/(double)(m+2))/1.000001;
        mable_ph_m(gama, p, dm, x, y, y2, N, x0,  ell, ddell, EPS, MAXIT, &prg, conv, res);
        lk[i]=ell[0];  
        for(j=0;j<=m+1;j++) phat[j+tmp]=p[j];
        tmp += m+2;
        for(j=0; j<dm[0]; j++) ghat[dm[0]*i+j]=gama[j];
//        Rprintf("lk[%d]=%f\n",i, lk[i]);
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        R_CheckUserInterrupt();
        pct +=2*(i+1)/ttl;
        if(*progress==1) ProgressBar(pct,""); 
        i++;
    }
    if(*progress==1){
        ProgressBar(1.0,"");
        Rprintf("\n");}
    if(m==M[1]){
        conv[0]+=1; 
        warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);}
//    else conv[0]=0;
    M[1]=m;
    tmp=cp[0]*(M[0]*2+(cp[0]+3))/2;
    dm[1]=cp[0]+M[0];
    m=dm[1];
    for(j=0;j<=m+1;j++) p[j]=phat[tmp+j];
    for(j=0; j<dm[0]; j++) gama[j]=ghat[dm[0]*cp[0]+j];
    if(*progress==1) Rprintf("\n");
    Free(phat);
    Free(ghat);
    Free(ell); 
//    Free(egx);
    Free(cp);
    Free(res);
}
/*///////////////////////////////////////////*/
/*       Without Covariate                   */
/*///////////////////////////////////////////*/
/* mable of p for interval censored data no covariate */
/*  ss: step-size epsilon:  p = (1-ss)p+ss Psi(p),    */
/*        default ss=1 so that p = Psi(p)             */
/*////////////////////////////////////////////////////*/
void mablem_ic(double *p, int m, int n0, int n1, double *egx, double *BSy, double *BSy2, 
        double *llik, double eps, int maxit, int prog, int *conv, double *delta){
    int i, j, n=n0+n1, mp2=m+2, it=0;
    double del=1.0, Sp, Sp2;//*egx, ,  llk=llik[0]
    double *Tmp, *Tmp2, *pnu, llik_nu;
//    egx=Calloc(n, double);
    Tmp=Calloc(n*mp2, double);
    Tmp2=Calloc(n*mp2, double);
    pnu=Calloc(mp2, double);
//    for(i=0;i<n;i++)  egx[i]=1.0;
    llik[0]=log_blik_ph(p, m, egx, n0, n1, BSy, BSy2); 
    while(del>eps && it<maxit){// || (llik[0]<=llk)){
        for(j=0;j<mp2;j++) pnu[j]=0.0;
        // p = p *w(p)/n
        for(i=0; i<n0;i++){
            Sp=0.0; Sp2=0.0; 
            for(j=0;j<=m;j++){
                Tmp[i+n*j]=BSy[i+n*j]*p[j];//unused
                Sp+=Tmp[i+n*j];//unused
                Tmp2[i+n*j]=BSy2[i+n*j]*p[j];
                Sp2+=Tmp2[i+n*j];
            }
            for(j=0;j<=m;j++){
                pnu[j]+=Tmp2[i+n*j]/Sp2;
            }
        }
        for(i=n0; i<n;i++){
            Sp=0.0; Sp2=0.0;
            for(j=0;j<mp2;j++){
                Tmp[i+n*j]=BSy[i+n*j]*p[j];
                Sp+=Tmp[i+n*j];
                Tmp2[i+n*j]=BSy2[i+n*j]*p[j];
                Sp2+=Tmp2[i+n*j];
            }
            for(j=0;j<mp2;j++){
                pnu[j]+=(Tmp[i+n*j]-Tmp2[i+n*j])/(Sp-Sp2);
            }
        }  
        for(j=0;j<mp2;j++) pnu[j] /= (double)n;
        //pnu<-(1-ss)*p+ss*pnu
        llik_nu=log_blik_ph(pnu, m, egx, n0, n1, BSy, BSy2); 
        del=fabs(llik[0]-llik_nu);
        it++;  llik[0]=llik_nu;
        for(j=0;j<mp2;j++) p[j]=pnu[j];
        R_CheckUserInterrupt();
    }
    if(prog==0) Rprintf("EM: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del, llik[0]);
    conv[0]=0;
    delta[0]=del;
    if(it==maxit) conv[0]+=1;
    Free(Tmp);  Free(Tmp2);  Free(pnu);//Free(Sp);  Free(Sp2); Free(egx); 
}
/*////////////////////////////////////////////////////////////*/
/*  Find optimal degree m and phat for interval censored      */
/*  data without covariate                                    */
/* M: set of positive integers as candidate degrees of        */ 
/*       Bernstein poly model                                 */
/*////////////////////////////////////////////////////////////*/
void mable_ic(int *M, double *pi0, double *y, double *y2, int *N, double *lk,
        double *lr, double *p, double *eps, int *maxit, int *progress, double *pval, 
        double *bic, int *chpts, int *optim, double *level, int *conv, double *delta){
    int d, i,j, n0=N[0], n1=N[1], n=n0+n1, k=M[1]-M[0], itmp=0;
    int lp, m=M[1],  mp1=m+1, mp2=m+2, tmp, *cp;//
    double tini=.000001, pct=0.0, ttl, *egx; //max_bic, 
    double *ell, *BSy, *BSy2, *phat, *Sy, *Sy2, *res; 
    lp=M[0]*(k+1)+(k+1)*(k+4)/2;
    cp = Calloc(1, int);
    res = Calloc(1, double);
    egx=Calloc(n, double);
    phat=Calloc(lp, double);
    ell = Calloc(1, double);
    BSy = Calloc(n*mp2, double);  
    BSy2 = Calloc(n*mp2, double); 
    Sy = Calloc(n, double);
    Sy2 = Calloc(n, double);
    if(*progress==1) {Rprintf("\n Fitting interval censored data ... \n");
        ProgressBar(0.0,""); }
    ttl = (double)((k+2)*(k+1));
    for(i=0;i<n;i++)  egx[i]=1.0;
    m=M[0]; 
    mp1=m+1;
    mp2=m+2;
    Bdata(y, m, 0, n, BSy);
    Bdata(y2, m, n0, n1, BSy2);
    for(i=0;i<=m;i++) p[i] = *pi0/(double)mp1;
    p[mp1] = 1.0-*pi0;
//    ell[0]=log_blik_ph(p, m, egx, n0, n1, BSy, BSy2); 
    mablem_ic(p, m, n0, n1, egx, BSy, BSy2, ell, eps[0], *maxit, *progress, conv, delta);
    itmp+=conv[0];
    tmp=mp2;//i=0; tmp=(i+1)*M[0]+(i+1)*(i+4)/2-1;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    lk[0]=ell[0]; 
    pval[0]=1.0;
    d=0;
    for(j=0;j<=mp1;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;
    bic[0]=lk[0]-.5*d*log(n);
//    max_bic=bic[0];
//    optim[1]=0;
    pct += 2/ttl;
    if(*progress==1) ProgressBar(pct,""); 
    i=1;
    while(i<=k && pval[i-1]>*level){
        p[mp2]=p[mp1];
        p[mp1] = mp1*p[m]/(double)mp2;
        for(j=m; j>0; j--) p[j] = (p[j-1]*j+p[j]*(mp1-j))/(double)mp2;
        p[0] = mp1*p[0]/(double)mp2;
        m=M[0]+i;
        mp1=m+1;
        mp2=m+2;
        for(j=0;j<=mp1;j++) p[j]=(p[j]+tini/(double) mp2)/(1.0+tini);
        Bdata(y, m, 0, n, BSy);
        Bdata(y2, m, n0, n1, BSy2);
        mablem_ic(p, m, n0, n1, egx, BSy, BSy2, ell, eps[0], *maxit, *progress, conv, delta);
        lk[i]=ell[0];  
        for(j=0;j<=mp1;j++) phat[j+tmp]=p[j];
        tmp += mp2;
//        Rprintf("lk[%d]=%f\n",i, lk[i]);
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        // Calculate BIC
        d=0;
        for(j=0;j<=mp1;j++) d+=1*(p[j]>=eps[1]);
        d=d-1;// bic not working
        bic[i]=lk[i]-.5*d*log(n);
//        if(bic[i]>max_bic){
//            max_bic=bic[i];
//            optim[1]=i;
//        }
        R_CheckUserInterrupt();
        pct +=2*(i+1)/ttl;
        if(*progress==1) ProgressBar(pct,""); 
        i++; 
        itmp+=conv[0];
    }
    if(*progress==1){
        ProgressBar(1.0,"");
        Rprintf("\n");}
    M[1]=m;
    tmp=cp[0]*(M[0]*2+(cp[0]+3))/2;
    optim[0]=cp[0]+M[0];
    m=optim[0];
    for(j=0;j<=m+1;j++) p[j]=phat[tmp+j];
    if(itmp>0) conv[0]=1; 
    else conv[0]=0;
    if(k>0){
        if(m==M[1]){
            conv[0]+=1; 
            warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);
            delta[0]=res[0];
        }
//        tmp=optim[1]*(M[0]*2+(optim[1]+3))/2;
//        optim[1]+=M[0];
//        m=optim[1];
//        for(j=0;j<=m+1;j++) p[optim[0]+2+j]=phat[tmp+j];
    }
    Free(cp);  
    Free(res);  
    Free(egx);  
    Free(phat);  
    Free(ell);  
    Free(BSy); 
    Free(BSy2); 
    Free(Sy); 
    Free(Sy2); 
}
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*      Maximum Bernstein likelihood Deconvolution        */
/*                     EM-Algorithm                       */
/*                                                        */
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/*                                                            */
/*  Reference:                                                */
/*   Zhong Guan,  Fast Nonparametric Maximum Likelihood       */
/*     Density Deconvolution Using Bernstein Polynomial       */
/*                                                            */
/*////////////////////////////////////////////////////////////*/
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
/*//////////////////////////////////////////////////////////////////*/
/*     Approximate Bernstein Log-Likelihood of y=x+e                */
/*  x: n-vector, sample from distribution F with support [0,1]      */
/*  e: n-vector, errors with known density g                        */
/*  p=(p0, p1, ..., pm),                                            */
/*  gBeta: (m+1)x n matrix of convolution of g and Beta(i+1, m-i+1) */
/*//////////////////////////////////////////////////////////////////*/
double approx_bern_lik(double *p, double *gBeta, int m, int n){
    int i,j;
    double llik, fx;
    llik = 0.0;
    for(i=0; i<n; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*gBeta[i+n*j];
//            Rprintf("  p = %g\n", p[j]);
        }
//        llik *= fx;
        llik += log(fx);
    }
//    Rprintf("  lik = %g\n", llik);
//    return(llik+n*log(m+1.));
    return(llik);
}
/*////////////////////////////////////////////////////////////*/
/*  mable density deconvolution with general error density    */
/*      using .External, code modified from 'integrate.c'     */
/*////////////////////////////////////////////////////////////*/
/* called via .External(.) :*/
SEXP C_mable_decon(SEXP args);

typedef struct int_struct
{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
} int_struct, *IntStruct;
typedef struct mable_struct
{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
    int m;
    int j;
    double y;
} mable_struct, *MableStruct;

/*//////////////////////////////////////////////////////////////////////*/
/*          The integrand eta(x,y;j,m)=g(y-x)*beta_mj(x) of the           */
/*       convolution of beta_{mj} and the error density g(x,...)        */
/*      returns a function of x for fixed y, ex = (env, m, j, y)        */
/*           code modified from 'Rintfn()' in 'integrate.c'             */
/*//////////////////////////////////////////////////////////////////////*/
static void eta_mj(double *x, int n, void *ex)
{
    SEXP args, resultsxp, tmp;
    int i, j, m;
    double y, *z;
    z = Calloc(n, double);
    MableStruct MS = (MableStruct) ex;
    m = MS->m;
    j = MS->j;
    y = MS->y;
    PROTECT(args = allocVector(REALSXP, n));
    for(i = 0; i < n; i++) REAL(args)[i] = y - x[i];

    PROTECT(tmp = lang2(MS->f, args));
    PROTECT(resultsxp = eval(tmp, MS->env));

    if(length(resultsxp) != n)
	error("evaluation of function gave a result of wrong length");
    if(TYPEOF(resultsxp) == INTSXP) {
	   resultsxp = coerceVector(resultsxp, REALSXP);
    } 
    else if(TYPEOF(resultsxp) != REALSXP)
	   error("evaluation of error density gave a result of wrong type");
    for(i = 0; i < n; i++){
	   z[i] = REAL(resultsxp)[i];
	   x[i] = z[i]*dbeta(x[i], j+1,  m-j+1, FALSE);
	   if(!R_FINITE(x[i]))
	       error("non-finite error density value");
    }
    UNPROTECT(3);
    Free(z);
    return;
}
/*////////////////////////////////////////////////////////////////////////*/
/*         Calculate the convolution psi_m(y)=beta_{mj}*g(y),             */  
/*                          Input: (m, y)                                 */
/*////////////////////////////////////////////////////////////////////////*/

void convol_beta_g(double *y, double *psi_m, int m, int n, void *ex){
    int i, j;
    mable_struct ms;
    IntStruct IS = (IntStruct) ex;
    double l=.0, u=1.0, epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400]; 
    int lenw=400, last=0, neval=0, ier=0, iwork[100]; 
    int limit=100; 
    ms.f = IS->f;
    ms.env = IS->env;
    ms.m = m;
    for(i=0;i<n;i++){
        ms.y = y[i];
        for(j=0; j<=m; j++){
            ms.j = j;
            Rdqags(eta_mj, (void*)&ms, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
            psi_m[i+j*n] = result;
        }
    }
}
/*//////////////////////////////////////////////////*/
/*          EM Algorithm for a fixed m              */
/*//////////////////////////////////////////////////*/
/* EM Method for mixture of g*beta(i+1, m+1-i), i=0,...,m. n=sample size,  */
void em_gBeta_mix(double *y, double *p, int m, int n, int maxit, double eps, double *llik, void *ex){
    int i, j, it;
    double *gBeta, del, llik_nu, *p_gBeta,  *fp; 
    gBeta = Calloc((m+1)*n, double);
    p_gBeta = Calloc((m+1)*n, double);
    fp = Calloc(n, double);
    convol_beta_g(y, gBeta, m, n, ex);
    llik[0] = 0.0;
    for(i=0; i<n; i++){
        fp[i] = 0.0;
        for(j=0; j<=m; j++) {
            p_gBeta[i+n*j] = p[j]*gBeta[i+n*j];
            fp[i] += p_gBeta[i+n*j];
        }
        llik[0]+=log(fp[i]);
    }
    del = 10.0;
    it = 1;
    while(del>eps && it<maxit){
        for(j=0; j<=m; j++){
            p[j] = 0.0;
            for(i=0; i<n; i++) p[j] += p_gBeta[i+n*j]/fp[i];
            p[j] /= (double)n;
        }
        llik_nu=0.0;
        for(i=0; i<n; i++){
            fp[i] = 0.0;
            for(j=0; j<=m; j++) {
                p_gBeta[i+n*j] = p[j]*gBeta[i+n*j];
                fp[i] += p_gBeta[i+n*j];
            }
            llik_nu+=log(fp[i]);
        }
        del = fabs(llik[0]-llik_nu);
        it += 1;
        llik[0] = llik_nu;
//        Rprintf("Iteration: %d, Del = %g\n", it-1, del);
    }
//    Rprintf("Number of EM Iterations = %d, Del = %g\n", it-1, del);
    Free(gBeta);
    Free(p_gBeta);
    Free(fp);
}
/* end function em_gBeta_mix */

/*////////////////////////////////////////////////////////////////*/
/*     Maximum Approximate Bernstein Likelihood Estimation        */
/*  for Density Deconvolution with an additive measurement error  */
/*                 and a known error distribution g               */
/*////////////////////////////////////////////////////////////////*/
// p=(1,...,1)/(mu+1)
SEXP mable_decon(SEXP args)
{
    int_struct is;
    int n, maxit_em, progress;
    SEXP ans, ansnames, yy;
    int d, i, j, m, *M,  nm, *chpts, lp, tmp, *cp; //*optim,
    double *y, *p, *phat, eps_em, pct, *res;//max_bic=0, 
    double *lr, *lk,  *llik, ttl, *pval, *bic, level, eps;  
    Rprintf("\n Deconvolution is runing. This may take several minutes.\n\n");
    args = CDR(args);
    is.f = CAR(args); args = CDR(args);
    is.env = CAR(args); args = CDR(args);
    yy = CAR(args); args = CDR(args);
//    p = REAL(CAR(args)); args = CDR(args);
    M = INTEGER(CAR(args)); args = CDR(args);
//    mu = asInteger(CAR(args)); args = CDR(args);
    eps_em = asReal(CAR(args)); args = CDR(args);
    maxit_em = asInteger(CAR(args)); args = CDR(args);
    progress = asInteger(CAR(args)); args = CDR(args);
    level = asReal(CAR(args)); args = CDR(args);
    eps = asReal(CAR(args)); args = CDR(args);
    n = length(yy);
    y= REAL(yy);
    nm =M[1]-M[0];//  
    ttl = (double)(nm*(nm-1));
//    Rprintf("  dim:  %d\n", (k+1)*(*d+1));
    llik = Calloc(1, double);
    lk = Calloc(nm+1, double);
    lr = Calloc(nm, double);
    pval = Calloc(nm+1, double);
    bic = Calloc(nm+1, double);
    chpts = Calloc(nm+1, int); 
    cp = Calloc(1, int); 
    res = Calloc(1, double); 
    lp = M[0]*(nm+1)+(nm+1)*(nm+2)/2;
    p = Calloc(M[1]+1, double); // Calloc(2*M[1]+2, double);
    phat = Calloc(lp, double);  
    m = M[0]; 
    for(j=0; j<=m; j++) p[j]=1.0/(double) (m+1);
    em_gBeta_mix(y, p, m, n, maxit_em, eps_em, llik, (void*)&is);
    for(j=0; j<=m; j++) phat[j]=p[j];
    tmp = m+1;
    lk[0] = *llik;
    pval[0]=1.0;
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps);
    d=d-1;
    bic[0]=lk[0]-.5*d*log(n);
    i=1;
    while(i<=nm && pval[i-1]>level){
        m = M[0]+i;
        for(j=0; j<=m; j++) p[j]=1.0/(double)(m+1.0);
        em_gBeta_mix(y, p, m, n, maxit_em, eps_em, llik, (void*)&is);
        lk[i] = *llik;
        for(j=0;j<=m;j++) phat[j+tmp]=p[j];
        tmp += m+1;
        //   lr: exponential LR of change-point
        if(i>=3){
            cp[0]=i;
            chpt_exp(lk, lr, res, cp);
            pval[i]=res[0];
            chpts[i]=cp[0];
        }
        else{            
            pval[i]=1.0;
            chpts[i]=i;
        }
        // Calculate BIC
        d=0;
        for(j=0;j<=m;j++) d+=1*(p[j]>=eps);
        d=d-1;
        bic[i]=lk[i]-.5*d*log(n);
        pct = i*(i+1)/ttl;
        if(progress==1) {
            ProgressBar(pct,"");}
//      clockProgress(pct);
        i++;
    }
    if(progress==1){
        ProgressBar(1.0,"");
        Rprintf("\n");}
    if(m==M[1]){
        Rprintf("Warning: The maximum candidate degree has been reached.\n");}
    M[1]=m;
    nm = M[1]-M[0];
    Free(llik); 
    if(progress==1) Rprintf("\n");
    PROTECT(ans = allocVector(VECSXP, 8));
    PROTECT(ansnames = allocVector(STRSXP, 8));
    SET_STRING_ELT(ansnames, 0, mkChar("lk"));
    SET_STRING_ELT(ansnames, 1, mkChar("lr"));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, nm+1));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, nm));
    for(i=0;i<nm;i++){
        REAL(VECTOR_ELT(ans, 0))[i] = lk[i];
        REAL(VECTOR_ELT(ans, 1))[i] = lr[i];
    }
    REAL(VECTOR_ELT(ans, 0))[nm] = lk[nm];
    SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, cp[0]+M[0]+1));
    SET_STRING_ELT(ansnames, 2, mkChar("p"));
    tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
    m=cp[0]+M[0];
    for(i=0;i<=m;i++){
        REAL(VECTOR_ELT(ans, 2))[i] = phat[tmp+i];
    }
    SET_STRING_ELT(ansnames, 3, mkChar("m"));
    SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(ans, 3))[0] = m;
    SET_STRING_ELT(ansnames, 4, mkChar("pval"));
    SET_STRING_ELT(ansnames, 5, mkChar("bic"));
    SET_STRING_ELT(ansnames, 6, mkChar("chpts"));
    SET_VECTOR_ELT(ans, 4, allocVector(REALSXP, nm+1));
    SET_VECTOR_ELT(ans, 5, allocVector(REALSXP, nm+1));
    SET_VECTOR_ELT(ans, 6, allocVector(INTSXP, nm+1));
    for(i=0;i<=nm;i++){
        REAL(VECTOR_ELT(ans, 4))[i] = pval[i];
        REAL(VECTOR_ELT(ans, 5))[i] = bic[i];
        INTEGER(VECTOR_ELT(ans, 6))[i] = chpts[i];
    }
    SET_STRING_ELT(ansnames, 7, mkChar("M"));
    SET_VECTOR_ELT(ans, 7, allocVector(INTSXP, 2));
    for(i=0;i<=1;i++){
        INTEGER(VECTOR_ELT(ans, 7))[i] = M[i];
    }
    setAttrib(ans, R_NamesSymbol, ansnames);
    UNPROTECT(2);
    Free(lk); Free(lr); Free(phat); Free(p);
    Free(pval); Free(bic); Free(chpts); Free(cp); Free(res); 
    return ans;
}
// End of file mable-decon.c
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*    Maximum Approximate Bernstein likelihood Estimate   */
/*            Multivariate Density Functions              */
/*                                                        */
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/*                                                            */
/*  Reference:                                                */
/*  Wang, T. and Guan, Z.,(2019) Bernstein Polynomial Model   */
/*     for Nonparametric Multivariate Density, Statistics,    */
/*              Vol. 53, no. 2, 321-338                       */
/*////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/* Multivariate Bernstein base polynomials */
/* Multivar_Beta: Returns a K x n matrix, K=(m1+1)...(md+1) */
/* km[1]=m1+1,km[2]=(m1+1)(m2+1),...,km[d]=(m1+1)...(md+1)=K,*/
/* beta(x1j, i1+1, m1+1-i1)...beta(xdj, id+1, md+1-id), 0<=ik<=mk, 1<=k<=d,  
         the (i1+km[1]*i2+...+km[d-1]*id,j)th element */
void Multivar_dBeta(double *x, int *m, int n, int d, int *km, double *dBta) {
    int i, j, jj, k, r, K, it;
    K=km[d-1];
    for(it=0; it<K; it++)
        for(j=0; j<n; j++) dBta[it+K*j]=1.0; //get initial value of dBta of K*n metrix
    for(j=0; j<n; j++){
        for(k=0; k<d; k++){ //d-var desity function
            it=0;
            for(jj =1; jj*(m[k]+1)<=km[k]; jj++){
               for(i=0; i<=m[k]; i++){
                    for(r=1; r*km[k]<=K; r++){
                        dBta[it+K*j]*=dbeta(x[j+n*k], i+1, m[k]+1-i, FALSE);
//   Rprintf("it=%d, k=%d, i= %d, j= %d, x[%d,%d]=%f\n", it, k, i, j, j, k, x[j+n*k]);
                        it+=1;
                    }
                }
            }
        }
    }
}
/*//////////////////////////////////////////////////////////////////*/
/* Multivariate Beta cdfs                                           */
/* Multivar_Beta_CDF: Returns a K x n matrix, K=(m1+1)...(md+1)     */
/* km[1]=m1+1,km[2]=(m1+1)(m2+1),...,km[d]=(m1+1)...(md+1)=K,       */
/* pbeta(x1j, i1+1, m1+1-i1)...pbeta(xdj, id+1, md+1-id),           */
/* 0<=ik<=mk, 1<=k<=d,  the (i1+km[1]*i2+...+km[d-1]*id,j)th element*/
/*//////////////////////////////////////////////////////////////////*/
void Multivar_pBeta(double *x, int *m, int n, int d, int *km, 
                double *pBta) {
    int i, j, jj, k, r, K, it;
    K=km[d-1];
    for(it=0; it<K; it++)
        for(j=0; j<n; j++) pBta[it+K*j]=1.0;
    for(j=0; j<n; j++){
        for(k=0; k<d; k++){
            it=0;
            for(jj =1; jj*(m[k]+1)<=km[k]; jj++){
               for(i=0; i<=m[k]; i++){
                    for(r=1; r*km[k]<=K; r++){
                        pBta[it+K*j]*=pbeta(x[j+n*k], i+1, m[k]-i+1, TRUE, FALSE);
//Rprintf("it=%d, k=%d, i= %d, j= %d, x[%d,%d]=%f\n", it, k, i, j, j, k, x[j+n*k]);
                        it+=1;
                    }
                }
            }
        }
    }
}

/*//////////////////////////////////////////////////////*/
/*   Log-Bernstein-Likelihood of multivariate sample    */
/*  p=p(i1,...,id), 0<=ij<=mj, j=1,...,d,               */
/*         an array of dim=m+1=(m1+1,...,md+1),         */
/*  x: n x d matrix, sample from d-dimensional          */
/*        distribution F with support [0,1]^d           */
/*//////////////////////////////////////////////////////*/
double loglik_bern_multivar(double *p, int K, int n, double *Bta){
    int i,j;
    double llik, fx;
    llik = 1.0;
    for(j=0; j<n; j++){
        fx = 0.0;
        for(i=0; i<K; i++) fx += p[i]*Bta[i+K*j];
        llik *= fx;
    }
//   Rprintf("  lik = %g\n", llik);
    return(log(llik));
}

/*/////////////////////////////////////////////////////*/
/* EM Method for mixture of                            */
/* beta(x1, i1+1, m1+1-i1),...,beta(xd, id+1, md+1-id),*/ 
/*               0<=ik<=mk, 1<=k<=d                    */
/*/////////////////////////////////////////////////////*/
void em_multivar_beta_mix(double *p, double *x, int *m, 
        int n, int d, int *km, int maxit, double eps, 
        double *llik, int progress, int *conv){
    int i, j, it, K;
    double del, llik_nu, *pBeta, *Bta, *fp, *pnu;
    double  ttl;
    ttl=(double) maxit;
    conv[0] = 0;
    K=km[d-1];
//    Rprintf("  d = %d\n", d);
//    for(i=0; i<d;i++) Rprintf("  km[%d] = %d\n", i, km[i]);
//    Rprintf("  K = %d\n", K);
//    for(i=0; i<d;i++) for(j=0; j<n;j++) Rprintf("x[%d,%d]=%f\n",  j, i, x[j+n*i]);
    Bta = Calloc(K*n, double);
    pBeta = Calloc(K*n, double);
    fp = Calloc(n, double);
    pnu = Calloc(K, double);
    Multivar_dBeta(x, m, n, d, km, Bta);
    llik[0] = loglik_bern_multivar(p, K, n, Bta);
    del = 10.0;
    it = 1;
    while(del>eps && it<maxit){
        for(j=0; j<n; j++){
            fp[j] = 0.0;
            for(i=0; i<K; i++) {
                pBeta[i+K*j] = p[i]*Bta[i+K*j];
                fp[j] += pBeta[i+K*j];
            }
        }
        for(i=0; i<K; i++){
            pnu[i] = 0.0;
            for(j=0; j<n; j++) pnu[i] += pBeta[i+K*j]/fp[j];
            pnu[i] /= (double)n;
        }
        llik_nu = loglik_bern_multivar(pnu, K, n, Bta);
        del = fabs(llik[0]-llik_nu);
        it += 1;
        for(i=0; i<K; i++) p[i] = pnu[i];
        llik[0] = llik_nu;
        R_CheckUserInterrupt();
//        Rprintf("  llik = %g\n", llik[0]);
        if(progress==1) ProgressBar(it/ttl,"");
    }
    if(progress==1){
        ProgressBar(1.0,"");
        Rprintf("\n");}
    if(it==maxit){
        conv[0]+=1; 
        warning("The maximum iteration has been reached \n with del %f.\n", del);}
    Free(Bta);
    Free(pBeta);
    Free(fp);
    Free(pnu);
}
/* end function em_beta_mix */

void mable_mvar(int *m, int *n, int *d, int *km, double *p, double *x, 
        int *maxit,  double *eps, double *llik, int *progress, int *conv){
    Rprintf("\n Mable fit of multivariate data. This may take several minutes.\n\n");
    em_multivar_beta_mix(p, x, m, *n, *d, km, *maxit, *eps, llik, *progress, conv);
}
/* t=(t1,...,td) is in [0,1]^d*/
void mable_mvdf(int *d, int *m, int *km, int *n, double *t, double *p, 
            double *mvdf, int *density){
    int i, j, K;
    double *tmp;
    K=km[*d-1];
    tmp = Calloc((*n)*K, double);
    if(*density==0) Multivar_pBeta(t, m, *n, *d, km, tmp);
    if(*density==1) Multivar_dBeta(t, m, *n, *d, km, tmp);
    for(i=0;i<*n;i++){
        mvdf[i]=0.0;
        for(j=0;j<K;j++) mvdf[i]+=p[j]*tmp[j+K*i];
    }
    Free(tmp);
}
// end of mable-multivar.c
