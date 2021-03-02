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
void ProgressBar (double percentage, char *txt){
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    Rprintf("\r%s%3d%% [%.*s%*s]",txt, val, lpad, PRGRSS, rpad, "");
//    Rprintf("\r%3d%% [%.*s%*s]\n%s", val, lpad, PRGRSS, rpad, "", txt);
    R_FlushConsole();
}

//void clockProgress (double percentage)
void clockProgress (int val, char *txt)
{
    const char *pcl[8]={"|","/","-","\\","|","/","-","\\"};
//    int val = (int) (percentage * 100);
    int agl = (int) (val % 8);
    Rprintf ("\r%s %s", pcl[agl], txt);
    R_FlushConsole();
}

/*////////////////////////////*/
/* Ordinal indicator function */
/*////////////////////////////*/
const char *Ord(int i){
    static char ord[3];
    int it=i%10;
    if(i!=11 && i!=12 && i!=13){
        if(it==1) strcpy(ord, "st");
        if(it==2) strcpy(ord, "nd"); 
        if(it==3) strcpy(ord, "rd");
    }
    else strcpy(ord, "th"); 
    return ord;
}

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
/*///////////////////////////////////////////////////////////////////*/
/*  Calculate the integral of B_{m1,i}(x)*B_{m2,j}(x) from 0 to 1,   */  
/*  for i in 0:m1, j in 0:m2; B_{mi} is the CDF of beta(i+1, m-i+1)  */
/*///////////////////////////////////////////////////////////////////*/
/* B_{m1,i}(x)*B_{m2,j}(x) */
static void Bm1ixBm2j(double *x, int n, void *ex)
{
    int i, j, m1, m2, k, *par;
    par = (int *) ex;
    m1 = par[0];
    m2 = par[1];
    i = par[2];
    j = par[3];
    for(k = 0; k < n; k++){
	   x[k] = pbeta(x[k],i+1,m1-i+1,TRUE,FALSE)*pbeta(x[k],j+1,m2-j+1,TRUE,FALSE);
    }
    return;
}
/* integral_0^1 B_{m1,i}(x)*B_{m2,j}(x) dx */
void int_Bm1xBm2(int m1, int m2, double *B){
    int i, j, m1p1=m1+1, *ex;
    double lo=0.0, up=1.0, epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400]; 
    int lenw=400, last=0, neval=0, ier=0, iwork[100]; 
    int limit=100; 
    ex = Calloc(4, int);
    ex[0] = m1;
    ex[1] = m2;
    for(i=0; i<=m1; i++){
        ex[2] = i; 
        for(j=0; j<=m2; j++){
            ex[3] = j;
            Rdqags(Bm1ixBm2j, ex, &lo, &up, &epsabs, &epsrel, 
                &result, &abserr, &neval, &ier, &limit, &lenw, &last, 
                iwork, work);
            B[i+j*m1p1] = result;
            //Rprintf("i=%d, j=%d, m1=%d, m2=%d, B=%g\n",i,j,m1,m2, result);
        }
    }
    Free(ex);
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
    return llik;
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
    return llik;
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
      int *convergence, double *delta, double *level, int *vb){
  int d, i, j, l, m, *cp, lp, tmp=0, *diverge, k=M[1]-M[0], cp0=0, cp1=1, i0=0, i1=0; 
  double *phat, *Bta, *llik, *res, pv0=1.0, pv1=1.0;//, max_bic; 
  double pct=0.0, ttl;
  lp=M[0]*(k+1)+(k+1)*(k+2)/2;
  phat = Calloc(lp, double);
  diverge = Calloc(1, int);
  cp = Calloc(1, int);
  res = Calloc(1, double);
  llik = Calloc(1, double);
  Bta = Calloc((M[1]+1)*(*n), double);
  ttl=(double) (k+2)*(k+1);

  if(*vb==-1 || *vb==2) i0=1;
  if(*vb== 1 || *vb==2) i1=1;
  m=M[0];
  if(m>0){
    //for(j=0; j<=m; j++) p[j]=1.0/(double) (m+1);
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
      for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(*vb));
      if(*vb==-1 || *vb==2) p[0]=0.0;
      if(*vb== 1 || *vb==2) p[m]=0.0;
    }
    dBeta(x, m, *n, Bta);
    em_beta_mix(p, Bta, m, *n,  *maxit, eps[0], llik, diverge, delta);
    convergence[0]+=diverge[0];
    tmp=m+1;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    lk[0] = *llik;
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;
    bic[0]=lk[0]-.5*d*log(*n);
  }
  else{
    phat[0]=1.0;
    lk[0]=0.0;
    bic[0]=0.0;
  }
  pval[0]=1.0;
  chpts[0]=0;
  pct = fmax2(*level/pval[0], pct);
  pct = fmax2(pct, 2.0/ttl);
  if(*progress==1) ProgressBar(fmin2(1.0,pct),"");
  i=1;
  //stop if the p-value for change-point is small 
  while(i<=k && pval[i-1]>*level){
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
      for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(*vb));
      if(*vb==-1 || *vb==2) p[0]=0.0;
      if(*vb== 1 || *vb==2) p[m]=0.0;
    }
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
    if(m>3) for(j=i0; j<=m-i1; j++) p[j] =(1.0-*tini)*p[j]+ *tini/(double)(m+1-abs(*vb));
    //for(j=0; j<=m; j++) p[j] =(p[j]+ *tini/(double)(m+1.0))/(1.0+*tini);
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
      chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
      cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
      cp0=cp1;
      pv0=pv1;
    }
    else pv0=pval[i];
    R_CheckUserInterrupt();
    // Calculate BIC
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;
    bic[i]=lk[i]-.5*d*log(*n);
    pct = fmax2(pct, *level/pval[i]);
    pct = fmax2(pct, i*(i+1)/ttl);
    if(*progress==1) ProgressBar(fmin2(1.0,pct),"");
    i++; 
  }
  if(convergence[0]>0) convergence[0]=1;
  if(*progress==1){
    ProgressBar(1.0,"");
    Rprintf("\n");}
  if(m==M[1] && pval[i-1]>*level){
    convergence[0]+=1; 
    Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value,  %f, of the change-point is returned.\n", pv0);
//        warning("The maximum candidate degree has been reached \n with a p-value of the change-point %f.\n", res[0]);
  }
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  optim[0]=cp0+M[0];
  m=optim[0];
  for(j=0;j<=m;j++) p[j]=phat[tmp+j];
  Free(phat); Free(Bta); Free(llik); 
  Free(cp); Free(res);
}
/*//////////////////////////////////////////////////////////*/
void mable_optim_group(int *M, int *N, double *p, double *t, int *n, int *maxit, double *eps, 
      double *lk, double *lr, int *optim, int *progress, int *convergence, double *delta, 
      double *tini, double *bic, double *pval, int *chpts, double *level, int *vb){
  int d, nn=0, i, j, m, k, lp, tmp=0, *cp, *diverge, cp0=0, cp1=1, i0=0, i1=0; 
  double *phat, *dBta, pct=0.0, *mlik, ttl, *res, pv0=1.0, pv1=1.0;//, max_bic;//maxlr, tmp, 
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
  if(*vb==-1 || *vb==2) i0=1;
  if(*vb== 1 || *vb==2) i1=1;
  m=M[0];
  if(m>0){
    //for(j=0; j<=m; j++) p[j]=1.0/(double)(m+1);
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
      for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(*vb));
      if(*vb==-1 || *vb==2) p[0]=0.0;
      if(*vb== 1 || *vb==2) p[m]=0.0;
    }
    cpBeta(t, m, *N, dBta);
    em_beta_mix_group(p, dBta, *N, m, n, *maxit, eps[0], mlik, diverge, delta);
    convergence[0]+=diverge[0];
    tmp=m+1;
    for(i=0;i<tmp;i++) phat[i]=p[i];
    lk[0] = mlik[0];
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;
    bic[0]=lk[0]-.5*d*log(nn);
  }
  else{
    phat[0]=1.0;
    lk[0]=0.0;
    bic[0]=0.0;
  }
  pval[0]=1.0;
  chpts[0]=0;
  pct = fmax2(pct,*level/pval[0]);
  pct = fmax2(pct,2/ttl);
  if(*progress==1) ProgressBar(pct,""); 
  i=1;
  while(i<=k && pval[i-1]>*level){
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
      for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(*vb));
      if(*vb==-1 || *vb==2) p[0]=0.0;
      if(*vb== 1 || *vb==2) p[m]=0.0;
    }
    // updating p
    p[m+1] = (m+1)*p[m]/(double)(m+2);
    for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
    p[0] = (m+1.)*p[0]/(double)(m+2.);
    m=M[0]+i;
    cpBeta(t, m, *N, dBta);
    // make sure initial p is in the interior of the simplex
    if(m>3) for(j=i0; j<=m-i1; j++) p[j] =(1.0-*tini)*p[j]+ *tini/(double)(m+1-abs(*vb));
    //for(j=0; j<=m; j++) p[j] =(p[j]+ *tini/(double)(m+1.0))/(1.0+*tini);
    //for(j=0;j<=m;j++) p[j]=1.0/(double) (m+1);
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
      chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
      cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
      cp0=cp1;
      pv0=pv1;
    }
    else pv0=pval[i];
    R_CheckUserInterrupt();
    // Calculate BIC
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;
    bic[i]=lk[i]-.5*d*log(nn);
    pct = fmax2(pct,*level/pval[i]);
    pct = fmax2(pct,i*(i+1)/ttl);
    if(*progress==1)ProgressBar(fmin2(1.0,pct),"");
    i++;
  }
  if(*progress==1){
    ProgressBar(1.0,"");
    Rprintf("\n");}
  if(convergence[0]>0) convergence[0]=1;
  if(m==M[1]){
    convergence[0]+=1; 
    Rprintf("\n The maximum candidate degree has been reached \n with a p-value of the change-point %f.\n", res[0]);}
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  optim[0]=cp0+M[0];
  //    tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
  //    optim[0]=cp[0]+M[0];
  m=optim[0];
  for(j=0;j<=m;j++) p[j]=phat[tmp+j];
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
    for(j=0;j<*n;j++){
       tmp=.0;
       for(i=0;i<=*m;i++){
          tmp += Bta[j+*n*i]*p[i];
       }
       u[j]=tmp;
    }
    Free(Bta);
}
/*//////////////////////////////////////////*/
/*   Generating PRN from beta(i+1, m-i+1)   */
/*  n: sample size,                         */
/*  w: vector of i values between 0 and m   */
/*  v: the generated variates               */
/*//////////////////////////////////////////*/
void rbeta_mi(int *n, int *m, int *w, double *v){
    int i, j, mp2=*m+2;
    double tmp1, tmp2;
    for(j=0; j<*n; j++) {
        tmp1 = 1.0;
        tmp2 = 1.0;
        for(i=0; i<w[j]+1; i++)
            tmp1 *= unif_rand();
        tmp2 *= tmp1;
        for(i=w[j]+1; i<mp2; i++)
            tmp2 *= unif_rand();
        v[j] = log(tmp1)/log(tmp2);
    }
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
  return llkhd;
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
    //if (y2[k]<=1 && z2[k]>1) error("\n Error: z2 >1 \n");
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
  //for(i=0; i<d; i++)
     //for(j=0;j<d;j++) Rprintf("\n ddell[%d,%d]=%f", i,j, ell[i+d*j]);

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
  int i,j, d=dm[0], k=M[1]-M[0], *cp, cp0=1, cp1=1, n0=N[0], n1=N[1], n=n0+n1, tmp, itmp=0;
  int m=M[1], mp1=m+1,  mp2=m+2, lp=(k+1)*M[0]+(k+1)*(k+2)/2;//optim=m,
  double tini=.0001, pct=0.0, ttl, *z, *z2, *res, pv0=1.0, pv1=1.0;
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
  if(*progress==1) {Rprintf("\n Mable fit of AFT model with given regression coefficients ... \n");
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
      chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
      cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
      cp0=cp1;
      pv0=pv1;
    }
    else pv0=pval[i];
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
      Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value,  %f, of the change-point is returned.\n", pv0);
      //warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);
      delta[0]=res[0];
    }
  }
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  dm[1]=cp0+M[0];
  //tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
  //dm[1]=cp[0]+M[0];
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
    if(y2[i]<=1 && z2[i]>1){
      Rprintf("\n");
      error("Try another baseline 'x0' and/or a larger truncation time 'tau'.\n");}
  }
  Bdata(z, m, 0, n, BSz);
  Bdata(z2, m, n0, n1, BSz2);
  if(*progress==1){
    Rprintf("\n Mable fit of AFT model with a given degree ... \n"); 
    ProgressBar(pct,"");} 
  pofg_aft(p, m, gx, n0, n1, BSz, BSz2, ell, eps_em, maxit_em, prog, conv, delta);
  //Rprintf("\n ell=%f, ell0=%f\n",ell[0], ell[1]);
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
    //Rprintf("\n ell=%f, ell0=%f\n",ell[0], ell1[0]);
    //Rprintf("\n pmp1=%d, x0=%f,  %f\n", *pmp1, x0[0], x0[1]);
    //pct=fmax(1-fabs(del-eps)/(.00001+eps),1-fabs(ell[1]-ell[0])/fabs(ell[1]));
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
    //warning("\nThe maximum iterations were reached \nwith a delta = %f.\n", del);
  }
  //Rprintf("mable-m: it=%d, del=%f\n", it, del);
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
  int i, j, d=dm[0], k=M[1]-M[0], tmp=0,*cp, cp0=1, cp1=1;//, n0=N[0], n1=N[1], n=n0+n1
  int m=M[1], mp1=m+1, lp=(k+1)*M[0]+(k+1)*(k+2)/2, prg=1-*progress; 
  double *phat, *ghat, *ell, pct=0.0, ttl, lnn=-1.0e20, *res, pv0=1.0, pv1=1.0;  //maxLR=0.0, lr0, 
  cp = Calloc(1, int);
  res = Calloc(1, double);
  phat=Calloc(lp, double);
  ghat=Calloc(d*(k+1), double);
  ell=Calloc(2, double);
  //egx=Calloc(n, double);
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
      chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
      cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
      cp0=cp1;
      pv0=pv1;
    }
    else pv0=pval[i];
    R_CheckUserInterrupt();
    pct +=2*(i+1);
    if(*progress==1) ProgressBar(pct/ttl,""); 
    i++;
  }
  if(*progress==1){
    ProgressBar(1.00,"");
    Rprintf("\n");}
  // Rprintf("mable-aft done!\n"); 
  if(m==M[1]){
    conv[0]+=1; 
    Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value of the change-point %f is returned.\n", pv0); 
    //warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);
  }
  //else conv[0]=0;
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  dm[1]=cp0+M[0];
  //tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
  //dm[1]=cp[0]+M[0];
  m=dm[1];
  for(j=0;j<=m;j++) p[j]=phat[tmp+j];
  for(j=0; j<dm[0]; j++) gama[j]=ghat[dm[0]*cp0+j];
  //for(j=0; j<dm[0]; j++) gama[j]=ghat[dm[0]*cp[0]+j];
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
    return llkhd;
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
  int *cp, lp, tmp, m=M[1],  mp1=m+1, mp2=m+2, itmp=0, cp0=1, cp1=1; 
  double tini=.000001, pct, ttl, *res, pv0=1.0, pv1=1.0;
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
  if(*progress==1) {
    Rprintf("\n Mable fit of PH model with given regression coefficients ... \n");
    ProgressBar(0.0,""); }
  ttl = (double)((k+2)*(k+1));
  //egx_x0(gama, d, x, n, egx, x0);
  egxmx0(gama, d, x, n, egx, x0);
  // add check to see if any egx is less than 1
  for(i=0;i<n;i++) 
    if(egx[i]<1) {
      Rprintf("\n");
      error("Try another baseline 'x0'.\n");}
  m=M[0]; 
  mp1=m+1;
  mp2=m+2;
  for(i=0;i<=m;i++) p[i]=*pi0/(double) mp1; 
  p[mp1]=1.0-*pi0;
  if(m>0){
    Bdata(y, m, 0, n, BSy);
    Bdata(y2, m, n0, n1, BSy2);
    pofg_ph(p, m, egx, n0, n1, BSy, BSy2, ell, *eps, *maxit, *progress, conv, delta);
    itmp+=conv[0];
    lk[0]=ell[0]; 
  }
  else{
    lk[0]=0;
    for(i=0;i<n0;i++) lk[0]+=log(egx[i])+(egx[i]-1)*log(1.0-y[i]);
    for(i=n0;i<n;i++) lk[0]+=log(R_pow(1.0-y[i],egx[i])-R_pow(1.0-y2[i],egx[i]));
  }
  tmp=mp2;
  for(i=0;i<tmp;i++) phat[i]=p[i];
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
    //Rprintf("lk[%d]=%f\n",i, lk[i]);
    if(i>=2){
      cp[0]=i;
      chpt_exp(lk, lr, res, cp);
      pval[i]=res[0];
      chpts[i]=cp[0];
    }
    else{            
      pval[i]=1.0;
      chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
      cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
      cp0=cp1;
      pv0=pv1;
    }
    else pv0=pval[i];
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
      Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value,  %f, of the change-point is returned.\n", pv0);
    }
    delta[0]=res[0];
    delta[1]=pv0;
  }
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+3))/2;
  dm[1]=cp0+M[0];
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
  //egx_x0(gama, d, x, n, egx, x0);
  egxmx0(gama, d, x, n, egx, x0);
  // add check to see if any egx is less than 1
  for(i=0;i<n;i++) 
    if (egx[i]<1) {
      Rprintf("\n");
      error("Try another baseline 'x0'.\n");}
  for(j=0;j<d;j++)  gnu[j]=gama[j];
  //Rprintf("gnu=%f\n", gnu[0]);
  if(m>0) pofg_ph(p, m, egx, n0, n1, BSy, BSy2, ell, eps_em, maxit_em, prog, conv, delta);
  gofp_ph(gnu, d, p, m, x, x0, n0, n1, BSy, BSy2, ell, dell, ddell, eps_nt, maxit_nt, prog);
  del=0.0;
  for(i=0;i<d;i++){
    del+=fabs(gnu[i]-gama[i]);
    gama[i]=gnu[i];
  }
  if(m==0) del=0.0;
  if(*progress==1){
    Rprintf("\n Mable fit of PH model with a given degree ... \n"); 
    ProgressBar(pct,""); }
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
    //Rprintf("         mable-m: it=%d, del=%f\n", it, del);
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
  //Rprintf("mable-m: it=%d, del=%f\n", it, del);
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
  int m, *cp, tmp, lp, cp0=1, cp1=1;  
  double *ghat, *phat, *res, *ell, pct, ttl, pv0=1.0, pv1=1.0; 
  lp=M[0]*(k+1)+(k+1)*(k+4)/2;
  cp = Calloc(1, int);
  res = Calloc(1, double);
  phat=Calloc(lp, double);
  ghat=Calloc(d*(k+1), double);
  ell=Calloc(1, double);
  //egx=Calloc(n, double);
  if(*progress==1) {Rprintf("\n Mable fit of Cox PH regression model ... \n");
      ProgressBar(0.0,""); }
  ttl=(double)(k+2)*(k+1);
  m=M[0]; 
  for(i=0;i<=m;i++) p[i]=*pi0/(double)(m+1);
  p[m+1] = 1-*pi0;
  dm[1]=m;
  mable_ph_m(gama, p, dm, x, y, y2, N, x0, ell, ddell, EPS, MAXIT, &prg, conv, res);
  for(i=0;i<dm[0];i++) ghat[i]=gama[i];
  lk[0]=ell[0];
  tmp=m+2;
  for(i=0;i<tmp;i++) phat[i]=p[i];
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
    //Rprintf("lk[%d]=%f\n",i, lk[i]);
    if(i>=3){
      cp[0]=i;
      chpt_exp(lk, lr, res, cp);
      pval[i]=res[0];
      chpts[i]=cp[0];
    }
    else{            
      pval[i]=1.0;
      chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
      cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
      cp0=cp1;
      pv0=pv1;
    }
    else pv0=pval[i];
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
    Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value of the change-point %f is returned.\n", pv0);}
    //warning("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);}
  //else conv[0]=0;
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+3))/2;
  dm[1]=cp0+M[0];
  //tmp=cp[0]*(M[0]*2+(cp[0]+3))/2;
  //dm[1]=cp[0]+M[0];
  m=dm[1];
  for(j=0;j<=m+1;j++) p[j]=phat[tmp+j];
  //for(j=0; j<dm[0]; j++) gama[j]=ghat[dm[0]*cp[0]+j];
  for(j=0; j<dm[0]; j++) gama[j]=ghat[dm[0]*cp0+j];
  if(*progress==1) Rprintf("\n");
  Free(phat);
  Free(ghat);
  Free(ell); 
  //Free(egx);
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
  //egx=Calloc(n, double);
  Tmp=Calloc(n*mp2, double);
  Tmp2=Calloc(n*mp2, double);
  pnu=Calloc(mp2, double);
  //for(i=0;i<n;i++)  egx[i]=1.0;
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
  //ell[0]=log_blik_ph(p, m, egx, n0, n1, BSy, BSy2); 
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
  //max_bic=bic[0];
  //optim[1]=0;
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
    //Rprintf("lk[%d]=%f\n",i, lk[i]);
    if(i>=3){
      cp[0]=i;
      chpt_exp(lk, lr, res, cp);
      pval[i]=res[0];
      chpts[i]=cp[0];
    }
    else{            
      pval[i]=1.0;
      chpts[i]=0;
    }
    // Calculate BIC
    d=0;
    for(j=0;j<=mp1;j++) d+=1*(p[j]>=eps[1]);
    d=d-1;// bic not working
    bic[i]=lk[i]-.5*d*log(n);
    //if(bic[i]>max_bic){
    //    max_bic=bic[i];
    //    optim[1]=i;
    //}
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
        Rprintf("\nThe maximum candidate degree has been reached \nwith a p-value of the change-point %f.\n", res[0]);
        delta[0]=res[0];
    }
    //tmp=optim[1]*(M[0]*2+(optim[1]+3))/2;
    //optim[1]+=M[0];
    //m=optim[1];
    //for(j=0;j<=m+1;j++) p[optim[0]+2+j]=phat[tmp+j];
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

typedef struct int_struct{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
} int_struct, *IntStruct;
typedef struct mable_struct{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
    int m;
    int j;
    double y;
} mable_struct, *MableStruct;

/*//////////////////////////////////////////////////////////////////////*/
/*          The integrand eta(x,y;j,m)=g(y-x)*beta_mj(x) of the         */
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
            Rdqags(eta_mj, (void*)&ms, &l, &u, &epsabs, &epsrel, &result, 
                &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
            psi_m[i+j*n] = result;
        }
    }
}
/*//////////////////////////////////////////////////*/
/*          EM Algorithm for a fixed m              */
/*//////////////////////////////////////////////////*/
/* EM Method for mixture of g*beta(i+1, m+1-i), i=0,...,m. n=sample size,  */
void em_gBeta_mix(double *y, double *p, int m, int n, int maxit, double eps, 
        double *llik, void *ex){
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
    if(m>0) del = 10.0;
    else del=0.0;
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
SEXP mable_decon(SEXP args){
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
    M = INTEGER(CAR(args)); args = CDR(args);
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
    d=0;
    for(j=0;j<=m;j++) d+=1*(p[j]>=eps);
    d=d-1;
    bic[0]=lk[0]-.5*d*log(n);
    pval[0]=1.0;
    chpts[0]=0;
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
            chpts[i]=0;
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
/*//////////////////////////////////////////////////////////////////////*/
/*       Density Estimation based on data with measurement errors       */
/*//////////////////////////////////////////////////////////////////////*/
/*//////////////////////////////////////////////////////////////////////*/
/*   The integrand eta(x,y;i,j,m,k)=beta_mi(x)*beta_kj((y-x-c)/(d-c))   */
/*    of the convolution of beta_{mi} and beta_{kj}((y-x-c)/(d-c))      */
/*    returns a function of x for fixed y, ex = (m, k, i, j, y, c, d)   */
/*           code modified from 'Rintfn()' in 'integrate.c'             */
/*//////////////////////////////////////////////////////////////////////*/
static void eta_mkij(double *x, int n, void *ex)
{
    int i, j, m, k, l;
    double c, d, dmc, y, z, *par;
    par =(double*) ex;
    m = (int)par[0];
    k = (int)par[1];
    i = (int)par[2];
    j = (int)par[3];
    y = par[4];
    c = par[5];
    d = par[6];
    dmc=d-c;
    for(l = 0; l < n; l++){
	   z = (y-x[l]-c)/dmc;
	   x[l] = dbeta(x[l], i+1,  m-i+1, FALSE)*dbeta(z, j+1,  k-j+1, FALSE);
    }
    return;
}
/*////////////////////////////////////////////////////////////////////////*/
/*       Calculate the convolution gam_mk(y)=beta_{mj}*beta(y),           */  
/*                          Input: (m, y)                                 */
/*////////////////////////////////////////////////////////////////////////*/

void gamma_mk(double *y, double c, double d, double *gam, 
        int m, int k, int n){
    int i, j, l, mp1=m+1, kp1=k+1;
    double lo, up, *ex, dmc=d-c, epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400]; 
    int lenw=400, last=0, neval=0, ier=0, iwork[100]; 
    int limit=100; 
    ex = Calloc(7, double);
    ex[0] = m;
    ex[1] = k;
    ex[5] = c;
    ex[6] = d;
    for(l=0; l<n; l++){
        ex[4] = y[l];
        lo = fmax2(0.0, y[l]-d); 
        up = fmin2(1.0, y[l]-c);
        if(lo<up){
            for(i=0; i<=m; i++){
                ex[2] = i; 
                for(j=0; j<=k; j++){
                    ex[3] = j;
                    Rdqags(eta_mkij, ex, &lo, &up, &epsabs, &epsrel, 
                        &result, &abserr, &neval, &ier, &limit, &lenw, &last, 
                        iwork, work);
                    gam[i+j*mp1+l*mp1*kp1] = result/dmc;
                    //Rprintf("i=%d, j=%d, m=%d, k=%d, gam=%g\n", i, j,m,k, result/dmc);
                }
            }
        }
        else
            for(i=0; i<=m; i++)
                for(j=0; j<=k; j++)
                    gam[i+j*mp1+l*mp1*kp1] = 0.0;
    }
    Free(ex);
}
// psi(y; p, q): a mixture of betas
void psi_pq(double *gam, double *p, double *q, double *psi, 
        double *spg, double *sqg, int n, int m, int k){
    int mp1=m+1, kp1=k+1, i, j, l;
    for(i=0; i<=m; i++){
        for(l=0; l<n; l++){
            sqg[i+l*mp1]=0.0;
            for(j=0; j<=k; j++)
                sqg[i+l*mp1] += q[j]*gam[i+j*mp1+l*mp1*kp1];
        }
    }
    for(j=0; j<=k; j++){
        for(l=0; l<n; l++){
            spg[j+l*kp1]=0.0;
            for(i=0; i<=m; i++) 
                spg[j+l*kp1] += p[i]*gam[i+j*mp1+l*mp1*kp1];
        }
    }
    for(l=0; l<n; l++){
        psi[l]=0.0;
        for(i=0; i<=m; i++)
            psi[l] += p[i]*sqg[i+l*mp1];
    }
    //if(l==0) Rprintf("psi[%d]=%g\n", l, psi[l]);
}
// Update p or q: p[s+1]=p[s]*mean(sqg/psi), q[s+1]=q[s]*mean(spg/psi),
void new_pq(double *p, double *psi, double *sqg, int n, int m){
    int mp1=m+1, i, l;
    double w;
    for(i=0; i<=m; i++){
        w = 0.0;
        for(l=0; l<n; l++) w += sqg[i+l*mp1]/psi[l];
        p[i] *= w/(double)n;
    }
}
// Find Lagrange multiplier of mean constraint
// sum((c+(d-c)*(j+1)/(k+2)-mu)q[j])=0
double lgrg_mltpl(double eta, double *q, double *Ck, int k, 
        double eps, int maxit){
    int it=0, j;
    double del, h=0.0, dh=0.0, wk, tmp;
    for(j=0; j<=k; j++){
        wk = Ck[j]/(1.0+eta*Ck[j]);//??
        h += wk*q[j];
        dh -= q[j]*wk*wk;
    }
    del = fabs(h);
    while(it<maxit && del>eps){
        tmp = h/dh;
        del = fabs(tmp);
        eta -= tmp;
        h = 0.0;
        dh = 0.0;
        for(j=0; j<=k; j++){
            wk = Ck[j]/(1.0+eta*Ck[j]);//??
            h += wk*q[j];
            dh -= q[j]*wk*wk;
        }
        del += fabs(h);
        it++;
        // Rprintf("Lagrange multiplier: it=%d, del=%g, eta=%g\n", it,  del,  eta);
    }
    //Rprintf("Lagrange multiplier: it=%d, del=%g, eta=%g\n", it,  del,  eta);
    return eta;
}
// Squared difference between estimated convolutions 
//   with and without deconvolution
static void sq_diff_psi(double *x, int n, void *ex){
    int i, j, l, m, k, v, mp1, kp1;
    double c, d, *pi, *p, *q, dmcp1, z, *gam, tmp1, tmp2, *par;
    par = (double*) ex;
   m = (int)par[0];
   k = (int)par[1];
   v = (int)par[2];
   c = par[3];
   d = par[4];
    mp1 = m+1;
    kp1 = k+1;
    gam=Calloc(mp1*kp1*n, double);
    pi=Calloc(v+1, double);
    p=Calloc(mp1, double);
    q=Calloc(kp1, double);
    dmcp1=d-c+1.0;
   for(i = 0; i <= v; i++) pi[i] = par[5+i];
   for(i = 0; i <= m; i++) p[i] = par[v+6+i];
   for(i = 0; i <= k; i++) q[i] = par[v+m+7+i];
    gamma_mk(x, c, d, gam, m, k, n);
    for(l = 0; l < n; l++){
        tmp1 = 0.0;
        tmp2 = 0.0;
        for(i = 0; i <= m; i++)
            for(j = 0; j <= k; j++)
                tmp1 += p[i]*q[j]*gam[i+mp1*j+mp1*kp1*l];
        for(i = 0; i <= v; i++){
	        z = (x[l]-c)/dmcp1;
            tmp2 += pi[i]*dbeta(z, i+1,  v-i+1, FALSE);
        }
        tmp2 /= dmcp1;
	    x[l] = (tmp1-tmp2)*(tmp1-tmp2);
    }
    Free(gam);
    Free(pi);
    Free(p);
    Free(q);
    return;
}   
// L2 distance between estimated convolutions 
//   with and without deconvolution
double d_btwn_psi(double *pi, int v, double *p, int m, double *q, int k, 
        double c, double d){
    int i;
    double lo, up, *ex;
    double epsabs=.00001, epsrel=.00001, result=0.0, abserr=0.0, work[400]; 
    int lenw=400, last=0, neval=0, ier=0, iwork[100]; 
    int limit=100; 
    ex = Calloc(8+v+m+k, double);
    ex[0] = m;
    ex[1] = k;
    ex[2] = v;
    ex[3] = c;
    ex[4] = d;
   for(i = 0; i <= v; i++) ex[5+i] = pi[i];
   for(i = 0; i <= m; i++) ex[v+6+i] = p[i];
   for(i = 0; i <= k; i++) ex[v+m+7+i] = q[i];
    lo = c; 
    up = d+1.0;
    Rdqags(sq_diff_psi, ex, &lo, &up, &epsabs, &epsrel, &result, 
        &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);
    //Rprintf("result=%g\n", result);   
    Free(ex); 
    return result;
} 
// EM for deconvolution with unknown error density
void mablem_decon(double *gam, int n, double *interval, int m, int k, double *lk,
        double *p, double *q, int mConstr, double ybar, double eps, int maxit, 
        double eps_nt, int maxit_nt){
    int it=0, i, j, l, mp1=m+1, kp1=k+1;
    double c=interval[0], d=interval[1], eta=0.0, zeta=0.0, del=0.0;
    double *Am, *Ck, *psi, *spg, *sqg, *p1, *q1;
    Am = Calloc(mp1, double);//??
    Ck = Calloc(kp1, double);//??
    p1 = Calloc(mp1, double);
    q1 = Calloc(kp1, double);
    psi = Calloc(n, double);
    spg = Calloc(n*kp1, double);
    sqg = Calloc(n*mp1, double);
//Rprintf("c=%g, d=%g\n", c, d);
    for(i=0; i<=m; i++){
//        p[i]=1.0/(double)mp1;
        Am[i] = (i+1)/(double)(m+2)-ybar;
        p1[i]=p[i];
    }
    for(j=0; j<=k; j++){
        Ck[j]=c+(d-c)*(j+1)/(double)(k+2);
        q1[j]=q[j];
    }
    psi_pq(gam, p, q, psi, spg, sqg, n, m, k);        
    new_pq(p1, psi, sqg, n, m);
    new_pq(q1, psi, spg, n, k);
    if(mConstr==1) zeta=lgrg_mltpl(zeta, p1, Am, m, eps_nt, maxit_nt);
    eta=lgrg_mltpl(eta, q1, Ck, k, eps_nt, maxit_nt);
//    Rprintf("eta=%g, zeta=%g\n", eta, zeta);
    for(i=0; i<=m; i++){
        p1[i]=p1[i]/(1.0+Am[i]*zeta);
        del+=fabs(p1[i]-p[i]);
        p[i]=p1[i];
    }
    for(j=0; j<=k; j++){
        q1[j]=q1[j]/(1.0+Ck[j]*eta);
        del+=fabs(q1[j]-q[j]);
        q[j]=q1[j];
    }
    while(it<maxit && del>eps){
        psi_pq(gam, p, q, psi, spg, sqg, n, m, k);        
        new_pq(p1, psi, sqg, n, m);
        new_pq(q1, psi, spg, n, k);
        if(mConstr==1) zeta=lgrg_mltpl(zeta, p1, Am, m, eps_nt, maxit_nt);
        eta=lgrg_mltpl(eta, q1, Ck, k, eps_nt, maxit_nt);
        del=0.0;
        for(i=0; i<=m; i++){
            p1[i]=p1[i]/(1.0+Am[i]*zeta);
            del+=fabs(p1[i]-p[i]);
            p[i]=p1[i];
        }
        for(j=0; j<=k; j++){
            q1[j]=q1[j]/(1.0+Ck[j]*eta);
            del+=fabs(q1[j]-q[j]);
            q[j]=q1[j];
        }
        it++;
        R_CheckUserInterrupt();
//    Rprintf("eta=%g, zeta=%g\n", eta, zeta);
        //Rprintf("it=%d, del=%g\n", it, del);
    }
    //Rprintf("m=%d, k=%d, it=%d, del=%g\n", m, k, it, del);
    //lk[0]=ell_pq(gam, p, q, n, m, k);
    lk[0]=0.0;
    for(l=0;l<n;l++) lk[0]+=log(psi[l]);
    Free(Am);
    Free(Ck);
    Free(p1);
    Free(q1);
    Free(psi);
    Free(spg);
    Free(sqg);
}
// Optimal Degrees for density deconvolution for r=1
// with unknown error density
SEXP optim_decon(SEXP args){
    int n, v, progress, vanished, mConstr=0, maxit_em, maxit_nt;
    SEXP ans, ansnames, yy;
    int *np, i, j, m, k, nij, *M,  itmp, *nu_d, *nu_aic, *nu_bic, kpm1; 
    double *llik, *pi, *interval, pct=0.0, ttl, eps_em, eps_nt, *p, *q, eps_p; 
    double *y, *gam, ybar=0.0, *lk, *D, min_D, *p_d, *q_d, c, d, tini=1.0;
    double *aic, *bic, max_aic, max_bic, *p_aic, *q_aic, *p_bic, *q_bic;  
    Rprintf("\n Deconvolution is runing. This may take several minutes.\n\n");
    args = CDR(args);
    yy = CAR(args); args = CDR(args);
    interval = REAL(CAR(args)); args = CDR(args);
    vanished = asInteger(CAR(args)); args = CDR(args);
    M = INTEGER(CAR(args)); args = CDR(args);
    pi = REAL(CAR(args)); args = CDR(args);
    v = asInteger(CAR(args)); args = CDR(args);
    eps_em = asReal(CAR(args)); args = CDR(args);
    maxit_em = asInteger(CAR(args)); args = CDR(args);
    progress = asInteger(CAR(args)); args = CDR(args);
    eps_nt = asReal(CAR(args)); args = CDR(args);
    maxit_nt = asInteger(CAR(args)); args = CDR(args);
    eps_p=1.0e-8;//????
    n = length(yy);
    y= REAL(yy);
    for(i=0; i<n; i++) ybar+=y[i];
    ybar /= (double)n;
    int m0=M[0], k0=M[1], m1=M[2], k1=M[3], mp1=m1+1, kp1=k1+1, nr=m1-m0+1, nc=k1-k0+1;
    gam = Calloc(n*mp1*kp1, double);
    llik = Calloc(1, double);
    np=Calloc(nr*nc, int);
    nu_d=Calloc(2, int);
    nu_aic=Calloc(2, int);
    nu_bic=Calloc(2, int);
    lk=Calloc(nr*nc, double);
    D=Calloc(nr*nc, double);
    aic=Calloc(nr*nc, double);
    bic=Calloc(nr*nc, double);
    p=Calloc(mp1, double);
    q=Calloc(kp1, double);
    p_d=Calloc(mp1, double);
    p_aic=Calloc(mp1, double);
    p_bic=Calloc(mp1, double);
    q_d=Calloc(kp1, double);
    q_aic=Calloc(kp1, double);
    q_bic=Calloc(kp1, double);
    min_D=1.0e+100;
    max_aic=-1.0e+100;
    max_bic=max_aic;
    c = interval[0];
    d = interval[1];
    //Rprintf("n=%d, m=%d, k=%d, m0=%d, k0=%d, nr=%d, nc=%d\n", n, m1, k1, m0, k0, nr, nc);
    ttl=(double)(nr*nc);
    //Rprintf("ok\n");
    for(k=k0; k<=k1; k++){
        if(vanished==1){
            q[0]=0.0;
            q[k]=0.0;
            kpm1=k-1;
        }
        else{
            kpm1=k+1;
            q[0]=1.0/(double)kpm1;
            q[k]=1.0/(double)kpm1;
        }
        for(j=1; j<k; j++)
            q[j]=1.0/(double)kpm1;
        mp1=m0+1;
        for(i=0; i<=m0; i++) 
            p[i]=1.0/(double)mp1;
        for(m=m0; m<=m1; m++){
            gamma_mk(y, c, d, gam, m, k, n);
            mablem_decon(gam, n, interval, m, k, llik, p, q, mConstr, 
                ybar, eps_em, maxit_em, eps_nt, maxit_nt);
            nij=m-m0+(k-k0)*nr;
            lk[nij]=*llik;
            itmp=0;
            for(i=0;i<=m;i++) itmp+=1*(p[i]>eps_p);
            np[nij]=imax2(1, itmp-1);
            itmp=0;
            for(j=0;j<=k;j++) itmp+=1*(q[j]>eps_p);
            np[nij]+=imax2(0, itmp-2);
            aic[nij]=lk[nij]-np[nij];
            bic[nij]=lk[nij]-np[nij]*.5*log(n);
            if(aic[nij]>max_aic){
                max_aic=aic[nij];
                nu_aic[0]=m;
                nu_aic[1]=k;
                for(i=0;i<=m;i++) p_aic[i]=p[i];
                for(j=0;j<=k;j++) q_aic[j]=q[j];
            }
            if(bic[nij]>max_bic){
                max_bic=bic[nij];
                nu_bic[0]=m;
                nu_bic[1]=k;
                for(i=0;i<=m;i++) p_bic[i]=p[i];
                for(j=0;j<=k;j++) q_bic[j]=q[j];
            }
            D[nij]=d_btwn_psi(pi, v, p, m, q, k, c, d);
            if(D[nij]<min_D){
                min_D=D[nij];
                nu_d[0]=m;
                nu_d[1]=k;
                for(i=0;i<=m;i++) p_d[i]=p[i];
                for(j=0;j<=k;j++) q_d[j]=q[j];
            }
            p[m+1] = (m+1)*p[m]/(double)(m+2);
            for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
            p[0] = (m+1.)*p[0]/(double)(m+2.);
            for(j=0; j<=m+1; j++) p[j] =(p[j]+tini/(double)(m+2.0))/(1.0+tini);
            if(vanished==1){
                q[0]=0.0;
                q[k]=0.0;
                for(j=1; j<k; j++) q[j] =(q[j]+tini/(double)(k-1.0))/(1.0+tini);
            }
            else
                for(j=0; j<=k; j++) q[j] =(q[j]+tini/(double)(k+1.0))/(1.0+tini);
            pct = (nij+1)/ttl;
            if(progress==1){
                ProgressBar(pct,"");}
        }
//        q[k+1] = (k+1)*q[k]/(double)(k+2);
//        for(j=k; j>=1; j--) q[j] = (q[j-1]*j+q[j]*(k-j+1))/(double)(k+2);
//        q[0] = (k+1.)*q[0]/(double)(k+2.);
//        for(j=0; j<=k+1; j++) q[j] =(q[j]+tini/(double)(k+1.0))/(1.0+tini);
    }
    if(progress==1) Rprintf("\n");
    PROTECT(ans = allocVector(VECSXP, 14));
    PROTECT(ansnames = allocVector(STRSXP, 14));
    SET_STRING_ELT(ansnames, 0, mkChar("lk"));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, nr*nc));
    SET_STRING_ELT(ansnames, 1, mkChar("np"));
    SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, nr*nc));
    SET_STRING_ELT(ansnames, 2, mkChar("aic"));
    SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, nr*nc));
    SET_STRING_ELT(ansnames, 3, mkChar("bic"));
    SET_VECTOR_ELT(ans, 3, allocVector(REALSXP, nr*nc));
    SET_STRING_ELT(ansnames, 4, mkChar("nu_aic"));
    SET_VECTOR_ELT(ans, 4, allocVector(INTSXP, 2));
    SET_STRING_ELT(ansnames, 5, mkChar("nu_bic"));
    SET_VECTOR_ELT(ans, 5, allocVector(INTSXP, 2));
    for(i=0; i<nr;i++){
        for(j=0; j<nc;j++){
            nij=i+j*nr;
            REAL(VECTOR_ELT(ans, 0))[nij] = lk[nij];
            INTEGER(VECTOR_ELT(ans, 1))[nij] = np[nij];
            REAL(VECTOR_ELT(ans, 2))[nij] = aic[nij];
            REAL(VECTOR_ELT(ans, 3))[nij] = bic[nij];
        }
    }
    for(i=0;i<2;i++){
        INTEGER(VECTOR_ELT(ans, 4))[i] = nu_aic[i];
        INTEGER(VECTOR_ELT(ans, 5))[i] = nu_bic[i];
    }    
    SET_STRING_ELT(ansnames, 6, mkChar("p_aic"));
    SET_STRING_ELT(ansnames, 7, mkChar("q_aic"));
    SET_VECTOR_ELT(ans, 6, allocVector(REALSXP, nu_aic[0]+1));
    SET_VECTOR_ELT(ans, 7, allocVector(REALSXP, nu_aic[1]+1));
    for(i=0;i<=nu_aic[0];i++)
        REAL(VECTOR_ELT(ans, 6))[i] = p_aic[i];
    for(j=0;j<=nu_aic[1];j++)
        REAL(VECTOR_ELT(ans, 7))[j] = q_aic[j];
    SET_STRING_ELT(ansnames, 8, mkChar("p_bic"));
    SET_STRING_ELT(ansnames, 9, mkChar("q_bic"));
    SET_VECTOR_ELT(ans, 8, allocVector(REALSXP, nu_bic[0]+1));
    SET_VECTOR_ELT(ans, 9, allocVector(REALSXP, nu_bic[1]+1));
    for(i=0;i<=nu_bic[0];i++)
        REAL(VECTOR_ELT(ans, 8))[i] = p_bic[i];
    for(j=0;j<=nu_bic[1];j++)
        REAL(VECTOR_ELT(ans, 9))[j] = q_bic[j];
    SET_STRING_ELT(ansnames, 10, mkChar("D"));
    SET_VECTOR_ELT(ans, 10, allocVector(REALSXP, nr*nc));
    for(i=0; i<nr;i++){
        for(j=0; j<nc;j++){
            nij=i+j*nr;
            REAL(VECTOR_ELT(ans, 10))[nij] = D[nij];
        }
    }
    SET_STRING_ELT(ansnames, 11, mkChar("nu_d"));
    SET_VECTOR_ELT(ans, 11, allocVector(INTSXP, 2));
    for(i=0;i<2;i++){
        INTEGER(VECTOR_ELT(ans, 11))[i] = nu_d[i];
    }    
    SET_STRING_ELT(ansnames, 12, mkChar("p_d"));
    SET_STRING_ELT(ansnames, 13, mkChar("q_d"));
    SET_VECTOR_ELT(ans, 12, allocVector(REALSXP, nu_d[0]+1));
    SET_VECTOR_ELT(ans, 13, allocVector(REALSXP, nu_d[1]+1));
    for(i=0;i<=nu_d[0];i++)
        REAL(VECTOR_ELT(ans, 12))[i] = p_d[i];
    for(j=0;j<=nu_d[1];j++)
        REAL(VECTOR_ELT(ans, 13))[j] = q_d[j];
    setAttrib(ans, R_NamesSymbol, ansnames);
    UNPROTECT(2);
    return ans;
//
    Free(np);
    Free(gam);
    Free(llik);
    Free(p);
    Free(q);
    Free(nu_d);
    Free(nu_aic);
    Free(nu_bic);
    Free(p_d);
    Free(p_aic);
    Free(p_bic);
    Free(q_d);
    Free(q_aic);
    Free(q_bic);
    Free(lk);
    Free(D);
    Free(aic);
    Free(bic);
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
/* km[0]=1, km[1]=m1+1,km[2]=(m1+1)(m2+1),...,km[d]=(m1+1)...(md+1)=K,*/
/*  Column-major layout  */
/* beta(x1j, i1+1, m1+1-i1)...beta(xdj, id+1, md+1-id), 0<=ik<=mk, 1<=k<=d,  
         the (i1+km[1]*i2+...+km[d-1]*id,j)th element */
void Multivar_dBeta(double *x, int *m, int n, int d, int *km, double *dBta){
  int i, j, jj, k, r, K, it;
  K=km[d];
  for(j=0; j<n; j++){
    for(it=0;it<K;it++){
      dBta[it+K*j]=1.0;
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        dBta[it+K*j]*=dbeta(x[j+n*k], i+1, m[k]+1-i, FALSE);
        r = jj;
        //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
      }
      dBta[it+K*j]*=dbeta(x[j], r+1, m[0]+1-r, FALSE); 
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
  K=km[d];
  for(j=0; j<n; j++){
    for(it=0;it<K;it++){
      pBta[it+K*j]=1.0;
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        pBta[it+K*j]*=pbeta(x[j+n*k], i+1, m[k]+1-i, TRUE, FALSE);
        r = jj;
      }
      pBta[it+K*j]*=pbeta(x[j], r+1, m[0]+1-r, TRUE, FALSE); 
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
  double lik, fx;
  lik = 1.0;
  for(j=0; j<n; j++){
    fx = 0.0;
    for(i=0; i<K; i++) fx += p[i]*Bta[i+K*j];
    lik *= fx;
  }
  //Rprintf("  lik = %g\n", lik);
  return log(lik);
}

/*/////////////////////////////////////////////////////*/
/* EM Method for mixture of                            */
/* beta(x1, i1+1, m1+1-i1),...,beta(xd, id+1, md+1-id),*/ 
/*               0<=ik<=mk, 1<=k<=d                    */
/*/////////////////////////////////////////////////////*/
void em_multivar_beta_mix(double *p, double *Bta, int *m, 
      int n, int d, int K, int maxit, double eps, 
      double *llik, int progress, int *conv){
  int i, j, it;
  double del, llik_nu, *pBeta, *fp, *pnu;
  double  ttl;
  ttl=(double) maxit;
  conv[0] = 0;
  pBeta = Calloc(K*n, double);
  fp = Calloc(n, double);
  pnu = Calloc(K, double);
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
    //Rprintf("  llik = %g\n", llik[0]);
    if(progress==1) ProgressBar(it/ttl,"");
  }
  if(progress==1){
    ProgressBar(1.0,"");
    Rprintf("\n");}
  if(it==maxit){
    conv[0]+=1; 
    if(progress==1) warning("\n The maximum iteration has been reached \n with del %f.\n", del);
  }
  Free(pBeta);
  Free(fp);
  Free(pnu);
}
/* end function em_beta_mix */
/*////////////////////////////////////////////////////////////*/
/* Calculate p of degree mt=m+e_k from p of degree m,         */
/* where ek=(0,...,0,1,0,...,0), the unit vector on xk-axis   */
/*   km[0] = 1,  km[1]=m1+1,  km[2]=(m1+1)(m2+1),...,         */
/*          km[d]=(m1+1)...(md+1)=K,                          */
/* p(i1, ..., id), 0<=ik<=mk, 1<=k<=d, are arranged by column-*/
/* major order:   i1+km[1]*i2+...+km[d-1]*id                  */
/*      For the p with m+ek, K1=K*(mk+2)/(mk+1)               */
/*////////////////////////////////////////////////////////////*/
void pm2pmpe_k(double *p, double *pt, int *mt, int d, int *kmt, int k){
  int i, j, l, r, K, K1, ii, ij, *I, itmp;
  I = Calloc(d, int);
  K1=kmt[d];
  K = K1*mt[k]/(mt[k]+1);
  for(i=0; i<K1; i++) p[i] = 0.0;
  for(ii=K-1; ii>=0; ii--){
    r = ii;
    ij = 0;
    for(l=d-1; l>k; l--){
      itmp = kmt[l]*mt[k]/(mt[k]+1);
      j = r%itmp;
      I[l] = (r-j)/itmp;
      ij += (r-j)*(mt[k]+1)/mt[k];
      r = j;
    }
    for(l=k; l>=0; l--){
      j = r%kmt[l];
      I[l] = (r-j)/kmt[l];
      ij += r-j;
      r = j;
    }
    //Rprintf("i=%d, r=%d\n", i, r);
    p[ij+kmt[k]] += (I[k]+1)*pt[ii]/(mt[k]+1);
    p[ij] += (mt[k]-I[k])*pt[ii]/(mt[k]+1);
  }
  Free(I);
}
// 
/*////////////////////////////////////////////////////////////*/
/* Calculate p(ij), the p's of the j-th marginal density from */
/*   p(i0, ..., i[d-1]), 0<=i[j]<=mj, 0<=j<d,  arranged in    */
/*   column-major order:   i0+km[1]*i1+...+km[d-1]*i[d-1]     */
/*   km[0] = 1,  km[1]=m0+1,  km[2]=(m0+1)(m1+1),...,         */
/*          km[d]=(m0+1)...(m[d-1]+1)=K,                      */
/*////////////////////////////////////////////////////////////*/
void p2pj(double *p, int *m, int d, int *km, double *pj, int j){
  int K=km[d], it, i, jj, k, r, *I;
  for(i=0;i<=m[j];i++) pj[i]=0.0;
  I=Calloc(d, int);
  it=0;
  while(it<K){
    r = it;
    for(k=d-1; k>0; k--){
      jj = r%km[k];
      i = (r-jj)/km[k];
      I[k]=i;
      r = jj;
      //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
    }
    I[0]=r;
    pj[I[j]]+=p[it];
    it++;
  }
  Free(I);
}
/*////////////////////////////////////////////////////////////*/
/* Calculate L2 distance D(F1,F2)=int_0^1[F1(x)-F2(x)]^2 dx   */
/*   between Bernstein CDFs F1 = F_{m1}(.; p1) and  F2 =      */
/*      F_{m2}(.; p2)                                            */
/*////////////////////////////////////////////////////////////*/
double L2_F1F2(double *p1, int m1, double *p2, int m2){
  int i, j, m1p1=m1+1, m2p1=m2+1, mmp1;
  double D=0.0, *B;
  mmp1=imax2(m1p1,m2p1);
  B=Calloc(mmp1*mmp1, double);
  int_Bm1xBm2(m1, m2, B);
  for(i=0;i<=m1;i++) 
    for(j=0;j<=m2;j++) D-=p1[i]*p2[j]*B[i+j*m1p1];
  D*=2.0;
  int_Bm1xBm2(m1, m1, B);
  for(i=0;i<=m1;i++){
    D+=p1[i]*p1[i]*B[i+i*m1p1];
    for(j=i+1;j<=m1;j++) D+=2*p1[i]*p1[j]*B[i+j*m1p1];
  }
  int_Bm1xBm2(m2, m2, B);
  for(i=0;i<=m2;i++){
    D+=p2[i]*p2[i]*B[i+i*m2p1];
    for(j=i+1;j<=m2;j++) D+=2*p2[i]*p2[j]*B[i+j*m2p1];
  }
  Free(B);
  return D;
}
/*////////////////////////////////////////////////////////////*/
/* Calculate L2 distance D(f1,f2)=int_0^1[f1(x)-f2(x)]^2 dx   */
/* between Bernstein densities f1 = f_{m1}(.; p1) and  f2 =   */
/*   f_{m2}(.; p2)                                            */
/*////////////////////////////////////////////////////////////*/
double L2_f1f2(double *p1, int m1, double *p2, int m2){
  int i, j;
  double D=0.0, dtmp;
  for(i=0;i<=m1;i++) {
    for(j=0;j<=m2;j++){
      dtmp=exp(lbeta(i+j+1, m1+m2-i-j+1)-lbeta(i+1,m1-i+1)-lbeta(j+1,m2-j+1));
      D-=p1[i]*p2[j]*dtmp;
    }
  }
  D*=2.0;
  for(i=0;i<=m1;i++) {
    dtmp=exp(lbeta(2*i+1,2*(m1-i)+1)-2*lbeta(i+1,m1-i+1));
    D+=p1[i]*p1[i]*dtmp;
    for(j=i+1;j<=m1;j++){
      dtmp=exp(lbeta(i+j+1,2*m1-i-j+1)-lbeta(i+1,m1-i+1)-lbeta(j+1,m1-j+1));
      D+=2*p1[i]*p1[j]*dtmp;
    }
  }
  for(i=0;i<=m2;i++) {
    dtmp=exp(lbeta(2*i+1,2*(m2-i)+1)-2*lbeta(i+1,m2-i+1));
    D+=p2[i]*p2[i]*dtmp;
    for(j=i+1;j<=m2;j++){
      dtmp=exp(lbeta(i+j+1,2*m2-i-j+1)-lbeta(i+1,m2-i+1)-lbeta(j+1,m2-j+1));
      D+=2*p2[i]*p2[j]*dtmp;
    }
  }
  return D;
}
/*////////////////////////////////////////////////////////////*/
/* Updating Multivariate Bernstein base polynomials when      */
/*     only the J-th degree changes by diff=1 or -1           */
/*                                                            */
/*////////////////////////////////////////////////////////////*/
void Update_dBeta(double *x, int *m, int n, int d, int J, int diff, int *km, double *dBta){
  int i, j, jj, k, r, K, it, jt, *I;
  K=km[d];
  I=Calloc(d, int);
  if(diff==1)
  for(j=n-1; j>=0; j--){
    it=K-1;
    while(it>=0){
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        I[k]=i;
        r = jj;
        //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
      }
      I[0]=r;
      jt=0;
      for(i=J+1;i<d;i++) jt+=I[i]*km[i];
      jt+=K*j;
      jt/=(m[J]+1);
      jt=it-jt+K*j;
      if(I[J]==m[J]) 
        dBta[it+K*j]=(m[J]+1)*x[j+n*J]*dBta[jt-km[J]]/(double)m[J];
      else dBta[it+K*j]=(m[J]+1)*(1.0-x[j+n*J])*dBta[jt]/(double)(m[J]-I[J]);
      it--;
    }
  }
  else
  for(j=0; j<n; j++){
    it=0;
    while(it<K){
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        I[k]=i;
        r = jj;
        //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
      }
      I[0]=r;
      jt=0;
      for(i=J+1;i<d;i++) jt+=I[i]*km[i];
      jt+=K*j;
      jt/=(m[J]+1);
      jt=it+jt+K*j;
      dBta[it+K*j]=((I[J]+1)*dBta[jt+km[J]]+(m[J]-I[J]+1)*dBta[jt])/(double)(m[J]+2);
      it++;
    }
  }
  Free(I);
}
/*////////////////////////////////////////////////////////////*/
/*     Initial Guess of p for EM iteration when               */
/*     only the J-th degree changes by diff=1 or -1           */
/*                                                            */
/*////////////////////////////////////////////////////////////*/
void Init_p(double *p, int *m, int d, int J, int diff, int *km){
  int i, jj, k, r, K, it, jt, *I;
  K=km[d];
  I=Calloc(d, int);
  if(diff==1){
    it=K-1;
    while(it>=0){
      r = it;
      for(k=d-1; k>0; k--){
          jj = r%km[k];
          i = (r-jj)/km[k];
          I[k]=i;
          r = jj;
          //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
      }
      I[0]=r;
      jt=0;
      for(i=J;i<d;i++) jt+=I[i]*km[i];
      jt/=(m[J]+1);
      jt=it-jt;
      if(I[J]==m[J]) 
          p[it]=m[J]*p[jt-km[J]]/(double)(m[J]+1.0);
      else if(I[J]>0) p[it]=(I[J]*p[jt-km[J]]+(m[J]-I[J])*p[jt])/(double)(m[J]+1.0);
      else p[it]=m[J]*p[jt]/(double)(m[J]+1.0);
      it--;
    }
  }
  else{
    double w=0.0;
    it=0;
    while(it<K){
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        I[k]=i;
        r = jj;
        //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
      }
      I[0]=r;
      jt=0;
      for(i=J;i<d;i++) jt+=I[i]*km[i];
      jt/=(m[J]+1);
      jt=it+jt;
      if(I[J]==0) 
        p[it]=p[jt+km[J]]+p[jt]/(double)(m[J]+1.0);
      else if(I[J]<m[J]) 
        p[it]=p[jt+km[J]]/(double)I[J]+p[jt]/(double)(m[J]-I[J]+1.0);
      else p[it]=p[jt+km[J]]/(double)(m[J]+1.0)+p[jt];
      w+=p[it];
      it++;
    }
    for(it=0;it<K;it++) p[it]/=w;
  }
  for(it=0;it<K;it++) p[it]=(p[it]+1.0/(double)(K*K))/(1.0/(double)K);
  Free(I);
}
// select degrees by minimizing distance between marginal densities
// estimated based on marginal data and joint data
void mable_mvar(int *M0, int *M, int *n, int *d, int *search, double *phat, int *mhat, 
      double *x, int *maxit,  double *eps, double *lk, int *progress,   
      int *conv, double *D, double *mlk, int *cdf){
  int i, j, k, l, prgrs, *km, *m, *mt, K, maxM;
  int N, *kN, itmp, *Itmp0, *Itmp1, jtmp=0, diff=0, Kt=0; 
  double *p, *pt, *pt1, *pt2, *llik, pct=0.0, ttl,  *Bta;
  double maxD, minimaxD=1e+5;//, Dt, tini=1e-5;
  km = Calloc(*d+1, int);
  kN = Calloc(*d+1, int);
  mt = Calloc(*d, int);
   m = Calloc(*d, int);
  Itmp0 = Calloc(*d, int);
  Itmp1 = Calloc(*d, int);
  llik = Calloc(1, double);
  km[0] = 1;
  maxM=0;
  // Rprintf("km[%d]=%d\n",0,km[0]);
  for(i=1; i<=*d; i++){
    km[i]=km[i-1]*(M[i-1]+1);
    maxM = imax2(maxM, M[i-1]);
    // Rprintf("km[%d]=%d\n",i,km[i]);
  }
  K=km[*d];
  Bta = Calloc(*n*K, double);
  p = Calloc(K, double);
  pt = Calloc(K, double);
  pt1 = Calloc(maxM+1, double);
  pt2 = Calloc(maxM+1, double);
  if(*progress==1)
    Rprintf("\n Mable fit of multivariate data. This may take several minutes.\n\n");
  if(*search==0){
    if(*progress==1) prgrs=1;
    else prgrs=0;
    Multivar_dBeta(x, M, *n, *d, km, Bta);
    for(i=0;i<K;i++) pt[i]=1.0/(double) K;
    em_multivar_beta_mix(pt, Bta, M, *n, *d, K, *maxit, *eps, lk, prgrs, conv);
    k = 0;
    maxD=0.0;
    for(i=0; i<*d; i++){
      for(j=0;j<=mhat[i];j++) pt1[j]=phat[k+j];
      k+=mhat[i]+1;
      p2pj(pt, M, *d, km, pt2, i);
      if(*cdf) maxD=fmax2(maxD, L2_F1F2(pt1, mhat[i], pt2, M[i]));
      else maxD=fmax2(maxD, L2_f1f2(pt1, mhat[i], pt2, M[i]));
      mt[i]=M[i];
    }
    minimaxD=maxD;
    mlk[0]=lk[0];
    Kt=K;
  }
  else{
    // Start with m=M0
    prgrs=0;
    kN[0]=1;
    km[0] = 1;
    for(i=0; i<*d; i++){
      M[i]=M[i]-M0[i];// candidate deggrees: M0[i]+0,...,M0[i]+M[i]
      kN[i+1]=kN[i]*(M[i]+1);
      Itmp0[i]=0;
      m[i] = M0[i];
      km[i+1]=km[i]*(m[i]+1);
    }
    N=kN[*d];// total possible combinations of degrees
    K = km[*d];
    ttl=N;
    for(j=0; j<K; j++) p[j]=1.0/(double) K;
    Multivar_dBeta(x, m, *n, *d, km, Bta);
      em_multivar_beta_mix(p, Bta, m, *n, *d, K, *maxit, *eps, llik, prgrs, conv);
    l = 1;
    while(l<N){
      k=0;
      for(i=0;i<*d;i++) k+=Itmp0[i]*kN[i]; 
      lk[k] = llik[0];// column-major order
      k = 0;
      maxD=0.0;
      for(i=0; i<*d; i++){
        for(j=0;j<=mhat[i];j++) pt1[j]=phat[k+j];
        k+=mhat[i]+1;
        p2pj(p, m, *d, km, pt2, i);
        if(*cdf) maxD=fmax2(maxD, L2_F1F2(pt1, mhat[i], pt2, m[i]));
        else maxD=fmax2(maxD, L2_f1f2(pt1, mhat[i], pt2, m[i]));
      }
      if(maxD<minimaxD){
        minimaxD=maxD;
        for(i=0; i<*d; i++) mt[i]=m[i];
        for(i=0; i<K; i++) pt[i]=p[i];
        mlk[0]=lk[k];
        Kt=K;
        //Rprintf("l=%d, minimaxD=%g\n",l,minimaxD);
      }
      //Rprintf("l=%d, minimaxD=%g\n",l,minimaxD);
      itmp = l;
      j = *d-1;
      while(j>=0){
        Itmp1[j] = itmp/kN[j];
        itmp = itmp-Itmp1[j]*kN[j];
        if(Itmp1[j]%2==1)
          itmp = kN[j]-itmp-1;
        if(Itmp1[j]-Itmp0[j]!=0){
          diff=Itmp1[j]-Itmp0[j];
          for(i=j+1; i<=*d;i++) km[i]=km[i]*(m[j]+1+diff)/(m[j]+1); 
          m[j] = m[j]+diff;                      
          jtmp=j;
        }
        j--;
      }
      for(i=0;i<*d;i++) Itmp0[i]=Itmp1[i];
      K=km[*d];
      /*
      Rprintf("m=(");
      for(i=0;i<*d;i++){
        Rprintf("%d", m[i]);
        if(i<*d-1)Rprintf(", ");
      }
      Rprintf("), ");
      Rprintf("D=%1.4f, m=(", maxD);
      for(i=0;i<*d;i++){ 
        Rprintf("%d", mt[i]);
        if(i<*d-1)Rprintf(", ");
      }
      Rprintf("), minD=%1.4f.\n", minimaxD);
      */
      Update_dBeta(x, m, *n, *d, jtmp, diff, km, Bta);   
      //Multivar_dBeta(x, m, *n, *d, km, Bta);
      if(diff==1) Init_p(p, m, *d, jtmp, diff, km);
      else for(j=0; j<K; j++) p[j]=1.0/(double) K;
      em_multivar_beta_mix(p, Bta, m, *n, *d, K, *maxit, *eps, llik, prgrs, conv);
      l++;
      pct += 1.0;
      if(*progress==1) ProgressBar(pct/ttl,"");
      R_CheckUserInterrupt();
    }
    if(*progress==1){
      ProgressBar(1.0,"");
      Rprintf("\n");
    }
  }
  for(j=0;j<*d;j++) mhat[j]=mt[j];
  for(i=0;i<Kt;i++) phat[i]=pt[i];
  D[0]=minimaxD;
  Free(Bta);
  Free(p);
  Free(pt);
  Free(pt1);
  Free(pt2);
  Free(km);
  Free(kN);
  Free(m);
  Free(mt);
  Free(Itmp0);
  Free(Itmp1);
  Free(llik);
}
/* t=(t1,...,td) is in [0,1]^d */
void mable_mvdf(int *d, int *m, int *km, int *n, double *t, double *p, 
        double *mvdf, int *density){
  int i, j, K;
  double *tmp;
  K=km[*d];
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
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*              C Program for EM-Algorithm                */
/*   Maximum Approximate Bernstein likelihood Density     */
/*    Estimation under two-sample density ratio model     */
/*                 r(x)=(1,x,...,x^d)                     */
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/*                                                            */
/*  Reference:                                                */
/*   Zhong Guan, Maximum Approximate Bernstein Likelihood     */
/*        Estimations of Density and ROC Curve Based on       */
/*        Grouped or Ungrouped Continuous Data                */
/*           Under a Density Ratio Model                      */
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
/*/////////////////////////////////////////////////////////////////*/
/*   MABLE for Density Ratio Model: f1(x)=f0(x)exp[alpha'r(x)]     */
/*    Maximum Approximate Bernstein Likelihood Estimation of       */
/*  density f0 and regression coefficients alpha base twe samples  */
/*  x0: n0-vector, sample from distribution f0 with support [0,1]  */
/*  x1: n1-vector, sample from distribution f1 with support [0,1]  */
/*/////////////////////////////////////////////////////////////////*/
/*/////////////////////////////////////////////////////////////////*/
/*              Log-Likelihood ell(alpha, p), where                */
/*        alpha=(alpha0,...,alpha_d), p=(p0, p1, ..., pm),         */
/*       z=(z1,...,zn)=(x01,...,x0n0,x11,...,x1n1), n=n0+n1,       */
/*                   rx1 = [r(x11),...,r(x1n1)]                    */
/*/////////////////////////////////////////////////////////////////*/
//loglik_bern(alpha, p, ry, beta_x, beta_y, m, nx, ny, d);
double loglik_bern(double *alpha, double *p, double *rx1, double *beta_x0, double *beta_x1, 
        int m, int n0, int n1, int d){
    int i,j;
    double llkhd, fx;
    llkhd = 0.0;
    for(i=0; i<n0; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*beta_x0[i+n0*j];
        }
        llkhd += log(fx);
    }
    for(i=0; i<n1; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*beta_x1[i+n1*j];
        }
        llkhd += log(fx);
    }
    for(i=0; i<n1; i++){
        fx = 0.0;
        for(j=0; j<=d; j++){
            fx += alpha[j]*rx1[i+n1*j];
        }
        llkhd += fx;
    }
    return(llkhd);
}
/*/////////////////////////////////////////////////////////////////*/
/*    Log-Likelihood for grouped data density ratio model          */
/*/////////////////////////////////////////////////////////////////*/
//loglik_bern_group(p, N, n0, n1, Pm, Pm_alfa, m, d);
double loglik_bern_group(double *p, int N, int *n0, int *n1,
            double *Pm, double *Pm_alfa, int m, int d){
    int i,j;
    double llkhd, fx;
    llkhd = 0.0;
    for(i=0; i<N; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*Pm[i+N*j];
        }
        llkhd += n0[i]*log(fx);
    }
    for(i=0; i<N; i++){
        fx = 0.0;
        for(j=0; j<=m; j++){
            fx += p[j]*Pm_alfa[i+N*j];
        }
        llkhd += n1[i]*log(fx);
    }
    return(llkhd);
}
/*////////////////////////////////////////////////////////////////////*/
/*   Integrand  r_j(x)*r_k(x)*beta_{mi}(x)*exp(alpha'r(x)), where     */
/*        r(x)=[1,r_1(x),...,r_d(x)], j,k=0,...,d, i=0,...,m.         */
/* ex = (m, i, j, k, d, a, b, alpha)                                  */
/*////////////////////////////////////////////////////////////////////*/
/* called via .External(.) :*/
//SEXP C_mable_dr(SEXP args);

//typedef struct int_struct
//{
//    SEXP f;    /* function */
//    SEXP env;  /* where to evaluate the calls */
//} int_struct, *IntStruct;

typedef struct mable_dr_struct
{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
    int m; int i; int j; int k; int d;
    double *alpha;
} mable_dr_struct, *MableDRStruct;
void func_ebeta_rjk(double *x, int n, void *ex)
{
    SEXP args, res_sxp, tmp;
    int i, j, k, m, d, ii, jj;
    double texp, *alpha;
    MableDRStruct MDS = (MableDRStruct) ex;
    m = MDS->m; i = MDS->i; j = MDS->j; k = MDS->k; 
    d = MDS->d; alpha = MDS->alpha;
    PROTECT(args = allocVector(REALSXP, n));
    for(ii = 0; ii < n; ii++) REAL(args)[ii] = x[ii];

    PROTECT(tmp = lang2(MDS->f, args));
    PROTECT(res_sxp = eval(tmp, MDS->env));

    if(length(res_sxp) != n*(d+1))
	error("evaluation of regression function(s) gave a result of wrong length");
    if(TYPEOF(res_sxp) == INTSXP) {
	res_sxp = coerceVector(res_sxp, REALSXP);
    } else if(TYPEOF(res_sxp) != REALSXP)
	error("evaluation of regression function(s) gave a result of wrong type");

    for(ii=0;ii<n;ii++) { 
        texp  = 0.0;
        for(jj=0;jj<=d;jj++) texp+=REAL(res_sxp)[ii+n*jj]*alpha[jj];
        x[ii] = REAL(res_sxp)[ii+n*j]*REAL(res_sxp)[ii+n*k]*dbeta(x[ii],i+1,m-i+1,FALSE)*exp(texp);
	    if(!R_FINITE(x[ii]))
	       error("non-finite r(x) value");
    }
    UNPROTECT(3);
    return;
}

/*/////////////////////////////////////////////////////////////////////////////*/
/*      Returns case cdf components  n x (m+1) matrix with (i,j)-element       */
/*  \int_0^z beta_{mj}(z[i])*exp(alpha[0]+alpha[1]*r1(x)+...+alpha[d]*rd(x)),  */
/*                      i=1,...,n; j=0,...,m, n=length(z)                      */
/*/////////////////////////////////////////////////////////////////////////////*/
// z =(z1,...,z_n)
void Betam_alpha(double *alpha, double *z, double *beta_alpha, int d, int m, int n, void *ex){
    mable_dr_struct mds;
    IntStruct IS = (IntStruct) ex;
    int i, j;
    double l=.0, u, epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400]; 
    int lenw=400, last=0, neval=0, ier=0, iwork[100]; 
    int limit=100; 
    mds.f = IS->f;
    mds.env = IS->env;
    mds.m = m; mds.d = d; 
    mds.j = 0; mds.k = 0;
    mds.alpha = alpha;
    for(j=0;j<n;j++){
       u = z[j];
       for(i=0; i<=m; i++){
          mds.i = i; 
          Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
          beta_alpha[j+i*n] = result;
       }
    }
}
/*///////////////////////////////////////////////////////////////////////*/
/*  Group prob:  P[t_{j-1}<X1<t_j | X1 ~ beta_{mi}(x)exp(alpha'r(x))]    */
/*       P(Ij; alpha)=(Pm0(Ij;alpha),...,Pmm(Ij;alpha)), j=1:N           */
/*   Pmi(Ij;alpha)=int_{t_{j-1}}^{t_j} beta_{mi}(x) exp(alpha'r(x))dx    */
/*     and derivatives w.r.t. alpha,  r(x)=(1, r1(x),...,rd(x))          */
/*///////////////////////////////////////////////////////////////////////*/
//  dPm_alpha(alpha, t, N, d, m, Pm, dPm, ddPm, ex);
// Returns
//   Pm: N x (m+1) matrix of P_{mi}(I_j;alpha): j=1:N, i=0:m
//  dPm: N x (m+1) x (d+1) array and 
// ddPm: N x (m+1) x (d+1) x (d+1) array
void dPm_alpha(double *alpha, double *t, int N, int d, int m, 
            double *Pm, double *dPm, double *ddPm, void *ex){
    mable_dr_struct mds;
    IntStruct IS = (IntStruct) ex;
    int i, j, k, ii, jj, Nm=N*(m+1), Nmd=Nm*(d+1);
    double l, u, epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400];
    int lenw=400, last=0, neval=0, ier=0, iwork[100];
    int limit=100;
    mds.f = IS->f;
    mds.env = IS->env;
    mds.m = m; mds.d = d; 
    mds.alpha = alpha;
    for(ii = 0; ii < N; ii++){
        l = t[ii]; u = t[ii+1];  
        for(i=0; i<=m; i++){
            mds.i = i; mds.j = 0; mds.k = 0;
            Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
            jj = ii+N*i;
            Pm[jj] = result;
            dPm[jj] = result;
            ddPm[jj] = result;
            for(j=1; j<=d; j++) {
                mds.j = j; mds.k = 0;
                Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
                dPm[jj+Nm*j] = result;
                ddPm[jj+Nm*j] = result;
                ddPm[jj+Nmd*j] = result;
                for(k=j; k<=d; k++){
                    mds.k = k;
                    Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
                    ddPm[jj+Nm*j+Nmd*k] = result;
                    ddPm[jj+Nm*k+Nmd*j] = result;
                }
            }
        }
    }
}
// no derivatitive
void Pm_alpha(double *alpha, double *t, int N, int d, int m, double *Pm,  void *ex){
    mable_dr_struct mds;
    IntStruct IS = (IntStruct) ex;
    int i, ii;
    double l, u, epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400];
    int lenw=400, last=0, neval=0, ier=0, iwork[100];
    int limit=100;
    mds.f = IS->f;
    mds.env = IS->env;
    mds.m = m; mds.d = d; 
    mds.j = 0; mds.k = 0;
    mds.alpha = alpha;
    for(ii = 0; ii < N; ii++){
        l = t[ii]; u = t[ii+1];  
        for(i=0; i<=m; i++){
            mds.i = i; 
            Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
            Pm[ii+N*i] = result;
        }
    }
}

/*///////////////////////////////////////////////////*/
/*  cdf of mixure of exponentially tilted betas      */
/*///////////////////////////////////////////////////*/
SEXP  mixtbeta_cdf(SEXP args){
    int i, j, d, m, n;
    int_struct is;
    SEXP ans, ansnames;
    double *alpha, *p, *x, *y, *beta_alpha, tmp;

    // *Read-only* R objects >
    args = CDR(args);
    is.f = CAR(args); args = CDR(args);
    is.env = CAR(args); args = CDR(args);
    alpha = REAL(CAR(args)); args = CDR(args);
    p = REAL(CAR(args)); args = CDR(args);
    x = REAL(CAR(args)); args = CDR(args);
    d = asInteger(CAR(args)); args = CDR(args);
    m = asInteger(CAR(args)); args = CDR(args);
    n = asInteger(CAR(args)); args = CDR(args);
    // *Read-only* R objects ^
    beta_alpha = Calloc(n*(m+1), double);
    y = Calloc(n, double);

    Betam_alpha(alpha, x, beta_alpha, d, m, n, (void*)&is);
    for(j=0;j<n;j++) {
        tmp = 0.0;
        for(i=0;i<=m;i++) tmp+=p[i]*beta_alpha[j+i*n];
        y[j] = tmp;
    }
    PROTECT(ans = allocVector(VECSXP, 2));
    PROTECT(ansnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(ansnames, 0, mkChar("x"));
    SET_STRING_ELT(ansnames, 1, mkChar("y"));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, n));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, n));
    for(i=0;i<n;i++){
        REAL(VECTOR_ELT(ans, 0))[i] = x[i];
        REAL(VECTOR_ELT(ans, 1))[i] = y[i];
    }
    Free(beta_alpha); Free(y); 
    setAttrib(ans, R_NamesSymbol, ansnames);
    UNPROTECT(2);
    return ans;
}

/*////////////////////////////////////////////////////*/
/* weight function w(alpha)=(w0(alpha),...,wm(alpha)) */
/*  wi(alpha)=int_0^1 beta_{mi}(x) exp(alpha'r(x))dx  */
/*               r(x)=(1, r1(x),...,rd(x))            */
/*////////////////////////////////////////////////////*/  
//    wt_alpha(alpha, d, m, wt, ex);
void wt_alpha(double *alpha, int d, int m, double *wt, void *ex){
    int i;
    mable_dr_struct mds;
    IntStruct IS = (IntStruct) ex;
    double l=.0, u=1., epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400]; 
    int lenw=400, last=0, neval=0, ier=0, iwork[100]; 
    int limit=100; 
    mds.f = IS->f;
    mds.env = IS->env;
    mds.m = m; mds.d = d; 
    mds.j = 0; mds.k = 0;
    mds.alpha = alpha;
    for(i=0; i<=m; i++){
       mds.i = i;
       Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
       wt[i] = result;
    }
}
/*/////////////////////////////////////*/
/*    Log-Likelihood when m=0          */
/*/////////////////////////////////////*/
double loglik_unif(double *alpha, double *rx1, int n1, int d){
  int i, j;
  double lk;
  lk=0.0;
  for(i=0; i<n1; i++){
    for(j=0; j<=d; j++){
      lk += alpha[j]*rx1[i+n1*j];
    }
  }
  return(lk);
}
double loglik_unif_group(double *t, int N, int *n0, int *n1, double *alpha, 
        int m, int d, void *ex){
    int i;
    double llkhd, fx, *Pm, *Pm_alfa;
    Pm = Calloc(N, double); 
    Pm_alfa = Calloc(N, double);
    llkhd = 0.0;
    cpBeta(t, 0, N, Pm);  
    Pm_alpha(alpha, t, N, d, 0, Pm_alfa, ex);

    for(i=0; i<N; i++){
        fx = Pm[i];
        llkhd += n0[i]*log(fx);
    }
    for(i=0; i<N; i++){
        fx = Pm_alfa[i];
        llkhd += n1[i]*log(fx);
    }
    Free(Pm); Free(Pm_alfa);
    return(llkhd);
}

/*////////////////////////////////////////////////////*/
/* weight function w(alpha)=(w0(alpha),...,wm(alpha)) */
/*  wi(alpha)=int_0^1 beta_{mi}(x) exp(alpha'r(x))dx  */
/*  and its derivatives w.r.t. alpha                  */
/*   r(x)=(1, r_1(x),...,r_d(x))                      */
/*////////////////////////////////////////////////////*/  
//    weights(alpha, d, m, wt, dwt, ddwt, ex);
void weights(double *alpha, int d, int m, 
        double *wt, double *dwt, double *ddwt, void *ex){
    int i, j, k, d1=d+1, m1=m+1, md=m1*d1;
    mable_dr_struct mds;
    IntStruct IS = (IntStruct) ex;
    double l=.0, u=1., epsabs=.00001, epsrel=.00001;
    double result=0.0, abserr=0.0, work[400]; 
    int lenw=400, last=0, neval=0, ier=0, iwork[100]; 
    int limit=100; 
    mds.f = IS->f;
    mds.env = IS->env;
    mds.m = m; mds.d = d; 
    mds.alpha = alpha;
    for(i=0; i<=m; i++){
        mds.i = i; mds.j = 0; mds.k = 0; 
        Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
        wt[i] = result;
        dwt[i] = result;
        ddwt[i] = result;
        for(j=1; j<=d;j++) {
            mds.j = j; mds.k = 0; 
            Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
            dwt[i+m1*j] = result;
            ddwt[i+m1*j] = result;
            ddwt[i+md*j] = result;
            for(k=j; k<=d; k++){
                mds.k = k;   
                Rdqags(func_ebeta_rjk, (void*)&mds, &l, &u, &epsabs, &epsrel, &result, &abserr, &neval, &ier, &limit, &lenw, &last, iwork,  work);
                ddwt[i+m1*j+md*k] = result;
                ddwt[i+m1*k+md*j] = result;
          }  
       }
    }
}
/*//////////////////////////////////////////////////*/
/*           Scores and Jacobian matrix             */
/*//////////////////////////////////////////////////*/       
//        score_alpha(alpha, ny, d, m, n, ry, Tk, H, Jac, wt, ex);
void score_alpha(double *alpha, int ny, int d, int m, int n, double *ry, 
           double *Tk, double *H, double *Jac, double *wt, void *ex){
    int i, j, k, nx=n-ny, d1=d+1, m1=m+1, md=m1*d1, mdd=md*d1;
    double *dwt, *ddwt; 
    dwt = Calloc(md, double);       
    ddwt =  Calloc(mdd, double);
    weights(alpha, d, m, wt, dwt, ddwt, ex);
    for(j=0;j<=d;j++){
        H[j] = 0.0;
        for(i=0;i<=m;i++) H[j] -= dwt[i+m1*j]*Tk[i]/(nx+ny*wt[i]);
        H[j] *= ny;
        for(i=0;i<ny;i++) {
            H[j] += ry[i+ny*j];
        }
        for(k=0;k<=d;k++){
           Jac[j+d1*k]=0.0;
           for(i=0;i<=m;i++)
              Jac[j+d1*k] += (ddwt[i+m1*j+md*k]*(nx+ny*wt[i])-ny*dwt[i+m1*j]*dwt[i+m1*k])*Tk[i]/((nx+ny*wt[i])*(nx+ny*wt[i]));
           Jac[j+d1*k] *= -ny;
        }  
    }
    Free(dwt); Free(ddwt);   
}
/*/////////////////////////////////////////*/
/*   Scores and Jacobian matrix of alpha   */
/*  for grouped data density ratio model   */
/*/////////////////////////////////////////*/
// T_zero: (m+1)-vector of T0[j]=sum{n0[i]*p[j]*Pm[i,j]/sum(p*Pm[i,]): i=1,...,N}, j=0,...,m
//  T_alpha: (m+1)-vector of T1[j]=sum{n1[i]*p[j]*Pm_alfa[i,j]/sum(p*Pm_alfa[i,]): i=1,...,N}, j=0,...,m
//  Tk=T_zero + T_alpha; 
// Pm[i,j]=Pm_alfa[i,j], with alpha=0, N x (m+1) matrix
// Pi_theta: N x (m+1) matrix of n1[i]*p[j]*Pm_alfa[i,j]/sum(p*Pm_alfa[i,]): i=1,...,N}, j=0,...,m
//   score_alpha_group(alpha, t, N, nx, ny, d, m, Pi_theta, T_zero, Tk, H, Jac, wt, ex);
void score_alpha_group(double *alpha, double *t, int N, int nx, int ny, int d, int m, 
        double *Pi_theta, double *T_zero, double *Tk, double *H, double *Jac, double *wt, void *ex){
    int i, j, k, l, d1=d+1, m1=m+1, md=m1*d1, Nm=N*m1, Nmd=Nm*d1;
    double *Pm_alfa, *dPm_alfa, *ddPm_alfa, *dwt, *ddwt; 
    Pm_alfa = Calloc(Nm, double);
    dPm_alfa = Calloc(Nmd, double);
    ddPm_alfa = Calloc(Nmd*d1, double);
    dwt = Calloc(md, double);
    ddwt =Calloc(md*d1, double);
    dPm_alpha(alpha, t, N, d, m, Pm_alfa, dPm_alfa, ddPm_alfa, ex);
    for(i=0;i<=m;i++){
        wt[i]=0.0; Tk[i]=T_zero[i];
        for(k=0;k<N;k++) {
           wt[i]+=Pm_alfa[k+N*i]; Tk[i]+=Pi_theta[k+N*i];
        }
        for(j=0;j<=d;j++){
           dwt[i+m1*j]=0.0;
           for(k=0;k<N;k++) dwt[i+m1*j]+=dPm_alfa[k+N*i+Nm*j];
           for(l=0;l<=d;l++){
               ddwt[i+m1*j+md*l]=0.0;
               for(k=0;k<N;k++) 
                  ddwt[i+m1*j+md*l]+=ddPm_alfa[k+N*i+Nm*j+Nmd*l];
           }
        }
    }
    for(j=0;j<=d;j++){
        H[j] = 0.0;
        for(i=0;i<=m;i++) H[j] -= dwt[i+m1*j]*Tk[i]/(nx+ny*wt[i]);
        H[j] *= ny;
        for(i=0;i<=m;i++)
           for(k=0;k<N;k++) H[j] += Pi_theta[k+N*i]*dPm_alfa[k+N*i+Nm*j]/Pm_alfa[k+N*i];
        for(k=0;k<=d;k++){
           Jac[j+d1*k]=0.0;
           for(i=0;i<=m;i++)
              Jac[j+d1*k] -= (ddwt[i+m1*j+md*k]*(nx+ny*wt[i])-ny*dwt[i+m1*j]*dwt[i+m1*k])*Tk[i]/((nx+ny*wt[i])*(nx+ny*wt[i]));
           Jac[j+d1*k] *= ny;
           for(i=0;i<=m;i++)
              for(l=0;l<N;l++) 
                Jac[j+d1*k] += Pi_theta[l+N*i]*(Pm_alfa[l+N*i]*ddPm_alfa[l+N*i+Nm*j+Nmd*k]-dPm_alfa[l+N*i+Nm*j]*dPm_alfa[l+N*i+Nm*k])/(Pm_alfa[l+N*i]*Pm_alfa[l+N*i]);
        }
    }
    Free(Pm_alfa); Free(dPm_alfa); Free(ddPm_alfa); Free(dwt); Free(ddwt);
}

// To be deleted
void checking_of_p(double *p, int m){
    int j;
    double sump=.0;
    for(j=0; j<=m; j++) sump+=p[j];
    Rprintf("sum of p = %g\n", sump);
    for(j=0; j<=m; j++) {
        if(p[j]<0 || p[j]>1)
	       error("wrong p[j] value");
        Rprintf("p[%d] = %g, ",j, p[j]);
    }
    Rprintf("\n");
}
/*//////////////////////////////////////////////////*/
/*  EM Method for profiling likelihood ell(alpha)   */
/*    finding the maximizer p of ell(alpha, p)      */
/*        where alpha is an MELE of alpha           */
/*//////////////////////////////////////////////////*/
/////////////////////////////////////////////////////////////////////////
//em_ptilde_dr(llik, alpha, p, x, y, hy, m, d, nx, ny, a, b, wt, 
//                    veps, eps_em, maxit_em, eps_nt, maxit_nt, ini_item)
// alpha: MELE of alpha: alpha-tilde
//     p: initial value for p
/////////////////////////////////////////////////////////////////////////
void em_ptilde_dr(double *llik, double *alpha, double *p, double *x, double *y, double *ry, 
         int m, int d, int nx, int ny, double *wt, //double veps,
         double eps_em, int maxit_em, double eps_nt, int maxit_nt, void *ex){
    int i, j, n, it_em, it_nt;
    double *beta_x, *beta_y, *fm, *Tk, *H, *Jac, llik_new, sump;  
    double del_em, del_nt, lam, tmp;   
    n = nx+ny;
    beta_x = Calloc(nx*(m+1), double);
    beta_y = Calloc(ny*(m+1), double);
    H = Calloc(d+1, double);
    Tk = Calloc(m+1, double);
    fm = Calloc(n, double);
    Jac = Calloc((d+1)*(d+1), double);
    dBeta(x, m, nx, beta_x);
    dBeta(y, m, ny, beta_y);
    wt_alpha(alpha, d, m, wt, ex);
    del_em=10.0;
    it_em=1;
    llik[0] = loglik_bern(alpha, p, ry, beta_x, beta_y, m, nx, ny, d);
    //EM iteration
    while(del_em>eps_em && it_em<maxit_em){
      R_CheckUserInterrupt();
        for(i=0;i<nx;i++){
           fm[i] = 0.0;
           for(j=0;j<=m;j++)  fm[i] += beta_x[i+nx*j]*p[j];
        }
        for(i=0;i<ny;i++){
           fm[nx+i] = 0.0;
           for(j=0;j<=m;j++)  fm[nx+i] += beta_y[i+ny*j]*p[j];
        }
        for(j=0;j<=m;j++){
           Tk[j] = 0.0;
           for(i=0; i<nx; i++)  Tk[j] += beta_x[i+nx*j]/fm[i];
           for(i=0; i<ny; i++)  Tk[j] += beta_y[i+ny*j]/fm[nx+i];
           Tk[j] *= p[j];
        }  
        // Newton method for finding lambda starting with lambda0=ny
        lam = ny;
        tmp = 0.0;
        for(j=0;j<=m;j++) tmp += Tk[j]*(wt[j]-1.0)/(n+lam*(wt[j]-1.0));
        del_nt = fabs(tmp);
        it_nt = 0;
        while(del_nt>eps_nt && it_nt<maxit_nt){
            del_nt = 0.0;    
            for(j=0;j<=m;j++) del_nt+= Tk[j]*(wt[j]-1.0)*(wt[j]-1.0)/((n+lam*(wt[j]-1.0))*(n+lam*(wt[j]-1.0))); 
            del_nt = tmp/del_nt;
            lam += del_nt;
            del_nt = fabs(del_nt);
            tmp = 0.0;
            for(j=0;j<=m;j++) tmp += Tk[j]*(wt[j]-1.0)/(n+lam*(wt[j]-1.0));
//            del_nt += fabs(tmp);
            del_nt = fabs(tmp);
            it_nt++;
        }
        for(j=0;j<=m;j++) p[j] = Tk[j]/(n+lam*(wt[j]-1.0));
    sump=.0;
    for(j=0; j<=m; j++) sump+=p[j]*wt[j];
        llik_new = loglik_bern(alpha, p, ry, beta_x, beta_y, m, nx, ny, d);
        del_em = fabs(llik_new-llik[0]);
        llik[0] = llik_new;
        it_em++;
    }
    Free(beta_x); Free(beta_y); Free(fm); Free(H); Free(Jac); Free(Tk); 
}

/*//////////////////////////////////////////////////*/
/*     EM Method for maximizing ell(alpha, p)       */
/*//////////////////////////////////////////////////*/
//  em_mable_dr(llik, alpha, p, x0, x1, rx1, m, d, n0, n1, wt, eps_em, eps, maxit_em, maxit,(void*)&is);
void em_mable_dr(double *llik, double *alpha, double *p, double *x, double *y,  
         double *ry, int m, int d, int nx, int ny, double *wt, //double veps,
         double eps_em, double eps_nt, int maxit_em, int maxit_nt,  
         void *ex, int progress){
    int i, j, n, it_nt, it_em;
    double *beta_x, *beta_y, *fm, *H, *Jac, *Tk, llik_new;  
    double del_em, del_nt, *tmp, sump; 
    n = nx+ny;
    beta_x = Calloc(nx*(m+1), double);
    beta_y = Calloc(ny*(m+1), double);
    H = Calloc(d+1, double);
    Tk = Calloc(m+1, double);
    fm = Calloc(n, double);
    Jac = Calloc((d+1)*(d+1), double);
    tmp = Calloc(d+1, double);
    // profiling likelihood at alpha
    dBeta(x, m, nx, beta_x);
    dBeta(y, m, ny, beta_y);
    del_em=10.0;
    it_em=1;
    llik[0] = loglik_bern(alpha, p, ry, beta_x, beta_y, m, nx, ny, d);
    //EM iteration
    while(del_em>eps_em && it_em<maxit_em){
      R_CheckUserInterrupt();
        for(i=0;i<nx;i++){
           fm[i] = 0.0;
           for(j=0;j<=m;j++)  fm[i] += beta_x[i+nx*j]*p[j];
        }
        for(i=0;i<ny;i++){
           fm[nx+i] = 0.0;
           for(j=0;j<=m;j++)  fm[nx+i] += beta_y[i+ny*j]*p[j];
        }
        for(j=0;j<=m;j++){
           Tk[j] = 0.0;
           for(i=0; i<nx; i++)  Tk[j] += beta_x[i+nx*j]/fm[i];
           for(i=0; i<ny; i++)  Tk[j] += beta_y[i+ny*j]/fm[nx+i];
           Tk[j] *= p[j];
        } 
        score_alpha(alpha, ny, d, m, n, ry, Tk, H, Jac, wt, ex);
        del_nt=0.0;
        for(i=0;i<=d;i++) del_nt += fabs(H[i]);
        it_nt=1;
        //Newton iteration:
        while(del_nt>eps_nt && it_nt<maxit_nt){
            minverse(Jac, d+1);  
            for(i=0;i<=d;i++){
               tmp[i] = 0.0;
               for(j=0;j<=d;j++) tmp[i] += Jac[i+(d+1)*j]*H[j];
            }
            del_nt = 0.0;
            for(i=0;i<=d;i++){
               alpha[i] -= tmp[i];
               del_nt += fabs(tmp[i]);
            }
            score_alpha(alpha, ny, d, m, n, ry, Tk, H, Jac, wt, ex);
            for(i=0;i<=d;i++) del_nt += fabs(H[i]);
            //if(progress==1) 
            //  Rprintf("m = %d,  Newton: it=%d, del=%g\n", m, it_nt, del_nt);
            it_nt++;
        }
        for(j=0;j<=m;j++) p[j] = Tk[j]/(nx+ny*wt[j]);
        llik_new = loglik_bern(alpha, p, ry, beta_x, beta_y, m, nx, ny, d);
        del_em = fabs(llik_new-llik[0]);
        llik[0] = llik_new;
        it_em++;
        if(progress==1) clockProgress((int)(it_em/70)," EM iteration is runing..."); 
        //  Rprintf("m = %d, EM: it_em=%d, del_em=%g\n", m, it_em, del_em);
    }
    if(progress==1) Rprintf("\n"); 
    sump=.0;
    for(j=0; j<=m; j++) sump+=p[j];
    sump=.0;
    for(j=0; j<=m; j++) sump+=p[j]*wt[j];
    Free(beta_x); Free(beta_y); Free(fm); Free(Jac); Free(Tk); Free(H); 
    Free(tmp);  
}
/*////////////////////////////////////////////*/
/*  Optimal degree m by change-point mehtod   */
/*      based on ell(alpha-hat, p-hat)        */
/*////////////////////////////////////////////*/
SEXP  C_mable_dr(SEXP args){
  //if(progress==1) 
  //  Rprintf("\n Program 'mable.dr' is runing. This may take several minutes.\n\n\n");
  int_struct is;
  int i, j, m, k, d, n0, n1, optim, *M, maxit, maxit_em, progress, emprog=0;
  int *chpts, cp0=0, cp1=1,*cp, itmp, lp, message, vb, i0=0, i1=0;
  SEXP ans, ansnames;
  double *alpha, *p, *lr, *lk, *x0, *x1, *wt, eps_em, eps, tini, tmp, pct;
  double *z, *rx1, *alpha_hat, *alpha_ini, *phat, *wt_ahat, *llik;
  double *pval, *res, pv0=1.0, pv1=1.0, level;  
  args = CDR(args);
  is.f = CAR(args); args = CDR(args);
  is.env = CAR(args); args = CDR(args);
  x0 = REAL(CAR(args)); args = CDR(args);
  x1 = REAL(CAR(args)); args = CDR(args);
  rx1 = REAL(CAR(args)); args = CDR(args);
  M = INTEGER(CAR(args)); args = CDR(args);
  alpha_ini = REAL(CAR(args)); args = CDR(args);
  d = asInteger(CAR(args)); args = CDR(args);
  n0 = asInteger(CAR(args)); args = CDR(args);
  n1 = asInteger(CAR(args)); args = CDR(args);
  eps_em = asReal(CAR(args)); args = CDR(args);
  eps = asReal(CAR(args)); args = CDR(args);
  maxit_em = asInteger(CAR(args)); args = CDR(args);
  maxit = asInteger(CAR(args)); args = CDR(args);
  tini = asReal(CAR(args)); args = CDR(args);
  progress = asInteger(CAR(args)); args = CDR(args);
  level = asReal(CAR(args)); args = CDR(args);
  message = asInteger(CAR(args)); args = CDR(args);
  vb = asInteger(CAR(args)); args = CDR(args);

  k =M[1]-M[0]; 
  if(progress==1 && k>0)  ProgressBar(0.0,""); 
  emprog=progress*(k==0);
  alpha_hat = Calloc((k+1)*(d+1), double);
  alpha = Calloc(d+1, double);
  lp=M[0]*(k+1)+(k+1)*(k+2)/2;
  phat = Calloc(lp, double);
  pval = Calloc(k+1, double);
  chpts = Calloc(k+1, int);
  p= Calloc(M[1]+1, double);
  wt = Calloc(M[1]+1, double);
  wt_ahat = Calloc(lp, double);
  llik = Calloc(1, double);
  lk = Calloc(k+1, double);
  lr = Calloc(k+1, double);
  z = Calloc(n0+n1, double);
  cp = Calloc(1, int);
  res = Calloc(1, double);

  if(vb==-1 || vb==2) i0=1;
  if(vb== 1 || vb==2) i1=1;
  
  for(i=0;i<n0;i++) z[i] = x0[i];
  for(i=0;i<n1;i++) {
     z[i+n0] = x1[i];
  } 
  lr[0]=.0;
  m = M[0];  
  //for(j=0;j<=m;j++) p[j]=1.0/(double)(m+1.0);
  if(m<=2) for(j=0;j<=m;j++) p[j]=1.0/(double)(m+1.0);
  else{
    for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
    if(vb==-1 || vb==2) p[0]=0.0;
    if(vb== 1 || vb==2) p[m]=0.0;
  }
  pval[0]=1.0;
  chpts[0]=0;
  itmp=m+1; 
  for(i=0; i<=d; i++) alpha[i] = alpha_ini[i];
  tmp =1.0+(double)(k*(k-1));
  pct = level/pval[0];
  if(m>0){
    //if(progress==1)  ProgressBar(fmax2(pct,1.0/tmp),""); 
    em_mable_dr(llik, alpha, p, x0, x1, rx1, m, d, n0, n1, wt, eps_em, eps, 
      maxit_em, maxit, (void*)&is, emprog);
    lk[0] = llik[0];
    for(i=0;i<itmp;i++) {
      phat[i]=p[i];
      wt_ahat[i]=wt[i];
    }       
  }
  else{
    phat[0]=1.0;
    lk[0]=loglik_unif(alpha, rx1, n1, d);
    wt_alpha(alpha, d, m, wt, (void*)&is);
    wt_ahat[0]=wt[0];
  }
  for(i=0;i<=d;i++) alpha_hat[i]=alpha[i];
  if(progress==1 && k>0)  ProgressBar(fmax2(pct,1.0/tmp),""); 
  i=1;
  while(i<=k && pval[i-1]>level){
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
      for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
      if(vb==-1 || vb==2) p[0]=0.0;
      if(vb== 1 || vb==2) p[m]=0.0;
    }
    p[m+1] = (m+1)*p[m]/(double)(m+2);
    for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
    p[0] = (m+1.)*p[0]/(double)(m+2.);
    m = M[0]+i;
    //make sure initial p is in the interior of the simplex
    if(m>3){
      for(j=i0; j<=m-i1; j++) p[j]=(1.0-tini)*p[j]+ tini/(double)(m+1-abs(vb));
    }       
    //for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    for(j=0; j<=d; j++) alpha[j] = alpha_ini[j];
    em_mable_dr(llik, alpha, p, x0, x1, rx1, m, d, n0, n1, wt, eps_em, eps, 
      maxit_em, maxit, (void*)&is, emprog);
    lk[i] = llik[0];
    for(j=0;j<=d;j++) alpha_hat[j+i*(d+1)]=alpha[j];
    for(j=0;j<=m;j++){
      phat[j+itmp]=p[j];
      wt_ahat[j+itmp]=wt[j];
    }
    itmp += m+1;
    if(i>=3){
      cp[0]=i;
      chpt_exp(lk, lr, res, cp);
      pval[i]=res[0];
      chpts[i]=cp[0];
    }
    else{            
      pval[i]=1.0;
      chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
      cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
      cp0=cp1;
      pv0=pv1;
    }
    else pv0=pval[i];
    R_CheckUserInterrupt();
    pct = fmax2(pct, level/pval[i]);
    pct = fmax2((1.0+i*(i+1))/tmp, pct);
    if(progress==1){
        ProgressBar(fmin2(1.0,pct)," ");}
    i++;
  }
  if(m==M[1] && pval[i-1]>level && message){
    //convergence[0]+=1; 
    warning("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value,  %f, of the change-point is returned.\n", pv0);
  }
  M[1]=m;
  itmp=cp0*(M[0]*2+(cp0+1))/2;
  optim=cp0+M[0];
  m = optim;
  k = M[1]-M[0];
  if(progress==1){
    ProgressBar(1.0," ");
    Rprintf("\n");}
  //Rprintf("\n 'mable.dr' done. \n\n\n");
  //if(progress==1) Rprintf("\n 'mable.dr' done. \n\n\n");
  PROTECT(ans = allocVector(VECSXP, 9));
  PROTECT(ansnames = allocVector(STRSXP, 9));
  SET_STRING_ELT(ansnames, 0, mkChar("lk"));
  SET_STRING_ELT(ansnames, 1, mkChar("lr"));
  SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, m+1));
  SET_STRING_ELT(ansnames, 2, mkChar("p"));
  SET_STRING_ELT(ansnames, 3, mkChar("m"));
  SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, 1));
  INTEGER(VECTOR_ELT(ans, 3))[0] = optim;
  SET_VECTOR_ELT(ans, 4, allocVector(REALSXP, m+1));
  SET_STRING_ELT(ansnames, 4, mkChar("wt"));
  for(i=0;i<=m;i++){
      REAL(VECTOR_ELT(ans, 2))[i] = phat[itmp+i];
      REAL(VECTOR_ELT(ans, 4))[i] = wt_ahat[itmp+i];
  }
  SET_VECTOR_ELT(ans, 5, allocVector(REALSXP, d+1));
  SET_STRING_ELT(ansnames, 5, mkChar("alpha"));
  for(i=0;i<=d;i++){
      REAL(VECTOR_ELT(ans, 5))[i] = alpha_hat[i+cp0*(d+1)];
  }
  SET_VECTOR_ELT(ans, 6, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 7, allocVector(INTSXP, k+1));
  SET_STRING_ELT(ansnames, 6, mkChar("pval"));
  SET_STRING_ELT(ansnames, 7, mkChar("chpts"));
  for(i=0;i<=k;i++){
      REAL(VECTOR_ELT(ans, 0))[i] = lk[i];
      REAL(VECTOR_ELT(ans, 1))[i] = lr[i];
      REAL(VECTOR_ELT(ans, 6))[i] = pval[i];
      INTEGER(VECTOR_ELT(ans, 7))[i] = chpts[i]+M[0];
  }
  SET_STRING_ELT(ansnames, 8, mkChar("M"));
  SET_VECTOR_ELT(ans, 8, allocVector(INTSXP, 2));
  for(i=0;i<=1;i++){
      INTEGER(VECTOR_ELT(ans, 8))[i] = M[i];
  }
  Free(alpha_hat); Free(alpha); Free(phat); Free(pval);    
  Free(chpts); Free(p); Free(wt); Free(wt_ahat);  
  Free(llik); Free(lk); Free(lr); Free(z); Free(cp); Free(res);

  setAttrib(ans, R_NamesSymbol, ansnames);
  UNPROTECT(2);
  return ans;
}


/*////////////////////////////////////////////////////*/
/*  Model degree m selection by change-point method   */
/*        based on ell(alpha-tilde, p-tilde)          */
/*////////////////////////////////////////////////////*/
// alpha: MELE of alpha
//     p: initial value for m=m0
SEXP  maple_dr(SEXP args){
  //if(progress==1) 
  //  Rprintf("\n Program 'maple.dr' is runing. This may take several minutes.\n\n\n");
  int_struct is;
  SEXP ans, ansnames;
  int i, j, m, k, d, n0, n1, itmp, *chpts, optim, *M, maxit_nt, maxit_em, progress;
  int cp0=0, cp1=1,*cp, lp, message, vb, i0=0, i1=0;
  double *alpha, *p, *lr, *lk, *x0, *x1, *wt, eps_em, eps_nt, tini, tmp, pct;
  double *z, *rx1, *phat, *wt_ahat, *llik, *pval, *res, level, pv0=1.0, pv1=1.0;  
  args = CDR(args);
  is.f = CAR(args); args = CDR(args);
  is.env = CAR(args); args = CDR(args);
  x0 = REAL(CAR(args)); args = CDR(args);
  x1 = REAL(CAR(args)); args = CDR(args);
  rx1 = REAL(CAR(args)); args = CDR(args);
  M = INTEGER(CAR(args)); args = CDR(args);
  alpha = REAL(CAR(args)); args = CDR(args);
//  p = REAL(CAR(args)); args = CDR(args);
  d = asInteger(CAR(args)); args = CDR(args);
  n0 = asInteger(CAR(args)); args = CDR(args);
  n1 = asInteger(CAR(args)); args = CDR(args);
//  veps = asReal(CAR(args)); args = CDR(args);
  eps_em = asReal(CAR(args)); args = CDR(args);
  eps_nt = asReal(CAR(args)); args = CDR(args);
  maxit_em = asInteger(CAR(args)); args = CDR(args);
  maxit_nt = asInteger(CAR(args)); args = CDR(args);
  tini = asReal(CAR(args)); args = CDR(args);
  progress = asInteger(CAR(args)); args = CDR(args);
  level = asReal(CAR(args)); args = CDR(args);
  message = asInteger(CAR(args)); args = CDR(args);
  vb = asInteger(CAR(args)); args = CDR(args);
  
  if(vb==-1 || vb==2) i0=1;
  if(vb== 1 || vb==2) i1=1;
  k =M[1]-M[0]; 
  lp=M[0]*(k+1)+(k+1)*(k+2)/2;
  phat = Calloc(lp, double);
  pval = Calloc(k+1, double);
  chpts = Calloc(k+1, int);
  tmp = 1.0+(double)(k*(k-1));
  wt = Calloc(M[1]+1, double);
  wt_ahat = Calloc(lp, double);
  p= Calloc(M[1]+1, double);
  llik = Calloc(1, double);
  lk = Calloc(k+1, double);
  lr = Calloc(k+1, double);
  z = Calloc(n0+n1, double);
  cp = Calloc(1, int);
  res = Calloc(1, double);

  for(i=0;i<n0;i++) z[i] = x0[i];
  for(i=0;i<n1;i++) z[i+n0] = x1[i];
  lr[0]=.0;
  m = M[0];  
  //for(j=0;j<=m;j++) p[j]=1.0/(double)(m+1.0);
  if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
  else{
    for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
    if(vb==-1 || vb==2) p[0]=0.0;
    if(vb== 1 || vb==2) p[m]=0.0;
  }
  pval[0]=1.0;
  chpts[0]=0;
  itmp=m+1;
  if(progress==1 && k>0)  ProgressBar(0.0,"");
  if(m>0){
    for(i=0; i<=m; i++){
       //p[i] = pini[j];
       p[i] = 1.0/(double)(m+1.0);
    }
    em_ptilde_dr(llik, alpha, p, x0, x1, rx1, m, d, n0, n1, wt, eps_em,  
      maxit_em, eps_nt, maxit_nt, (void*)&is);
    lk[0] = llik[0];
    for(i=0;i<itmp;i++){
      phat[i]=p[i];
      wt_ahat[i]=wt[i];}
  }
  else{
    phat[0]=1.0;
    lk[0]=loglik_unif(alpha, rx1, n1, d);
    wt_alpha(alpha, d, m, wt, (void*)&is);
    wt_ahat[0]=wt[0];
  }
  pct = level/pval[0];
  if(progress==1)  ProgressBar(fmax2(pct,1.0/tmp),"");
  i=1;
  while(i<k && pval[i-1]>level){
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
      for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
      if(vb==-1 || vb==2) p[0]=0.0;
      if(vb== 1 || vb==2) p[m]=0.0;
    }
    p[m+1] = (m+1)*p[m]/(double)(m+2);
    for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
    p[0] = (m+1.)*p[0]/(double)(m+2.);
    m = M[0]+i;
    //make sure initial p is in the interior of the simplex
    if(m>3) for(j=i0; j<=m-i1; j++) p[j] =(1.0-tini)*p[j]+ tini/(double)(m+1-abs(vb));
    em_ptilde_dr(llik, alpha, p, x0, x1, rx1, m, d, n0, n1, wt, eps_em, 
      maxit_em, eps_nt, maxit_nt, (void*)&is);
    lk[i] = llik[0];
    for(j=0;j<=m;j++){
      phat[j+itmp]=p[j];
      wt_ahat[j+itmp]=wt[j];}
    itmp += m+1;
    if(i>=3){
        cp[0]=i;
        chpt_exp(lk, lr, res, cp);
        pval[i]=res[0];
        chpts[i]=cp[0];
    }
    else{            
        pval[i]=1.0;
        chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
        cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
        cp0=cp1;
        pv0=pv1;
    }
    else pv0=pval[i];
    R_CheckUserInterrupt();
    pct = fmax2(pct, level/pval[i]);
    pct = fmax2((1.0+i*(i+1))/tmp, pct);
    if(progress==1){
        ProgressBar(fmin2(1.0,pct)," ");}
    i++;
  }
  if(m==M[1] && pval[i-1]>level && message){
     //convergence[0]+=1; 
     warning("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value,  %f, of the change-point is returned.\n", pv0);
  }
  M[1]=m;
  itmp=cp0*(M[0]*2+(cp0+1))/2;
  optim=cp0+M[0];
  m = optim;
  k = M[1]-M[0];
  if(progress==1){
    ProgressBar(1.0," ");
    Rprintf("\n");}
  //if(progress==1) Rprintf("\n 'maple.dr' done. \n\n\n");
  PROTECT(ans = allocVector(VECSXP, 9));
  PROTECT(ansnames = allocVector(STRSXP, 9));
  SET_STRING_ELT(ansnames, 0, mkChar("lk"));
  SET_STRING_ELT(ansnames, 1, mkChar("lr"));
  SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, m+1));
  SET_STRING_ELT(ansnames, 2, mkChar("p"));
  SET_STRING_ELT(ansnames, 3, mkChar("m"));
  SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, 1));
  INTEGER(VECTOR_ELT(ans, 3))[0] = m;
  SET_VECTOR_ELT(ans, 4, allocVector(REALSXP, m+1));
  SET_STRING_ELT(ansnames, 4, mkChar("wt"));
  for(i=0;i<=m;i++){
    REAL(VECTOR_ELT(ans, 2))[i] = phat[i+itmp];
    REAL(VECTOR_ELT(ans, 4))[i] = wt_ahat[i+itmp];
  }
  SET_VECTOR_ELT(ans, 5, allocVector(REALSXP, d+1));
  SET_STRING_ELT(ansnames, 5, mkChar("alpha"));
  for(i=0;i<=d;i++){
    REAL(VECTOR_ELT(ans, 5))[i] = alpha[i];
  }
  SET_VECTOR_ELT(ans, 6, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 7, allocVector(INTSXP, k+1));
  SET_STRING_ELT(ansnames, 6, mkChar("pval"));
  SET_STRING_ELT(ansnames, 7, mkChar("chpts"));
  for(i=0;i<=k;i++){
    REAL(VECTOR_ELT(ans, 0))[i] = lk[i];
    REAL(VECTOR_ELT(ans, 1))[i] = lr[i];
    REAL(VECTOR_ELT(ans, 6))[i] = pval[i];
    INTEGER(VECTOR_ELT(ans, 7))[i] = chpts[i]+M[0];
  }
  SET_STRING_ELT(ansnames, 8, mkChar("M"));
  SET_VECTOR_ELT(ans, 8, allocVector(INTSXP, 2));
  for(i=0;i<=1;i++){
    INTEGER(VECTOR_ELT(ans, 8))[i] = M[i];
  }
  
  Free(lk); Free(lr); Free(wt); Free(wt_ahat); Free(p); Free(phat); 
  Free(chpts); Free(pval); Free(res); Free(llik); Free(cp); Free(z);
  setAttrib(ans, R_NamesSymbol, ansnames);
  UNPROTECT(2);
  return ans;
}



// MABLE DR MODEL for Grouped Data
/*//////////////////////////////////////////////////////////*/
/*      EM Algorithm for (alpha, p) for grouped  data       */
/*//////////////////////////////////////////////////////////*/
// mablem_dr_group(llik, alpha, p, t, n0, n1, nx, ny, N, m, d, wt,
//       eps_em, eps_nt, maxit_em, maxit_nt, progress, ex);
// T_zero: (m+1)-vector of T0[j]=sum{n0[i]*p[j]*Pm[i,j]/sum(p*Pm[i,]): i=1,...,N}, j=0,...,m
void mablem_dr_group(double *llik, double *alpha, double *p, double *t, int *n0, int *n1,
       int nx, int ny, int N, int m, int d, double *wt, double *se, double eps_em,
       double eps_nt, int maxit_em, int maxit_nt, int progress, void *ex){
  int i, j, k, l, it_nt, it_em, d1=d+1, m1=m+1, Nm=N*m1,Nmd=Nm*d1;
  double *Pm, *Fm, *Fm_s, *H, *Jac, llik_new;
  double del_em, del_nt, *tmp, *h, *Sig, *Er, *dFm_s, *ddFm_s;//, eps0=10.0, diff0=10.0;
  double *Tk, *Pm_s, *dPm_s, *ddPm_s, *Pi_theta, *T_zero;
  T_zero = Calloc(m1,double); Tk = Calloc(m1,double);
  Pi_theta = Calloc(Nm,double); tmp = Calloc(d1,double);
  Pm = Calloc(Nm,double); Pm_s = Calloc(Nm,double);
  dPm_s = Calloc(Nmd,double); ddPm_s = Calloc(Nmd*d1,double);
  Fm = Calloc(N,double); Fm_s = Calloc(N,double);
  H = Calloc(d1,double); Jac = Calloc(d1*d1,double);

  h = Calloc(d1*d1,double); Sig = Calloc(d1*d1,double);  
  dFm_s = Calloc(N*d1,double); ddFm_s = Calloc(N*d1*d1,double);
  Er = Calloc(d1,double);
  
  cpBeta(t, m, N, Pm); //cdf_Beta(t, m, N, Pm); 
  Pm_alpha(alpha, t, N, d, m, Pm_s, ex);
  del_em=10.0; it_em=0;
  llik[0]=loglik_bern_group(p, N, n0, n1, Pm, Pm_s, m, d);
  while(del_em>eps_em && it_em<maxit_em){
    R_CheckUserInterrupt();
    for(i=0;i<N;i++){
       Fm[i] = 0.0; Fm_s[i] = 0.0;
       for(j=0;j<=m;j++){
          Fm[i] += Pm[i+N*j]*p[j];
          Fm_s[i] += Pm_s[i+N*j]*p[j];
       }
       for(j=0;j<=m;j++) Pi_theta[i+N*j] = n1[i]*Pm_s[i+N*j]*p[j]/Fm_s[i];
    }
    for(j=0;j<=m;j++){
       T_zero[j]=0.0;
       for(i=0;i<N;i++) T_zero[j]+=n0[i]*Pm[i+N*j]*p[j]/Fm[i];
    }
    score_alpha_group(alpha, t, N, nx, ny, d, m, Pi_theta, T_zero, Tk, H, Jac, wt, ex);
    del_nt=0.0;
    for(i=0;i<=d;i++) del_nt += fabs(H[i]);
    it_nt=1;
    while(del_nt>eps_nt && it_nt<maxit_nt){
      minverse(Jac, d1);
      for(i=0;i<=d;i++){
         tmp[i] = 0.0;
         for(j=0;j<=d;j++) tmp[i] += Jac[i+d1*j]*H[j];
      }
      del_nt = 0.0;
      for(i=0;i<=d;i++){
         alpha[i] -= tmp[i];
         del_nt += fabs(tmp[i]);
         //Rprintf(" alpha[%d] = %g\n",i, alpha[i]);
      }
      dPm_alpha(alpha, t, N, d, m, Pm_s, dPm_s, ddPm_s, ex);
      for(i=0;i<N;i++){
         Fm_s[i] = 0.0;
         for(j=0;j<=m;j++)  Fm_s[i] += Pm_s[i+N*j]*p[j];
         for(j=0;j<=m;j++) Pi_theta[i+N*j] = n1[i]*Pm_s[i+N*j]*p[j]/Fm_s[i];
      }
      score_alpha_group(alpha, t, N, nx, ny, d, m, Pi_theta, T_zero, Tk, H, Jac, wt, ex);
      for(i=0;i<=d;i++) del_nt += fabs(H[i]);
      //Rprintf("m = %d,  Newton: it_nt=%d, del_nt=%g\n", m, it_nt, del_nt);
      it_nt++;
    }
    for(j=0;j<=m;j++) p[j] = Tk[j]/(nx+ny*wt[j]);
    //for(j=0;j<=m;j++) p[j] = (1.0-veps)*p[j]+veps*Tk[j]/(nx+ny*wt[j]);
    dPm_alpha(alpha, t, N, d, m, Pm_s, dPm_s, ddPm_s, ex);
    llik_new = loglik_bern_group(p, N, n0, n1, Pm, Pm_s, m, d);
    del_em = fabs(llik_new-llik[0]);
    llik[0] = llik_new;
    it_em++;
    //if(progress==1) {
    //    if(it_em==1) {eps0 = del_em; diff0=eps0-eps_em;}
    //    pct = R_pow_di((eps0-del_em)/diff0, 1000);
    //    ProgressBar(pct," ");}
    if(progress==1) clockProgress((int)(it_em/50)," EM iteration is runing...");
    //Rprintf("m = %d,  Newton: it=%d, del=%g\n", m, it_nt, del_nt);
    //Rprintf("\n m = %d, EM: it_em=%d, del_em=%g\n", m, it_em, del_em);
  }
  if(progress==1) Rprintf("\n");
  //Rprintf("m = %d, EM: nit=%d, del=%g\n", m, it_em, del_em);  
  
  for(k=0;k<N;k++){
    Fm_s[k] = 0.0;
    for(i=0;i<=m;i++) Fm_s[k] += Pm_s[k+N*i]*p[i];
    for(i=0;i<=d;i++){
      dFm_s[k+N*i]=0.0;
      for(j=0;j<=m;j++) dFm_s[k+N*i]+=p[j]*dPm_s[k+N*j+Nm*i];
      for(j=0;j<=d;j++){
        ddFm_s[k+N*i+N*d1*j]=0.0;
        for(l=0;l<=m;l++) ddFm_s[k+N*i+N*d1*j]+=p[l]*ddPm_s[k+N*l+Nm*i+Nmd*j];
      }
    }
  }
  for(i=0;i<=d;i++){
    Er[i]=0.0;
    for(k=0;k<N;k++) Er[i]+=dFm_s[k+N*i];
    for(j=0;j<=i;j++){
      h[i+j*d1]=0.0;
      Sig[i+j*d1]=0.0;
      for(k=0;k<N;k++){ 
        h[i+j*d1]+=(1.0-Fm_s[k])*ddFm_s[k+N*i+N*d1*j]+dFm_s[k+N*i]*dFm_s[k+N*j];
        Sig[i+j*d1]+=dFm_s[k+N*i]*dFm_s[k+N*j]/Fm_s[k];
      }
      Sig[i+j*d1]-=Er[i]*Er[j];
      h[j+i*d1]=h[i+j*d1];
      Sig[j+i*d1]=Sig[i+j*d1];
    }
  }
  minverse(h, d1);
  for(k=0;k<=d;k++){
    se[k]=0.0;
    for(i=0;i<=d;i++){
      for(j=0;j<=d;j++)
        se[k]+=h[k+i*d1]*Sig[i+j*d1]*h[j+k*d1];
    }
    se[k]=sqrt(se[k]/(double) ny);
  }
  
  
  Free(H); Free(Jac); Free(Tk); 
  Free(Er); Free(h); Free(Sig); Free(dFm_s); Free(ddFm_s);
  Free(tmp); Free(T_zero); Free(Pi_theta); Free(Pm); Free(Pm_s);
  Free(dPm_s); Free(ddPm_s); Free(Fm); Free(Fm_s);
}
/*/////////////////////////////////////////////////////////////*/
/*     Model degree m selection by change-point method         */
/*  using approximate loglikelihood ell(alpha-tilde,p-tilde)   */
/*  alpha-tilde:  estimated logistic regression coefficient    */
/*      p-tilde: obtained by iteration p=A(alpha-tilde,p)      */
/*/////////////////////////////////////////////////////////////*/

/*//////////////////////////////////////////////////*/
/*      EM Method for p-tilde for grouped data      */
/*//////////////////////////////////////////////////*/
//em_ptilde_dr_group(llik, alpha, p, t, n0, n1, nx, ny, N, m, d, wt, eps_em, maxit_em, eps_nt, maxit_nt, ex);
// alpha: MELE of alpha: alpha-tilde
//     p: initial value for p
void em_ptilde_dr_group(double *llik, double *alpha, double *p, double *t, 
       int *n0, int *n1, int nx, int ny, int N, int m, int d, double *wt,  
       double eps_em, int maxit_em, double eps_nt, int maxit_nt, void *ex){
  int i, j, it_em, it_nt, n, m1=m+1, Nm=N*m1;
  double *Tk, *Pm, *Pm_s, *Fm, *Fm_s, llik_new;
  double del_em, del_nt, lam, tmp;
  
  n=nx+ny;
  Tk = Calloc(m1,double);
  Pm = Calloc(Nm,double); Pm_s = Calloc(Nm,double);
  Fm = Calloc(N,double); Fm_s = Calloc(N,double);

  wt_alpha(alpha, d, m, wt, ex);
  cpBeta(t, m, N, Pm); //cdf_Beta(t, m, N, Pm);  
  Pm_alpha(alpha, t, N, d, m, Pm_s, ex);
  del_em=10.0;
  it_em=1;
  llik[0]=loglik_bern_group(p, N, n0, n1, Pm, Pm_s, m, d);
  while(del_em>eps_em && it_em<maxit_em){
    R_CheckUserInterrupt();
    for(i=0;i<N;i++){
       Fm[i] = 0.0; Fm_s[i] = 0.0;
       for(j=0;j<=m;j++){
          Fm[i] += Pm[i+N*j]*p[j];
          Fm_s[i] += Pm_s[i+N*j]*p[j];
       }
    }
    for(j=0;j<=m;j++){
       Tk[j]=0.0;
       for(i=0;i<N;i++) Tk[j]+=n0[i]*Pm[i+N*j]*p[j]/Fm[i];
       for(i=0;i<N;i++) Tk[j]+=n1[i]*Pm_s[i+N*j]*p[j]/Fm_s[i];
    }
    // Newton method for finding lambda starting with lambda0=ny
    lam = ny; tmp = 0.0;
    for(j=0;j<=m;j++) tmp += Tk[j]*(wt[j]-1.0)/(n+lam*(wt[j]-1.0));
    del_nt = fabs(tmp); it_nt = 0;
    while(del_nt>eps_nt && it_nt<maxit_nt){
        del_nt = 0.0;
        for(j=0;j<=m;j++) del_nt+= Tk[j]*(wt[j]-1.0)*(wt[j]-1.0)/((n+lam*(wt[j]-1.0))*(n+lam*(wt[j]-1.0)));
        del_nt = tmp/del_nt;
        lam += del_nt;
        del_nt = fabs(del_nt);
        tmp = 0.0;
        for(j=0;j<=m;j++) tmp += Tk[j]*(wt[j]-1.0)/(n+lam*(wt[j]-1.0));
        del_nt += fabs(tmp);
        it_nt++;
    }
    //Rprintf("m = %d, Lambda=%g: it=%d, del=%g\n", m, lam, it_nt, del_nt);
/////////
    for(j=0;j<=m;j++) p[j] = Tk[j]/(n+lam*(wt[j]-1.0));
    llik_new = loglik_bern_group(p, N, n0, n1, Pm, Pm_s, m, d);
    del_em = fabs(llik_new-llik[0]);
    llik[0] = llik_new;
    it_em++;
    //Rprintf("m = %d, EM: it_em=%d, del_em=%g\n", m, it_em, del_em);
  }
  //Rprintf("m = %d, EM: nit=%d, del=%g\n", m, it_em, del_em);
  Free(Tk); Free(Pm); Free(Pm_s); Free(Fm); Free(Fm_s);
}
/*////////////////////////////////////////////////////*/
/*  Model degree m selection by change-point method   */
/*  with a fixed alpha-tilde based on grouped data    */
/*////////////////////////////////////////////////////*/
// alpha: MELE of alpha, alpha-tilde
//     p: initial value for m=mu
SEXP maple_dr_group(SEXP args){
  //if(progress==1) 
  //  Rprintf("\n Program 'maple.dr.group' is runing. This may take several minutes.\n\n\n");
  int_struct is;
  SEXP ans, ansnames;
  int i, j, m, k, d, nx, ny, *n0, *n1, N, optim, *M, maxit_nt, maxit_em, progress;
  int cp0=0, cp1=1, *cp, *chpts, itmp, lp, message, vb, i0=0, i1=0;
  double  *phat, *wt, *llik, *res, pv0=1.0, pv1=1.0, level, *pval; 
  double *alpha, *p, *lr, *lk, *t, eps_em, eps_nt, tmp, pct, tini=0.0001;
  
  args = CDR(args);
  is.f = CAR(args); args = CDR(args); 
  is.env = CAR(args); args = CDR(args);
  alpha = REAL(CAR(args)); args = CDR(args); 
  t = REAL(CAR(args)); args = CDR(args); 
  n0 = INTEGER(CAR(args)); args = CDR(args);
  n1 = INTEGER(CAR(args)); args = CDR(args); 
  nx = asInteger(CAR(args)); args = CDR(args);
  ny = asInteger(CAR(args)); args = CDR(args);
  N = asInteger(CAR(args)); args = CDR(args);
  M = INTEGER(CAR(args)); args = CDR(args);
  d = asInteger(CAR(args)); args = CDR(args);
  eps_em = asReal(CAR(args)); args = CDR(args); 
  eps_nt = asReal(CAR(args)); args = CDR(args); 
  maxit_em = asInteger(CAR(args)); args = CDR(args); 
  maxit_nt = asInteger(CAR(args)); args = CDR(args); 
  progress = asInteger(CAR(args)); args = CDR(args);
  level = asReal(CAR(args)); args = CDR(args);
  message = asInteger(CAR(args)); args = CDR(args);
  vb = asInteger(CAR(args)); args = CDR(args);
      
  if(vb==-1 || vb==2) i0=1;
  if(vb== 1 || vb==2) i1=1;
  k = M[1]-M[0]; 
  lp=M[0]*(k+1)+(k+1)*(k+2)/2;
  phat = Calloc(lp, double);
  pval = Calloc(k+1, double);
  chpts = Calloc(k+1, int);
  tmp = 1.0+(double)(k*(k+1));
  //Rprintf("  dim:  %d\n", (k+1)*(d+1));
  lk = Calloc(k+1, double);
  lr = Calloc(k+1, double);
  p = Calloc(M[1]+1, double);
  wt = Calloc(M[1]+1, double);
  llik = Calloc(1, double);
  cp = Calloc(1, int);
  res = Calloc(1, double);

  m = M[0];
  //for(i=0; i<=m; i++) p[i] = 1.0/(double)(m+1.0);
  if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
  if(m>=3){
    for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
    if(vb==-1 || vb==2) p[0]=0.0;
    if(vb== 1 || vb==2) p[m]=0.0;
  }
  lr[0] = .0;
  pval[0]=1.0;
  chpts[0]=0;
  itmp=m+1;
  if(m>0){
    //for(i=0; i<=d; i++) Rprintf(" alpha[%d] = %g\n",i, alpha[i]);
    em_ptilde_dr_group(llik, alpha, p, t, n0, n1, nx, ny, N, m, d, wt, 
        eps_em, maxit_em, eps_nt, maxit_nt, (void*)&is);
    lk[0] = llik[0];
    for(i=0;i<itmp;i++) phat[i]=p[i]; 
  }
  else{
      phat[0]=1.0;
      lk[0]=loglik_unif_group(t,N,n0,n1,alpha,m,d,(void*)&is);    
  }
  pct = level/pval[0];
  if(progress==1)  ProgressBar(fmax2(pct,1.0/tmp),""); 
  i=1;
  while(i<k && pval[i-1]>level){
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
      for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
      if(vb==-1 || vb==2) p[0]=0.0;
      if(vb== 1 || vb==2) p[m]=0.0;
    }
    p[m+1] = (m+1)*p[m]/(double)(m+2);
    for(j=m; j>=1; j--) phat[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
    phat[0] = (m+1.)*p[0]/(double)(m+2.);
    m = M[0]+i;
    //make sure initial p is in the interior of the simplex
    if(m>3) for(j=i0; j<=m-i1; j++) p[j] =(1.0-tini)*p[j]+tini/(double)(m+1-abs(vb));
    //for(j=0; j<=m; j++)  p[j]=1.0/(double)(m+1.0);
    em_ptilde_dr_group(llik, alpha, p, t, n0, n1, nx, ny, N, m, d, wt, 
        eps_em, maxit_em, eps_nt, maxit_nt, (void*)&is);
    //for(j=0; j<=m; j++) Rprintf("p[%d] = %g,\t",j, p[j]);
    //Rprintf("\n");
    lk[i] = llik[0];
   //Rprintf("                     lk[%d]= %g\n", i, lk[i]);
    for(j=0;j<=m;j++) phat[j+itmp]=p[j];
    itmp += m+1;
    if(i>=3){
        cp[0]=i;
        chpt_exp(lk, lr, res, cp);
        pval[i]=res[0];
        chpts[i]=cp[0];
    }
    else{            
        pval[i]=1.0;
        chpts[i]=0;
    }
    if(chpts[i]>chpts[i-1]){
        cp1=chpts[i];
    }
    if(cp0<cp1) pv1=pval[i];
    else pv0=pval[i];
    if(pv1<pv0){
        cp0=cp1;
        pv0=pv1;
    }
    else pv0=pval[i];
    if(progress==1){
        pct = fmax2(pct, level/pval[i]);
        pct = fmax2((1.0+i*(i+1))/tmp, pct);
        ProgressBar(fmin2(1.0,pct)," ");}
    i++;
  }
  if(m==M[1] && pval[i-1]>level && message){
      //convergence[0]+=1; 
      warning("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value,  %f, of the change-point is returned.\n", pv0);
  }
  M[1]=m;
  itmp=cp0*(M[0]*2+(cp0+1))/2;
  optim=cp0+M[0];
  m = optim;
  k = M[1]-M[0];
  if(progress==1){
    ProgressBar(1.0," ");
    Rprintf("\n");}
  //if(progress==1) Rprintf("\n 'maple.dr.group' done. \n\n\n");
  
  PROTECT(ans = allocVector(VECSXP, 7));
  PROTECT(ansnames = allocVector(STRSXP, 7));
  SET_STRING_ELT(ansnames, 0, mkChar("lk"));
  SET_STRING_ELT(ansnames, 1, mkChar("lr"));
  SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, m+1));
  SET_STRING_ELT(ansnames, 2, mkChar("p"));
  for(i=0;i<=m;i++){
      REAL(VECTOR_ELT(ans, 2))[i] = phat[i];
  }
  SET_STRING_ELT(ansnames, 3, mkChar("m"));
  SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, 1));
  INTEGER(VECTOR_ELT(ans, 3))[0] = m;

  SET_VECTOR_ELT(ans, 4, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 5, allocVector(INTSXP, k+1));
  SET_STRING_ELT(ansnames, 4, mkChar("pval"));
  SET_STRING_ELT(ansnames, 5, mkChar("chpts"));
  for(i=0;i<=k;i++){
      REAL(VECTOR_ELT(ans, 0))[i] = lk[i];
      REAL(VECTOR_ELT(ans, 1))[i] = lr[i];
      REAL(VECTOR_ELT(ans, 4))[i] = pval[i];
      INTEGER(VECTOR_ELT(ans, 5))[i] = chpts[i]+M[0];
  }
  SET_STRING_ELT(ansnames, 6, mkChar("M"));
  SET_VECTOR_ELT(ans, 6, allocVector(INTSXP, 2));
  for(i=0;i<=1;i++){
      INTEGER(VECTOR_ELT(ans, 6))[i] = M[i];
  }

  Free(llik); Free(wt); Free(phat); Free(p); Free(lk); Free(lr);
  Free(cp); Free(res); Free(cp); Free(res);
  setAttrib(ans, R_NamesSymbol, ansnames);
  UNPROTECT(2);
  return ans;
}

/*////////////////////////////////////////////////////*/
/*  Model degree m selection by change-point method   */
/*    with (alpha-hat, phat) based on grouped data    */
/*////////////////////////////////////////////////////*/
SEXP C_mable_dr_group(SEXP args){
//////
  //if(progress==1) 
  //  Rprintf("\n Program 'mable.dr.group' is runing. This may take several minutes.\n\n\n");
  int_struct is;
  SEXP ans, ansnames;
  int i, j, m, k, d, nx, ny, *n0, *n1, N, optim, *M, maxit_nt, maxit_em, progress;
  int *chpts, itmp, cp0=0, cp1=1, *cp, lp, message, vb, i0=0, i1=0;
  double *alpha_hat, *alpha_ini, *phat, *wt, *se, *se_hat, *llik, pv0=1.0, pv1=1.0, *pval, level; 
  double *alpha, *p, *lr, *lk, *t, eps_em, eps_nt, tmp, pct, *res, tini=0.0001;
  
  args = CDR(args);
  is.f = CAR(args); args = CDR(args); 
  is.env = CAR(args); args = CDR(args);
  alpha_ini = REAL(CAR(args)); args = CDR(args); 
  //p = REAL(CAR(args)); args = CDR(args);
  t = REAL(CAR(args)); args = CDR(args); 
  n0 = INTEGER(CAR(args)); args = CDR(args);
  n1 = INTEGER(CAR(args)); args = CDR(args); 
  nx = asInteger(CAR(args)); args = CDR(args);
  ny = asInteger(CAR(args)); args = CDR(args);
  N = asInteger(CAR(args)); args = CDR(args);
  M = INTEGER(CAR(args)); args = CDR(args);
  d = asInteger(CAR(args)); args = CDR(args);
  //veps = asReal(CAR(args)); args = CDR(args); 
  eps_em = asReal(CAR(args)); args = CDR(args); 
  eps_nt = asReal(CAR(args)); args = CDR(args); 
  maxit_em = asInteger(CAR(args)); args = CDR(args); 
  maxit_nt = asInteger(CAR(args)); args = CDR(args); 
  progress = asInteger(CAR(args)); args = CDR(args);
  level = asReal(CAR(args)); args = CDR(args); 
  message = asInteger(CAR(args)); args = CDR(args);
  vb = asInteger(CAR(args)); args = CDR(args);
     
  if(vb==-1 || vb==2) i0=1;
  if(vb== 1 || vb==2) i1=1;
  k = M[1]-M[0];// tini=mu*.000001;
  lp=M[0]*(k+1)+(k+1)*(k+2)/2;
  phat = Calloc(lp, double);
  pval = Calloc(k+1, double);
  chpts = Calloc(k+1, int);
  tmp = 1.0+(double)(k*(k-1));
  //Rprintf("  dim:  %d\n", (k+1)*(d+1));
  lk = Calloc(k+1, double);
  lr = Calloc(k+1, double);
  alpha_hat = Calloc((k+1)*(d+1), double);
  alpha = Calloc(d+1, double);
  se = Calloc(d+1, double);
  se_hat = Calloc((k+1)*(d+1), double);
  p= Calloc(M[1]+1, double);
  wt = Calloc(M[1]+1, double);
  llik = Calloc(1, double);
  cp = Calloc(1, int);
  res = Calloc(1, double);
  for(i=0; i<=d; i++) alpha[i] = alpha_ini[i]; 
  m = M[0];
  //for(j=0;j<=m;j++) p[j]=1.0/(double)(m+1.0);
  if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
  if(m>=3){
    for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
    if(vb==-1 || vb==2) p[0]=0.0;
    if(vb== 1 || vb==2) p[m]=0.0;
  }
  lr[0] = .0;
  pval[0]=1.0;
  chpts[0]=0;
  itmp=m+1;
  pct = level/pval[0];
  if(progress==1)  ProgressBar(fmax2(pct,1.0/tmp),"");
  if(m>0){
    //for(i=0; i<=d; i++) Rprintf(" alpha[%d] = %g\n",i, alpha[i]);
    //for(i=0; i<=d; i++)  alpha[i] = alpha_ini[i];
    mablem_dr_group(llik, alpha, p, t, n0, n1, nx, ny, N, m, d, wt, se,
        eps_em, eps_nt, maxit_em, maxit_nt, progress*(k==0), (void*)&is);
    lk[0] = llik[0];
    for(i=0;i<itmp;i++) phat[i]=p[i]; 
  }
  else{
      phat[0]=1.0;
      lk[0]=loglik_unif_group(t,N,n0,n1,alpha,m,d,(void*)&is);    
  }
  //if(progress==1)  ProgressBar(1.0/tmp,"");
  for(j=0; j<=d; j++){
    alpha_hat[j] = alpha[j];
    se_hat[j] = se[j];
  }
  //Rprintf("                     lk[%d]= %g\n", 0, lk[0]);
  i=1;
  while(i<k && pval[i-1]>level){
    if(m<=2) for(j=0; j<=m; j++) p[j] = 1.0/(double)(m+1);
    if(m==3){
       for(j=i0; j<=m-i1; j++) p[j] = 1.0/(double)(m+1-abs(vb));
       if(vb==-1 || vb==2) p[0]=0.0;
       if(vb== 1 || vb==2) p[m]=0.0;
    }
    p[m+1] = (m+1)*p[m]/(double)(m+2);
    for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m-j+1))/(double)(m+2);
    phat[0] = (m+1.)*phat[0]/(double)(m+2.);
    m = M[0]+i;
    //make sure initial p is in the interior of the simplex
    if(m>3) for(j=i0; j<=m-i1; j++) p[j] =(1.0-tini)*p[j]+tini/(double)(m+1-abs(vb));
    for(j=0; j<=d; j++) alpha[j] = alpha_ini[j];
    //for(j=0; j<=m; j++) p[j]=1.0/(double)(m+1.0);      
    mablem_dr_group(llik, alpha, p, t, n0, n1, nx, ny, N, m, d, wt, se, eps_em, eps_nt,
          maxit_em, maxit_nt, 0, (void*)&is);
      //for(j=0; j<=m; j++) Rprintf("p[%d] = %g,\t",j, phat[j]);
      //Rprintf("\n");
      lk[i] = llik[0];
      //Rprintf("                     lk[%d]= %g\n", i, lk[i]);
      for(j=0;j<=m;j++) phat[j+itmp]=p[j];
      for(j=0; j<=d; j++){
        alpha_hat[j+i*(d+1)] = alpha[j]; 
        se_hat[j+i*(d+1)] = se[j];
      } 
      itmp += m+1;
      if(i>=3){
          cp[0]=i;
          chpt_exp(lk, lr, res, cp);
          pval[i]=res[0];
          chpts[i]=cp[0];
      }
      else{            
          pval[i]=1.0;
          chpts[i]=0;
      }
      if(chpts[i]>chpts[i-1]){
          cp1=chpts[i];
      }
      if(cp0<cp1) pv1=pval[i];
      else pv0=pval[i];
      if(pv1<pv0){
          cp0=cp1;
          pv0=pv1;
      }
      else pv0=pval[i];
      R_CheckUserInterrupt();
      pct = fmax2(pct, level/pval[i]);
      pct = fmax2(pct, (1.0+i*(i+1))/tmp);
      if(progress==1) {
          ProgressBar(fmin2(1.0,pct)," ");}
      i++;
  }
  if(m==M[1] && pval[i-1]>level && message){
      //convergence[0]+=1; 
      warning("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value,  %f, of the change-point is returned.\n", pv0);
  }
  M[1]=m;
  itmp=cp0*(M[0]*2+(cp0+1))/2;
  optim=cp0+M[0];
  m = optim;
  k = M[1]-M[0];
  if(progress==1){
    ProgressBar(1.0," ");
    Rprintf("\n");}
  //if(progress==1) Rprintf("\n 'mable.dr.group' done. \n\n\n");
  
  PROTECT(ans = allocVector(VECSXP, 9));
  PROTECT(ansnames = allocVector(STRSXP, 9));
  SET_STRING_ELT(ansnames, 0, mkChar("lk"));
  SET_STRING_ELT(ansnames, 1, mkChar("lr"));
  SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 2, allocVector(REALSXP, d+1));
  SET_STRING_ELT(ansnames, 2, mkChar("alpha"));
  SET_VECTOR_ELT(ans, 3, allocVector(REALSXP, m+1));
  SET_STRING_ELT(ansnames, 3, mkChar("p"));
  for(i=0;i<=m;i++){
      REAL(VECTOR_ELT(ans, 3))[i] = phat[i+itmp];
  }
  SET_STRING_ELT(ansnames, 4, mkChar("m"));
  SET_VECTOR_ELT(ans, 4, allocVector(INTSXP, 1));
  INTEGER(VECTOR_ELT(ans, 4))[0] = m;

  SET_VECTOR_ELT(ans, 5, allocVector(REALSXP, k+1));
  SET_VECTOR_ELT(ans, 6, allocVector(INTSXP, k+1));
  SET_STRING_ELT(ansnames, 5, mkChar("pval"));
  SET_STRING_ELT(ansnames, 6, mkChar("chpts"));
  for(i=0;i<=k;i++){
      REAL(VECTOR_ELT(ans, 0))[i] = lk[i];
      REAL(VECTOR_ELT(ans, 1))[i] = lr[i];
      REAL(VECTOR_ELT(ans, 5))[i] = pval[i];
      INTEGER(VECTOR_ELT(ans, 6))[i] = chpts[i]+M[0];
  }
  SET_STRING_ELT(ansnames, 7, mkChar("M"));
  SET_VECTOR_ELT(ans, 7, allocVector(INTSXP, 2));
  for(i=0;i<=1;i++){
      INTEGER(VECTOR_ELT(ans, 7))[i] = M[i];
  }
  SET_VECTOR_ELT(ans, 8, allocVector(REALSXP, d+1));
  SET_STRING_ELT(ansnames, 8, mkChar("se"));
  for(i=0;i<=d;i++){
      REAL(VECTOR_ELT(ans, 2))[i] = alpha_hat[i+cp0*(d+1)];
      REAL(VECTOR_ELT(ans, 8))[i] = se_hat[i+cp0*(d+1)];
  }

  Free(llik); Free(wt); Free(p); Free(alpha_hat); Free(alpha); Free(phat); 
  Free(lk); Free(lr); Free(chpts); Free(pval); Free(cp); Free(res); 
  Free(se); Free(se_hat);
  setAttrib(ans, R_NamesSymbol, ansnames);
  UNPROTECT(2);
  return ans;
}


