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
//#include <windows.h>
//#include <Defn.h>
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
#include <R_ext/Print.h>	/* for Rprintf */
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
/* Given a matrix A, this routine replaces it by the LU decomposition of A rowwise
permutation of itself. A and n are input. A is output, arranged as in A=(L-I)+U;
indx is an output vector that records the row permutation effected by the partial
pivoting; d is output as ±1 depending on whether the number of row interchanges 
was even or odd, respectively. This routine is used in combination with lubksb 
to solve linear equations or invert a matrix. */
void ludcmp(double *A, int n, int *indx, double *d){
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv; /*vv stores the implicit scaling of each row.*/
  vv=R_Calloc(n, double);
  imax=0;
  *d=1.0; /*No row interchanges yet.*/
  for (i=0;i<n;i++) { /*Loop over rows to get the implicit scaling information.*/
    big=0.0;
    for (j=0;j<n;j++)
        if ((temp=fabs(A[i+j*n])) > big) big=temp;
    if (big == 0.0) error("\nSingular matrix in routine ludcmp\n");
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
  R_Free(vv);
 }
/*//////////////////////////////////////////////////*/
/*          Forward and backsubstitution            */
/*//////////////////////////////////////////////////*/
/*Solves the set of n linear equations A·X = B. Here A is input, 
not as the matrix A but rather as its LU decomposition, determined 
by the routine ludcmp. indx is input as the permutation vector returned 
by ludcmp. b[] is input as the right-hand side vector B, and returns with 
the solution vector X. A, n, and indx are not modified by this routine
and can be left in place for successive calls with different right-hand 
sides b. This routine takes into account the possibility that b will 
begin with many zero elements, so it is efficient for use 
in matrix inversion. */
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

/*//////////////////////////////////////////////*/
/*  Newton Iteration x1=x0-J(x0)^{-1}F(x0)      */
/*   (1) Solve J(x0)x=F(x0), (2) x1=x0-x        */
/*//////////////////////////////////////////////*/
// input: A=J(x0), c=x0, b=F(x0)
//output: c=x1, del=sum(abs(x1-x0))
void newton_iter(double *A, int N, double *b, double *c, double *del){
  int i, *indx;
  double d=0.0;
  indx = R_Calloc(N, int); 
  ludcmp(A,N,indx,&d);
  /* Decompose the matrix.*/
  lubksb(A,N,indx,b);
  del[0]=0.0;
  for(i=0;i<N;i++) {
    del[0]+=fabs(b[i]);
    c[i]-=b[i];}
  R_Free(indx);
}


/*//////////////////////////////////////////////////*/
/*               Matrix Inverse                     */
/*//////////////////////////////////////////////////*/

void minverse(double *A, int N){
  double *tmp, d, *col;
  int i,j,*indx;
  d=0.0;
  indx = R_Calloc(N, int); 
  col = R_Calloc(N, double);
  tmp = R_Calloc(N*N, double);
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
  R_Free(tmp); R_Free(col); R_Free(indx);
}

/*//////////////////////////////////////////////////*/
/*         Check if Matrix is Singular              */
/*//////////////////////////////////////////////////*/

int matrix_singular(double *A, int n){
  int i,j, singular=0;
  double big, temp;
  for (i=0;i<n;i++) { 
    big=0.0;
    for (j=0;j<n;j++)
        if ((temp=fabs(A[i+j*n])) > big) big=temp;
    if (big == 0.0) singular=1;
  }
  return singular;
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
/* Bernstein base polynomials */
/* dBeta: Returns n x (m+1) matrix, n=length(u) */
/* dbeta(u[i], j+1, m+1-j)= the j-th column, j=0, 1, ..., m */
void dBeta(double *u, int m, int n, double *Bta) {
  int i, j;
  for(i=0; i<n; i++) Bta[i]=(m+1)*R_pow_di(1.0-u[i], m);
  for(i=0; i<n; i++) {
    if(u[i]<1){
      j=0;
      while(j<m){
        Bta[i+n*(j+1)]=(m-j)*(u[i]/(1.0-u[i]))*Bta[i+n*j]/(double)(j+1.0);
//          Rprintf("  Bta = %g\n",  Bta[i+n*(j+1)]);
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
//          Rprintf("  Bta = %g\n",  pBta[i+n*j]);
    }
  }
}
/*////////////////////////////////////////////////////////////*/
/*          Group Probability of beta(j+1, m-j+1)             */
/*    cpBeta: Returns N x (m+1) matrix, N=length(u)-1         */
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
/* Bdata: Returns n x (m+2) matrix B[0:(n-1), 0:(m+1)], n=length(y)  */
/* The first n0 rows are                                             */
/*   B[i,j]=dbeta(y[i], j+1, m+1-j), j=0:m; B[i,m+1]=0, i=0:(n0-1)   */
/* The other n1=n-n0 rows are                                        */
/*   B[i,j]=1-pbeta(y[i], j+1, m+1-j), j=0:m; B[i,m+1]=1, i=n0:(n-1) */
/*///////////////////////////////////////////////////////////////////*/
void Bdata(double *y, int m, int n0, int n1, double *Bta) {
  int i, j, n=n0+n1, mp1=m+1, nmp1=n*mp1;
  for(i=0; i<n0; i++){
    for(j=0;j<=m;j++) Bta[i+n*j]=dbeta(y[i], j+1, mp1-j, FALSE);
    Bta[i+nmp1]=0.0;}
  for(i=n0; i<n; i++){ 
    if(y[i]<=1.0){
      for(j=0;j<=m;j++) Bta[i+n*j]=1-pbeta(y[i], j+1, mp1-j, TRUE, FALSE);
      Bta[i+nmp1]=1.0;}
    else for(j=0;j<=mp1;j++) Bta[i+n*j]=0.0;  
  }
}
/*////////////////////////////////////////////////////*/
/*          Calculate fm(y,p) and Sm(y,p)             */
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
  for(i=0;i<n;i++){
    egx[i]=0.0;
    for(j=0;j<d;j++) egx[i]+= x[i+n*j]*gama[j];
    egx[i]=exp(egx[i]-gx0);
  }
}
/*////////////////////////////////////////////////////*/
/*       Calculate exp(g*x-tilde) and find            */
/*      x0=x(gama)=argmin{gamma'x_i: i=1:n}           */
/*////////////////////////////////////////////////////*/
void egx_x0(double *gama, int d, double *x, int n, double *egx, double *x0){
  int i,j,imin=0;
  double egx0=0.0, min_gx;
  //Rprintf("\n * x0=%f, gamma=%f\n",x0[0], gama[0]);
  for(j=0;j<d;j++) egx0 += x0[j]*gama[j];
  //egx0=exp(egx0);
  min_gx=egx0;
  for(i=0;i<n;i++){
    egx[i]=0.0;
    for(j=0;j<d;j++) egx[i]+= x[i+n*j]*gama[j];
    //egx[i]=exp(egx[i]);
    if(egx[i]<min_gx){
      min_gx=egx[i];
      imin=i;
    }
  }
  //Rprintf("\n ** x0=%f, gamma=%f\n",x0[0], gama[0]);
  if(egx0>min_gx){
    egx0=min_gx;
    for(j=0;j<d;j++) x0[j] = x[imin+n*j];
  }
  for(i=0;i<n;i++) egx[i] = exp(egx[i]-egx0);
}

/*///////////////////////////////////////*/
/*    Exponential Change-point Method    */
/*   for choosing optimal model degree   */
/*///////////////////////////////////////*/
// chpt[0] input k, ouput change-point
void chpt_exp(double *lk, double *lr, double *pv, int *chpt){
  int i, k = chpt[0]; 
  double mLR=0.0, lr0, lnk=log(k), llnk=log(lnk);
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
  res = R_Calloc(3, double);
  x = R_Calloc(k, double);
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
  R_Free(x); R_Free(res);
}
//choosing optimal degree based on lk using gamma chage-point model
void optim_gcp(int *M, double *lk, double *lr, int *m, 
      double *pval, int *chpts){
  int i, *cp, k=M[1]-M[0];
  double *res; 
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
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
  R_Free(cp); R_Free(res);
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
  ex = R_Calloc(4, int);
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
  R_Free(ex);
}

// Next combination
// a = (a1, a2, .., ar) a proper subset of (1,2,..,n) which 
//    is not equal to (n-r+1,..,n)  with a1<a2<...<ar
void next_combo(int *a, int r, int n){
    int i, j, d=0;
    if(r==0) error("\nEmpty Combination\n"); 
    if(n<r) error("\nIncorrect input n value\n");
    for(i=0;i<r;i++) d+=(a[i]==n-r+1+i);
    if(d==r) error("\nInput 'a' is the last combination.\n");
    i=r-1;
    while(a[i]==n-r+1+i) i--;
    a[i]++;
    j=i+1;
    while(j<r){
        a[j]=a[i]+j-i;
        j++;
    }
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
    pBta = R_Calloc((m+1)*n, double);
    fp = R_Calloc(n, double);
    pnu = R_Calloc(m+1, double);
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
    R_Free(pBta); R_Free(fp); R_Free(pnu);
}
/* end function em_beta_mix */
void mable_em(int *m, int *n, double *p, double *x, int *maxit,  double *eps,
    double *llik, int *convergence, double *delta){
    double *Bta;
    Bta = R_Calloc((*m+1)*(*n), double);
    dBeta(x, *m, *n, Bta);
    em_beta_mix(p, Bta, *m, *n,  *maxit, *eps, llik, convergence, delta);
    R_Free(Bta);
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
    pBta = R_Calloc((m+1)*N, double);
    fp = R_Calloc(N, double);
    pnu = R_Calloc(m+1, double);
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
    R_Free(pBta); R_Free(fp); R_Free(pnu);
}
/* end function em_beta_mix_group */

void mable_em_group(int *m, int *n, int *N, double *p, double *t, int *maxit,
    double *eps, double *llik, int *convergence, double *delta){
    double *dBta;
    dBta = R_Calloc((*m+1)*(*N), double);
    cpBeta(t, *m, *N, dBta);
    em_beta_mix_group(p, dBta, *N, *m, n, *maxit, *eps, llik, convergence, delta);
    R_Free(dBta);
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
  phat = R_Calloc(lp, double);
  diverge = R_Calloc(1, int);
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  llik = R_Calloc(1, double);
  Bta = R_Calloc((M[1]+1)*(*n), double);
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
    em_beta_mix(p, Bta, m, *n, *maxit, eps[0], llik, diverge, delta);
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
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  optim[0]=cp0+M[0];
  if(m==M[1] && pval[i-1]>*level){
    convergence[0]+=1; 
    if(*progress==1) warning("\nThe maximum candidate degree has been reached. \nA model degree m=%d with p-value,  %f, of the change-point is returned.\n", optim[0], pval[i-2]);
//        warning("The maximum candidate degree has been reached \n with a p-value of the change-point %f.\n", res[0]);
  }
  M[1]=m;
  m=optim[0];
  for(j=0;j<=m;j++) p[j]=phat[tmp+j];
  R_Free(phat); R_Free(Bta); R_Free(llik); 
  R_Free(cp); R_Free(res);
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
  cp=R_Calloc(1, int);
  diverge=R_Calloc(1, int);
  dBta = R_Calloc((M[1]+1)*(*N), double);
  phat = R_Calloc(lp, double);
  res = R_Calloc(2, double);
  ttl = (double)((k+2)*(k+1));
  mlik = R_Calloc(1, double);
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
    if(*progress==1) warning("\n The maximum candidate degree has been reached: \nm=%d with a p-value of the change-point %f.\n", m, res[0]);}
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  optim[0]=cp0+M[0];
  //    tmp=cp[0]*(M[0]*2+(cp[0]+1))/2;
  //    optim[0]=cp[0]+M[0];
  m=optim[0];
  for(j=0;j<=m;j++) p[j]=phat[tmp+j];
  R_Free(dBta); R_Free(phat); R_Free(diverge); R_Free(cp); R_Free(res);
  R_Free(mlik);
}
/*//////////////////////////////////////////////////////////*/
/*//////////////////////////////////////////*/
/*   Bernstein Polynomial Approximation     */
/*  Returns n-dim row vector, n=length(u)   */
/*//////////////////////////////////////////*/
void mable_approx(double *u, double *p, int *m, int *n, int *cdf){
    int i, j;
    double *Bta, tmp;
    Bta = R_Calloc(*n*(*m+1), double);
    if(*cdf==0) dBeta(u, *m, *n, Bta);
    if(*cdf==1) pBeta(u, *m, *n, Bta);
    for(j=0;j<*n;j++){
       tmp=.0;
       for(i=0;i<=*m;i++){
          tmp += Bta[j+*n*i]*p[i];
       }
       u[j]=tmp;
    }
    R_Free(Bta);
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
    SEXP args, res, tmp;
    int i, j, m, pc = 0;
    double y, *z;
    z = R_Calloc(n, double);
    MableStruct MS = (MableStruct) ex;
    m = MS->m;
    j = MS->j;
    y = MS->y;
    PROTECT(args = allocVector(REALSXP, n));
    for(i = 0; i < n; i++) REAL(args)[i] = y - x[i];

    PROTECT(tmp = lang2(MS->f, args));
    PROTECT(res = eval(tmp, MS->env));

    if(length(res) != n)
	     error("evaluation of function gave a result of wrong length");
    if(TYPEOF(res) == INTSXP) {
	     PROTECT(res = coerceVector(res, REALSXP));
       pc = 1;
	     //res = coerceVector(res, REALSXP);
       //UNPROTECT(1); /* uprotect the old res */
       //PROTECT(res);
    } 
    else if(TYPEOF(res) != REALSXP)
	   error("evaluation of error density gave a result of wrong type");
    for(i = 0; i < n; i++){
	   z[i] = REAL(res)[i];
	   //x[i] = z[i]*dbeta(x[i], j+1,  m-j+1, FALSE);
	   x[i] = z[i]*(m+1)*dbinom_raw(j, m, x[i],1-x[i], FALSE);
	   if(!R_FINITE(x[i]))
	       error("non-finite error denity value");
    }
    R_Free(z);
    UNPROTECT(3+pc);
    //UNPROTECT(3);
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
    gBeta = R_Calloc((m+1)*n, double);
    p_gBeta = R_Calloc((m+1)*n, double);
    fp = R_Calloc(n, double);
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
    R_Free(gBeta);
    R_Free(p_gBeta);
    R_Free(fp);
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
    llik = R_Calloc(1, double);
    lk = R_Calloc(nm+1, double);
    lr = R_Calloc(nm, double);
    pval = R_Calloc(nm+1, double);
    bic = R_Calloc(nm+1, double);
    chpts = R_Calloc(nm+1, int); 
    cp = R_Calloc(1, int); 
    res = R_Calloc(1, double); 
    lp = M[0]*(nm+1)+(nm+1)*(nm+2)/2;
    p = R_Calloc(M[1]+1, double); // R_Calloc(2*M[1]+2, double);
    phat = R_Calloc(lp, double);  
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
    R_Free(llik); 
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
    R_Free(lk); R_Free(lr); R_Free(phat); R_Free(p);
    R_Free(pval); R_Free(bic); R_Free(chpts); R_Free(cp); R_Free(res); 
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
    ex = R_Calloc(7, double);
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
    R_Free(ex);
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
    gam=R_Calloc(mp1*kp1*n, double);
    pi=R_Calloc(v+1, double);
    p=R_Calloc(mp1, double);
    q=R_Calloc(kp1, double);
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
    R_Free(gam);
    R_Free(pi);
    R_Free(p);
    R_Free(q);
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
    ex = R_Calloc(8+v+m+k, double);
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
    R_Free(ex); 
    return result;
} 
// EM for deconvolution with unknown error density
void mablem_decon(double *gam, int n, double *interval, int m, int k, double *lk,
        double *p, double *q, int mConstr, double ybar, double eps, int maxit, 
        double eps_nt, int maxit_nt){
    int it=0, i, j, l, mp1=m+1, kp1=k+1;
    double c=interval[0], d=interval[1], eta=0.0, zeta=0.0, del=0.0;
    double *Am, *Ck, *psi, *spg, *sqg, *p1, *q1;
    Am = R_Calloc(mp1, double);//??
    Ck = R_Calloc(kp1, double);//??
    p1 = R_Calloc(mp1, double);
    q1 = R_Calloc(kp1, double);
    psi = R_Calloc(n, double);
    spg = R_Calloc(n*kp1, double);
    sqg = R_Calloc(n*mp1, double);
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
    R_Free(Am);
    R_Free(Ck);
    R_Free(p1);
    R_Free(q1);
    R_Free(psi);
    R_Free(spg);
    R_Free(sqg);
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
    gam = R_Calloc(n*mp1*kp1, double);
    llik = R_Calloc(1, double);
    np=R_Calloc(nr*nc, int);
    nu_d=R_Calloc(2, int);
    nu_aic=R_Calloc(2, int);
    nu_bic=R_Calloc(2, int);
    lk=R_Calloc(nr*nc, double);
    D=R_Calloc(nr*nc, double);
    aic=R_Calloc(nr*nc, double);
    bic=R_Calloc(nr*nc, double);
    p=R_Calloc(mp1, double);
    q=R_Calloc(kp1, double);
    p_d=R_Calloc(mp1, double);
    p_aic=R_Calloc(mp1, double);
    p_bic=R_Calloc(mp1, double);
    q_d=R_Calloc(kp1, double);
    q_aic=R_Calloc(kp1, double);
    q_bic=R_Calloc(kp1, double);
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
    R_Free(np);
    R_Free(gam);
    R_Free(llik);
    R_Free(p);
    R_Free(q);
    R_Free(nu_d);
    R_Free(nu_aic);
    R_Free(nu_bic);
    R_Free(p_d);
    R_Free(p_aic);
    R_Free(p_bic);
    R_Free(q_d);
    R_Free(q_aic);
    R_Free(q_bic);
    R_Free(lk);
    R_Free(D);
    R_Free(aic);
    R_Free(bic);
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
/* Multivar_Beta: Returns a n x K matrix, K=(m1+1)...(md+1) */
/* km[0]=1, km[1]=m1+1,km[2]=(m1+1)(m2+1),...,km[d]=(m1+1)...(md+1)=K,*/
/*  Column-major layout  */
/* beta(x1j, i1+1, m1+1-i1)...beta(xdj, id+1, md+1-id), 0<=ik<=mk, 1<=k<=d,  
         the (i1+km[1]*i2+...+km[d-1]*id,j)th element */
void MV_dBeta(double *x, int *m, int n, int d, int *km, double *dBta){
  int i, j, jj, k, r, K, it;
  K=km[d];
  for(j=0; j<n; j++){
    for(it=0;it<K;it++){
      dBta[j+it*n]=1.0;
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        dBta[j+it*n]*=dbeta(x[j+n*k], i+1, m[k]+1-i, FALSE);
        r = jj;
        //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
      }
      dBta[j+it*n]*=dbeta(x[j], r+1, m[0]+1-r, FALSE); 
    }
  }
}
void MVdBeta_One_Obs(double *x, int *m, int j, int n, int d, int *km, double *dBta){
  int i, jj, k, r, K, it;
  K=km[d];
  
  for(it=0;it<K;it++){
    dBta[it]=1.0;
    r = it;
    for(k=d-1; k>0; k--){
      jj = r%km[k];
      i = (r-jj)/km[k];
      dBta[it]*=dbeta(x[j+n*k], i+1, m[k]+1-i, FALSE);
      r = jj;
      //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
    }
    dBta[it]*=dbeta(x[j], r+1, m[0]+1-r, FALSE); 
  }

}
/*//////////////////////////////////////////////////////////////////*/
/* Multivariate Beta cdfs                                           */
/* Multivar_Beta_CDF: Returns a K x n matrix, K=(m1+1)...(md+1)     */
/* km[1]=m1+1,km[2]=(m1+1)(m2+1),...,km[d]=(m1+1)...(md+1)=K,       */
/* pbeta(x1j, i1+1, m1+1-i1)...pbeta(xdj, id+1, md+1-id),           */
/* 0<=ik<=mk, 1<=k<=d, the (i1+km[1]*i2+...+km[d-1]*id,j)th element */
/*//////////////////////////////////////////////////////////////////*/
void MV_pBeta(double *x, int *m, int n, int d, int *km, 
         double *pBta) {
  int i, j, jj, k, r, K, it;
  K=km[d];
  for(j=0; j<n; j++){
    for(it=0;it<K;it++){
      pBta[j+it*n]=1.0;
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        pBta[j+it*n]*=pbeta(x[j+n*k], i+1, m[k]+1-i, TRUE, FALSE);
        r = jj;
      }
      pBta[j+it*n]*=pbeta(x[j], r+1, m[0]+1-r, TRUE, FALSE); 
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
    for(i=0; i<K; i++) fx += p[i]*Bta[j+i*n];
    lik *= fx;
  }
  //Rprintf("  lik = %g\n", lik);
  return log(lik);
}

/*////////////////////////////////////////////////////////*/
/*     Log-Bernstein-Likelihood of multivariate sample    */
/*     reparameterized p=alpha^2/||alpha||^2, where       */
/*      alpha=alpha(i1,...,id), 0<=ij<=mj, j=1,...,d,     */
/*         an array of dim=m+1=(m1+1,...,md+1),           */
/*    x: n x d matrix, sample from d-dimensional          */
/*          distribution F with support [0,1]^d           */
/*////////////////////////////////////////////////////////*/
double loglik_alpha(double *alpha, int K, int n, double *dBta){
  int i,j;
  double tmp=1.0, v;
  //Rprintf("  n=%d\n", n);
  for(j=0; j<n; j++){
    v = 0.0;
    for(i=0; i<K; i++){
        v += dBta[j+i*n]*alpha[i]*alpha[i];
    }
    tmp *= v;  
  }
  v = log(tmp);
  tmp=0.0;
  for(i=0; i<K; i++) tmp+=alpha[i]*alpha[i];
  v -= n*log(tmp);
  return log(v);
}


/*////////////////////////////////////////////////////////*/
/*     Log-Bernstein-Likelihood of multivariate sample    */
/*     reparameterized p=alpha^2/||alpha||^2, where       */
/*      alpha=alpha(i1,...,id), 0<=ij<=mj, j=1,...,d,     */
/*         an array of dim=m+1=(m1+1,...,md+1),           */
/*    x: n x d matrix, sample from d-dimensional          */
/*          distribution F with support [0,1]^d           */
/*  also return "grad" and "hess=grad'", so that          */
/*    2 x alpha x grad = the gradient of log-likelihood   */
/*    2grad+2alpha*hess = diagonoal of Hessian matrix     */
/*////////////////////////////////////////////////////////*/
void log_blik_alpha(double *alpha, int K, int n, double *dBta, 
            double *llik, double *grad, double *hess){
  int i,j;
  double anorm=0.0, fx, tmp;
  //Rprintf("  n=%d\n", n);
  for(i=0; i<K; i++) anorm+=alpha[i]*alpha[i];
  llik[0] = 1.0;
  for(i=0; i<K; i++){
      grad[i]=0.0;
      hess[i]=0.0;
  }
  for(j=0; j<n; j++){
    fx = 0.0;
    for(i=0; i<K; i++){
        fx += dBta[j+i*n]*alpha[i]*alpha[i];
    }
    llik[0] *= fx;  
    for(i=0; i<K; i++){
        tmp=dBta[j+i*n]/fx;
        grad[i]+=tmp;
        hess[i]+=2.0*alpha[i]*tmp*tmp;
    }
  }
  llik[0] = log(llik[0])-n*log(anorm);
  tmp=n/anorm;
  for(i=0; i<K; i++){
     grad[i]-=tmp;
     hess[i]-=2.0*alpha[i]*tmp/anorm;
  }
  //Rprintf("  lik = %g\n", lik);
  //return log(lik);
}


/* Utilities for Quasi-Newton method for coordinate-ascent method */
// This algorithm is SLOW!!!???
/* Extended version of struct opt_struct */
typedef struct mable_mvar_struct
{
    int K, n, k;
    double *alpha,  *dBeta; 
    //double anorm_k; // ||alpha||^2-alpha[k]^2
} mable_mvar_struct, *MableMVARStruct;

typedef struct made_mvar_struct
{
    int K, n, k;
    double *alpha,  *pBeta, *Fn; 
} made_mvar_struct, *MadeMVARStruct;

// minus log-likelihood to be called by vmmin()
// with par = alpha[k], npar = 1;  

static double mll_alphak(int npar, double *par, void *ex){
  int i, j, n, k, K;
  double val=1.0, anorm=0.0, fx;
  MableMVARStruct PS = (MableMVARStruct) ex;
  n=PS->n; k=PS->k; K=PS->K;
  for(i=0; i<k; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  anorm += par[0]*par[0];
  for(i=k+1; i<K; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  for(j=0;j<n;j++){
    fx = 0.0;
    for(i=0; i<k; i++){
        fx += (PS->dBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    fx += (PS->dBeta[j+k*n])*par[0]*par[0];
    for(i=k+1; i<K; i++){
        fx += (PS->dBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    val *= fx;  
  }
  val = n*log(anorm)-log(val);
  return val;
}


// derivative of minus log-likelihood to be called by vmmin()
// with par = alpha[k], npar = 1;  

static void dmll_alphak(int npar, double *par, double *df, void *ex)
{
  int i, j, n, k, K;
  double fx, anorm=0.0;
  MableMVARStruct PS = (MableMVARStruct) ex;
  n=PS->n; k=PS->k; K=PS->K;
  df[0]=0.0;
  
  for(i=0; i<k; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  anorm += par[0]*par[0];
  for(i=k+1; i<K; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  for(j=0;j<n;j++){
    fx = 0.0;
    for(i=0; i<k; i++){
        fx += (PS->dBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    fx += (PS->dBeta[j+k*n])*par[0]*par[0];
    for(i=k+1; i<K; i++){
        fx += (PS->dBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    df[0]-=(PS->dBeta[j+k*n])/fx;
  }
  df[0]+=n/anorm;
  df[0]*=2.0*par[0];
}


// l_2 distance between empirical cdf and Bernstein model to be called by vmmin()
// with par = alpha[k], npar = 1;  

static double D_alphak(int npar, double *par, void *ex){
  int i, j, n, k, K;
  double val=0.0, anorm=0.0, Fm;
  MadeMVARStruct PS = (MadeMVARStruct) ex;
  n=PS->n; k=PS->k; K=PS->K;
  for(i=0; i<k; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  anorm += par[0]*par[0];
  for(i=k+1; i<K; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  if(fabs(anorm)<0.000001) error("zero alpha's\n.");
  for(j=0;j<n;j++){
    Fm = 0.0;
    for(i=0; i<k; i++){
        Fm += (PS->pBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    Fm += (PS->pBeta[j+k*n])*par[0]*par[0];
    for(i=k+1; i<K; i++){
        Fm += (PS->pBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    val += (Fm/anorm-(PS->Fn[j]))*(Fm/anorm-(PS->Fn[j]));  
  }
  return val;
}


// derivative of l_2 distance between empirical cdf and Bernstein model to be called by vmmin()
// with par = alpha[k], npar = 1;  

static void dD_alphak(int npar, double *par, double *df, void *ex)
{
  int i, j, n, k, K;
  double Fm, anorm=0.0;
  MadeMVARStruct PS = (MadeMVARStruct) ex;
  n=PS->n; k=PS->k; K=PS->K;
  df[0]=0.0;
  
  for(i=0; i<k; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  anorm += par[0]*par[0];
  for(i=k+1; i<K; i++){
      anorm += (PS->alpha[i])*(PS->alpha[i]);
  }
  for(j=0;j<n;j++){
    Fm = 0.0;
    for(i=0; i<k; i++){
        Fm += (PS->pBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    Fm += (PS->pBeta[j+k*n])*par[0]*par[0];
    for(i=k+1; i<K; i++){
        Fm += (PS->pBeta[j+i*n])*(PS->alpha[i])*(PS->alpha[i]);
    }
    df[0] += (Fm/anorm-(PS->Fn[j]))*((PS->pBeta[j+k*n])-Fm/anorm);  
  }
  df[0]=df[0]*4.0*par[0]/anorm;
}


/* end of utilities for quasi-Newton method */

/*////////////////////////////////////////////////////////*/
/*     Updating the k-th component of alpha for           */
/*     reparametrization p=alpha^2/||alpha||^2, where     */
/*      alpha=alpha(i1,...,id), 0<=ij<=mj, j=1,...,d,     */
/*         an array of dim=m+1=(m1+1,...,md+1),           */
/*    x: n x d matrix, sample from d-dimensional          */
/*          distribution F with support [0,1]^d           */
/*//////////////////////////////////////////////////////*/
void new_alpha_k(int k, double *alpha, int K, int n, double *dBta, 
            double *lk, int maxit, double eps){
  //int it=0;
  double *par, val = 0.0, abstol=-INFINITY, reltol=eps;
  int npar=1, *mask, trace=0, fncount = 0, grcount = 0, nREPORT=10;
  int ifail = 0;
  mask = R_Calloc(1, int);
  mask[0] = 1;
  MableMVARStruct PS;
  PS = (MableMVARStruct) R_alloc(1, sizeof(mable_mvar_struct));
  PS->K=K; PS->n=n; PS->k=k; 
  PS->alpha=alpha; PS->dBeta=dBta;  
  par = R_Calloc(1, double);
  par[0]=alpha[k];
  //Rprintf("  n=%d, eps = %g\n", n, eps);
  //Rprintf("  maxit=%d, eps = %g\n", maxit, eps);
  //Quasi-Newton
  vmmin(npar, par, &val, mll_alphak, dmll_alphak, maxit, trace, mask, abstol, 
	      reltol, nREPORT, (void *)PS, &fncount, &grcount, &ifail);
    lk[0]=-val;  
    alpha[k]=par[0];
  // Faster approach but not stable
  /* 
  alpha[k]=0.0;
  log_blik_alpha(alpha, K, n, dBta, lk, grad, hess);
  if(grad[k]<=0){
    alpha[k]=0.0;
  }
  else {
    alpha[k]=alpha_k;
    log_blik_alpha(alpha, K, n, dBta, lk, grad, hess);
    del=fabs(2*alpha[k]*grad[k]);
    ddll=2*(grad[k]+alpha[k]*hess[k]);
    while(it<maxit && del>eps && fabs(ddll)>tol){
        tmp=2*alpha[k]*grad[k]/ddll;
        del = fabs(tmp);
        alpha[k]-=tmp;
        log_blik_alpha(alpha, K, n, dBta, lk, grad, hess);
        del+=fabs(2*alpha[k]*grad[k]);
        ddll=2*(grad[k]+alpha[k]*hess[k]);
        it+=1;
        R_CheckUserInterrupt();
        //Rprintf("  grad[%d] = %g, hess[%d] = %g\n", k, grad[k], k, hess[k]);
    }
  }
  */
  R_Free(par); R_Free(mask); 
}

/*////////////////////////////////////////////////////////////////////*/
/*  Coordinate Descent optimization for mable with given degree m     */
/*  x = n-vector or n x d matrix of sample data in [0,1]^d            */
/*  m d-vector of model degree(s)                                     */
/*  alpha an K=prod(m+1)-vector so that the coefficients              */
/*       p_j=\alpha_j^2/\Vert\alpha\Vert^2,  j=0,1,\dots,K-1}.        */
/*////////////////////////////////////////////////////////////////////*/
void mable_m_cd(int *m, int *n, int *d, double *alpha, double *x, int *maxit,
      double *eps, double *lk){
  int k, it=0, K, *km, nn=n[0];
  double *alphat, *dBeta,  lkt, del=1.0,  tmp=.0;//*grad, *hess,
  km = R_Calloc(*d+1, int);
  km[0]=1;
  for(k=1;k<=*d;k++) km[k]=km[k-1]*(m[k-1]+1);
  K=km[*d];
  alphat = R_Calloc(K, double);
  dBeta = R_Calloc(nn*K, double);
  //grad = R_Calloc(K, double);
  //hess = R_Calloc(K, double);
  //Rprintf("  n=%d, d = %d, maxit= %d\n", *n, *d, *maxit);
  //for(k=0;k<K;k++) alpha[k]=1.0;
  //for(k=0;k<K;k++) Rprintf("  alpha[%d] = %g\n", k, alpha[k]);
  MV_dBeta(x, m, nn, *d, km, dBeta);
  //log_blik_alpha(alpha, K, nn, dBeta, lk, grad, hess);
  lk[0]=loglik_alpha(alpha, K, nn, dBeta);
  for(k=0;k<K;k++) alphat[k]=alpha[k]; 
  lkt=lk[0];
  while(it<*maxit && del>*eps){
    for(k=0;k<K;k++){
      new_alpha_k(k, alpha, K, nn, dBeta,  lk, *maxit, *eps);
    }
    del=fabs(lkt-lk[0]);
    for(k=0;k<K;k++){
        //del+=fabs(alpha[k]-alphat[k]);
        alphat[k]=alpha[k];
    }
    lkt=lk[0];
    it+=1;
    R_CheckUserInterrupt();
    //Rprintf("it=%d, del=%g\n", it, del);
  }
  for(k=0;k<K;k++){
    alpha[k]=alpha[k]*alpha[k]; 
    tmp+=alpha[k]; 
  } 
  for(k=0;k<K;k++){
    alpha[k]/=tmp;  // p=alpha
    //Rprintf("  alpha[%d] = %g\n", k, alpha[k]);
  }
  R_Free(km);
  R_Free(alphat); R_Free(dBeta); 
}

//  Minimum Approximate Distance Estimate
// Approximate Distance between empirical cdf and Bernstein model
double AD_alpha(int n, int K, double *alpha, double *pBeta, double *Fn){
  int i, j;
  double val=0.0, anorm=0.0, Fm;
  for(i=0; i<K; i++){
      anorm += alpha[i]*alpha[i];
  }
  for(j=0;j<n;j++){
    Fm = 0.0;
    for(i=0; i<K; i++){
        Fm += pBeta[j+i*n]*alpha[i]*alpha[i];
    }
    val += (Fm/anorm-Fn[j])*(Fm/anorm-Fn[j]);  
  }
  return val;
}

/*////////////////////////////////////////////////////////*/
/*     Updating the k-th component of alpha for           */
/*     reparametrization p=alpha^2/||alpha||^2, where     */
/*      alpha=alpha(i1,...,id), 0<=ij<=mj, j=1,...,d,     */
/*         an array of dim=m+1=(m1+1,...,md+1),           */
/*    x: n x d matrix, sample from d-dimensional          */
/*          distribution F with support [0,1]^d           */
/*////////////////////////////////////////////////////////*/
void update_alpha_k(int k, double *alpha, int K, int n, double *pBta, 
            double *Fn, double *Dn, int maxit, double eps){
  //int it=0;
  double *par, val = 0.0, abstol=-INFINITY, reltol=eps;
  int npar=1, *mask, trace=0, fncount = 0, grcount = 0, nREPORT=10;
  int ifail = 0;
  mask = R_Calloc(1, int);
  mask[0] = 1;
  MadeMVARStruct PS;
  PS = (MadeMVARStruct) R_alloc(1, sizeof(made_mvar_struct));
  PS->K=K; PS->n=n; PS->k=k; 
  PS->alpha=alpha; PS->pBeta=pBta; PS->Fn=Fn;  
  par = R_Calloc(1, double);
  par[0]=alpha[k];
  //Rprintf("  alpha[%d] = %g\n", k, alpha[k]);
  //Rprintf("  n=%d, eps = %g\n", n, eps);
  //Rprintf("  maxit=%d, eps = %g\n", maxit, eps);
  //Quasi-Newton
  vmmin(npar, par, &val, D_alphak, dD_alphak, maxit, trace, mask, abstol, 
	      reltol, nREPORT, (void *)PS, &fncount, &grcount, &ifail);
    Dn[0]=val;  
    alpha[k]=par[0];
  R_Free(par); R_Free(mask); 
}

/*////////////////////////////////////////////////////////////////////*/
/*  Coordinate Descent optimization for MADE with given degree m      */
/*  x = n-vector or n x d matrix of sample data in [0,1]^d            */
/*  m d-vector of model degree(s)                                     */
/*  alpha an K=prod(m+1)-vector so that the coefficients              */
/*       p_j=\alpha_j^2/\Vert\alpha\Vert^2,  j=0,1,\dots,K-1}.        */
/*////////////////////////////////////////////////////////////////////*/
void made_m_cd(int *m, int *n, int *d, double *alpha, double *x, double *Fn, 
      int *maxit, double *eps, double *Dn){
  int k, it=0, K, *km, nn=*n, dd=*d;
  double *alphat, *pBeta,  del=1.0,  tmp=.0, Dt; 
  km = R_Calloc(dd+1, int);
  Rprintf("maxit=%d, eps = %g, \n", *maxit, *eps);
  //Rprintf("n=%d, d = %d, \n", nn, dd);
  //for(k=0;k<nn;k++) Rprintf("  Fn[%d] = %g\n", k, Fn[k]);
  km[0]=1;
  for(k=1;k<=dd;k++) km[k]=km[k-1]*(m[k-1]+1);
  K=km[dd];
  alphat = R_Calloc(K, double);//??? necessary
  pBeta = R_Calloc(nn*K, double);
  MV_pBeta(x, m, nn, dd, km, pBeta);
  Dn[0]=AD_alpha(nn, K, alpha, pBeta, Fn);
  for(k=0;k<K;k++) alphat[k]=alpha[k]; 
  Dt=Dn[0];
  while(it<*maxit && del>*eps){
    for(k=0;k<K;k++){
      update_alpha_k(k, alpha, K, nn, pBeta, Fn, Dn, *maxit, *eps);
    }
    del=fabs(Dt-Dn[0]);
    for(k=0;k<K;k++){
        //del+=fabs(alpha[k]-alphat[k]);
        alphat[k]=alpha[k];
    }
    Dt=Dn[0];
    it+=1;
    R_CheckUserInterrupt();
    Rprintf("it=%d, del=%g\n", it, del);
  }
  for(k=0;k<K;k++){
    alpha[k]=alpha[k]*alpha[k]; 
    tmp+=alpha[k]; 
  } 
  for(k=0;k<K;k++){
    alpha[k]/=tmp;  // p=alpha
    //Rprintf("  alpha[%d] = %g\n", k, alpha[k]);
  }
  R_Free(km);
  R_Free(alphat); 
  R_Free(pBeta); 
}


/*/////////////////////////////////////////////////////*/
/* A slower version but takes less memory              */
/* EM Method for mixture of                            */
/* beta(x1, i1+1, m1+1-i1),...,beta(xd, id+1, md+1-id),*/ 
/*               0<=ik<=mk, 1<=k<=d                    */
/*/////////////////////////////////////////////////////*/
void em_mixmvbeta_lm(double *x, double *p, int *m, int *km,
      int n, int d, int K, int maxit, double eps, 
      double *llik, int progress, int *conv){
  int i, j, it;
  double del, llik_nu, *pb, *pt, fp=0.0, ttl;
  ttl=(double) maxit;
  conv[0] = 0;
  pb = R_Calloc(K, double);
  pt = R_Calloc(K, double);
  llik[0] = -n*log(n);
  del = 10.0;
  it = 0;
  while(del>eps && it<maxit){
    llik_nu=0.0;
    for(i=0; i<K; i++) pt[i]=0.0;
    for(j=0; j<n; j++){
      fp = 0.0;
      MVdBeta_One_Obs(x, m, j, n, d, km, pb);
      for(i=0; i<K; i++) {
        pb[i] = p[i]*pb[i];
        fp += pb[i];
      }
      for(i=0; i<K; i++) pt[i] += pb[i]/fp;
      llik_nu += log(fp);
    }
    for(i=0; i<K; i++){
      p[i] = pt[i]/(double)n;
    }
    del = fabs(llik[0]-llik_nu);
    it++;
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
  R_Free(pb);
  R_Free(pt);
}
/* end function em_mixmvbeta_lm */

/*/////////////////////////////////////////////////////*/
/* EM Method for mixture of                            */
/* beta(x1, i1+1, m1+1-i1),...,beta(xd, id+1, md+1-id),*/ 
/*               0<=ik<=mk, 1<=k<=d                    */
/*/////////////////////////////////////////////////////*/
void em_mixmvbeta(double *p, double *Bta, int *m, 
      int n, int d, int K, int maxit, double eps, 
      double *llik, int progress, int *conv){
  int i, j, it;
  double del, llik_nu, *fp;//*pBeta,  *pnu;
  double tmp, ttl;
  ttl=(double) maxit;
  conv[0] = 0;
  //pBeta = R_Calloc(K*n, double);//very BIG memory
  fp = R_Calloc(n, double);
  //pnu = R_Calloc(K, double);
  llik[0] = loglik_bern_multivar(p, K, n, Bta);
  del = 10.0;
  it = 1;
  while(del>eps && it<maxit){
    for(j=0; j<n; j++){
      fp[j] = 0.0;
      for(i=0; i<K; i++) {
        //pBeta[i+K*j] = p[i]*Bta[j+i*n];
        //fp[j] += pBeta[j+i*n];
        fp[j] += p[i]*Bta[j+i*n];
      }
    }
    for(i=0; i<K; i++){
      //p[i] = 0.0;
      //for(j=0; j<n; j++) p[i] += pBeta[j+i*n]/fp[j];
      //p[i] /= (double)n;
      tmp = 0.0; 
      for(j=0; j<n; j++) tmp += p[i]*Bta[j+i*n]/fp[j];
      p[i] = tmp/(double)n;
    }
    llik_nu = loglik_bern_multivar(p, K, n, Bta);
    del = fabs(llik[0]-llik_nu);
    it++;
    //for(i=0; i<K; i++) p[i] = pnu[i];
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
  //Free(pBeta);
  R_Free(fp);
  //Free(pnu);
}
/* end function em_beta_mix */
/*////////////////////////////////////////////////////////////*/
/* Calculate pt of degree mt=m+e_k from p of degree m,        */
/* where e_k is the k-th unit vector of standard basis of R^d */
/*   km[0] = 1,  km[1]=m[0]+1,  km[2]=(m[0]+1)(m[1]+1),...,   */
/*          km[d]=(m[0]+1)...(m[d-1]+1)=K,                    */
/* p(i[0], ..., i[d-1]), 0<=i[k]<=mk, 0<=k<d, are arranged by */
/* column-major order:   i[0]+km[1]*i[1]+...+km[d-1]*i[d-1]   */
/* input p, output pt, return it as pt                        */
/*////////////////////////////////////////////////////////////*/

void pm2pmpe_k(double *p, double *pt, int d, int *m, int *km, int k){
  int j, l, r, K, K1, it, *I;
  I = R_Calloc(d, int);
  K=km[d];
  K1 = km[d]*(m[k]+2)/(m[k]+1);
  for(j=0; j<K1; j++) pt[j] = 0.0;
  for(it=K-1; it>=0; it--){
    r = it;
    for(l=d-1; l>0; l--){
      j = r%km[l];
      I[l] = (r-j)/km[l];
      r = j;
    }
    I[0]=r;
    //Rprintf("i=%d, r=%d\n", i, r);
    j=0;
    r=0;
    for(l=0;l<d;l++){
        j+=(I[l]+(l==k))*km[l]*(m[k]+2*(l>k))/(m[k]+(l>k));
        r+=I[l]*km[l]*(m[k]+2*(l>k))/(m[k]+(l>k));
    }
    pt[j] += (I[k]+1)*p[it]/(m[k]+2.0);
    pt[r] += (m[k]+1-I[k])*p[it]/(m[k]+2.0);
  }
  // updating m and km
  m[k]+=1;
  for(j=0; j<=d; j++) km[j] = km[j]*(m[k]+2*(j>k))/(m[k]+(j>k));
  R_Free(I);  
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
  int K=km[d], it, i, jj, k, r;//, *I;
  for(i=0;i<=m[j];i++) pj[i]=0.0;
  //I=R_Calloc(d, int);
  it=0;
  while(it<K){
    r = it;
    for(k=d-1; k>0; k--){
      jj = r%km[k];
      i = (r-jj)/km[k];
      //I[k]=i;
      if(j==k) pj[i]+=p[it];
      r = jj;
      //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
    }
    //I[0]=r;
    //pj[I[j]]+=p[it];
    if(j==0) pj[r]+=p[it];
    it++;
  }
  //Free(I);
}

/*############################################################*/
// MABLE of multivariate density with model degrees given 
//  or preselected based on marginal data 
/*############################################################*/
void mable_m_mvar(int *M, int *n, int *d, double *phat, double *x, int *maxit,
      double *eps, double *lk, int *progress, int *conv, int *hd){
  int i, prgrs, *km, K;
  double *Bta;
  km = R_Calloc(*d+1, int);
  km[0] = 1;
  // Rprintf("km[%d]=%d\n",0,km[0]);
  for(i=1; i<=*d; i++){
    km[i]=km[i-1]*(M[i-1]+1);
    // Rprintf("km[%d]=%d\n",i,km[i]);
  }
  K=km[*d]; 
  //Rprintf("progrss=%d\n",*progress);
  if(*progress!=0)
  Rprintf("\n Mable fit of multivariate data. This may take several minutes.\n\n");
  if(*progress==1) prgrs=1;
  else prgrs=0;
  for(i=0;i<K;i++) phat[i]=1.0/(double) K;
  if(*hd==0){
    Bta = R_Calloc(*n*K, double);// very BIG memory
    MV_dBeta(x, M, *n, *d, km, Bta);
    em_mixmvbeta(phat, Bta, M, *n, *d, K, *maxit, *eps, lk, prgrs, conv);
    R_Free(Bta);
  }
  else{
    em_mixmvbeta_lm(x, phat, M, km, *n, *d, K, *maxit, *eps, lk, prgrs, conv);
  }
  R_Free(km);
}
// end of mablem_mvar 

/*############################################################*/
// MABLE of multivariate density with model degrees 
//  or to be selected by the method of change-point in [m, M]
/*############################################################*/
void mable_mvar(int *m, int *M, int *n, int *d, double *phat, double *x,  
      int *maxit,  double *eps, double *level, double *pval, double *lk, 
      double *lr, int *chpts, int *progress, int *conv,  int *hd){
  int i, j, k, prgrs, *km, K, maxK=1, *cp, max_k=0;
  int cp0=0, *Mhat, Khat=0, dd=d[0], *mt, *kmt; 
  double *p, *pt, *Phat, *llik, maxlik, pct=0.0, ttl,  *res;
  //double tini=1e-4, sump;//, Dt, tini=1e-4;
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);

  km = R_Calloc(dd+1, int);
  Mhat = R_Calloc(dd, int);
  mt = R_Calloc(dd, int);
  kmt = R_Calloc(dd+1, int);
  llik = R_Calloc(1, double);
  km[0] = 1;
  // Rprintf("km[%d]=%d\n",0,km[0]);
  for(i=1; i<=dd; i++){
    maxK*=(M[i-1]+1);
    km[i]=km[i-1]*(m[i-1]+1);
    max_k += M[i-1]-m[i-1];
    // Rprintf("km[%d]=%d\n",i,km[i]);
  }
  //max_k-=(dd+1);
  Phat = R_Calloc(maxK, double);  
  pt = R_Calloc(maxK, double);  
  p = R_Calloc(maxK, double);  
  K=km[dd]; 
  //Rprintf("progrss=%d\n",*progress);
  if(*progress!=0)
    Rprintf("\n Mable fit of multivariate data. This may take several minutes.\n\n");
  
  prgrs=0;
//  Rprintf("level=%g\n",*level);
  ttl= (double) max_k;
  //p = Calloc(K, double);  
  
  for(j=0;j<K;j++) p[j]=1.0/(double) K;
  mable_m_mvar(m, n, &dd, p, x, maxit, eps, llik, &prgrs, conv, hd);
//error("mablem:overflow.\n");
  lk[0]=llik[0];
  maxlik=lk[0];
  for(j=0;j<K;j++){
    Phat[j]=p[j];
    phat[j]=p[j];
    pt[j]=p[j];
//    Rprintf("pt[%d]=%g\n",j,pt[j]);
  }
  k = 1; 
  while(k<=max_k && pval[k]>*level){
      for(i=0; i<dd; i++){
      kmt[0]=1;
         for(j=0;j<dd;j++) {
            mt[j]=m[j];
            kmt[j+1]=kmt[j]*(m[j]+1);
         }
          pm2pmpe_k(p, pt, dd, mt, kmt, i);
          K=kmt[dd];
          mable_m_mvar(mt, n, &dd, pt, x, maxit, eps, llik, &prgrs, conv, hd);
          //Rprintf("k=%d, max_deg=%d, llik = %g\n", k, max_deg, llik[0]);
          if(i==0 || llik[0]>=maxlik){
              lk[k]=llik[0];
              maxlik=lk[k];
              Khat=K;
              for(j=0; j<dd; j++) Mhat[j]=mt[j];
              for(j=0; j<K; j++) Phat[j]=pt[j];
          }
      }
      K=Khat;
      for(i=0; i<dd; i++) m[i]=Mhat[i];
      for(j=0; j<K; j++) p[j]=Phat[j];
      if(k<=4){
        pval[k]=1.0;
        chpts[k]=0;
      }
      else{            
        cp[0]=k;
        chpt_exp(lk, lr, res, cp);
        pval[k]=res[0];
        chpts[k]=cp[0];
        
        if(chpts[k]!=cp0){
            cp0=chpts[k];
            for(i=0; i<dd; i++) M[i]=m[i];
            for(j=0;j<K;j++) phat[j]=p[j];             
        }
      } 
//      Rprintf("\n chpts[%d]=%d, pval[%d]=%g,  lk[%d]=%g\n",k, chpts[k], k, pval[k], k, lk[k]);
      pct += 1.0;
      if(*progress!=0) ProgressBar(fmin2(1.0,pct/ttl),"");
      R_CheckUserInterrupt();
      k++;
  }
  d[0]=k-1;      
       
  R_Free(p); R_Free(Phat); R_Free(pt); 
  R_Free(Mhat);
  R_Free(mt);
  R_Free(kmt);
  R_Free(km);
  R_Free(llik);
  R_Free(cp);
  R_Free(res);
}
// end of mable_mvar


/*///////////////////////////////////////////////////////////////*/
/*  Coordinate Descent optimization for mable with degree m      */
/*   selected between [m,M] by the method of change-point        */
/*  x = n-vector or n x d matrix of sample data in [0,1]^d       */
/*  m d-vector of model degree(s)                                */
/*  alpha an K=prod(m+1)-vector so that the coefficients         */
/*       p_j=\alpha_j^2/\Vert\alpha\Vert^2,  j=0,1,\dots,K-1}.   */
/*///////////////////////////////////////////////////////////////*/
 
void mable_cd(int *m, int *M, int *n, int *d, double *phat, double *x,  
      int *maxit, double *eps, double *level, double *pval, double *lk, 
      double *lr, int *chpts, int *progress){
  int i, j, k, *km, K, maxK=1, *cp, max_k=0;
  int cp0=0, *Mhat, Khat=0, dd=d[0], *mt, *kmt; 
  double *p, *pt, *Phat, maxlik, pct=0.0, ttl, *llik, *res;
  //double tini=1e-4, sump;//, Dt, tini=1e-4;

  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  llik = R_Calloc(1, double);

  km = R_Calloc(dd+1, int);
  Mhat = R_Calloc(dd, int);
  kmt = R_Calloc(dd+1, int);
  mt = R_Calloc(dd, int);
  km[0] = 1;
  // Rprintf("km[%d]=%d\n",0,km[0]);
  for(i=1; i<=dd; i++){
    maxK*=M[i-1]+1;
    km[i]=km[i-1]*(m[i-1]+1);
    max_k += M[i-1]-m[i-1];
    // Rprintf("km[%d]=%d\n",i,km[i]);
  }
  //max_k-=(dd+1);
//  Rprintf("maxK=%d\n",maxK);
  Phat = R_Calloc(maxK, double);  
  pt = R_Calloc(maxK, double);  
  p = R_Calloc(maxK, double);  
  K=km[dd]; 
//  Rprintf("K=%d\n",K);
  //Rprintf("progress=%d\n",*progress);
  if(*progress!=0)
    Rprintf("\n Mable fit of multivariate data. This may take several minutes.\n\n");
  
  //Rprintf("level=%g\n",*level);
  ttl= (double) max_k;
  for(j=0;j<K;j++) p[j]=1.0;
  mable_m_cd(m, n, &dd, p, x, maxit, eps, llik);
  lk[0]=llik[0];
  maxlik=lk[0];
  for(j=0;j<K;j++){
    Phat[j]=p[j];
    phat[j]=p[j];
    pt[j]=p[j];
  }
  k = 1; 
  while(k<=max_k && pval[k]>*level){
      for(i=0; i<dd; i++){
      kmt[0]=1;
         for(j=0;j<dd;j++) {
            mt[j]=m[j];
            kmt[j+1]=kmt[j]*(m[j]+1);
         }
          pm2pmpe_k(p, pt, dd, mt, kmt, i);
          K=kmt[dd];
          for(j=0; j<K; j++) pt[j]=sqrt(pt[j]);
          mable_m_cd(mt, n, &dd, pt, x, maxit, eps, llik);
          //Rprintf("k=%d, max_kdeg=%d, llik = %g\n", k, max_k, llik[0]);
          if(i==0 || llik[0]>=maxlik){
              lk[k]=llik[0];
              maxlik=lk[k];
              Khat=K;
              for(j=0; j<dd; j++) Mhat[j]=mt[j];
              for(j=0; j<K; j++) Phat[j]=pt[j];
          }
      }
    
      for(i=0; i<dd; i++) m[i]=Mhat[i];
      for(j=0; j<Khat; j++) p[j]=Phat[j];
      if(k<=4){
        pval[k]=1.0;
        chpts[k]=0;
      }
      else{            
        cp[0]=k;
        chpt_exp(lk, lr, res, cp);
        pval[k]=res[0];
        chpts[k]=cp[0];
        if(chpts[k]!=cp0){
            cp0=chpts[k];
            for(i=0; i<dd; i++) M[i]=m[i];
            for(j=0;j<Khat;j++) phat[j]=p[j];             
        }
      } 
      //Rprintf("\n chpts[%d]=%d, pval[%d]=%g,  lk[%d]=%g\n",k, chpts[k], k, pval[k], k, lk[k]);
      pct += 1.0;
      if(*progress!=0) ProgressBar(fmin2(1.0,pct/ttl),"");
      R_CheckUserInterrupt();
      k+=1;
  }
  d[0]=k-1; 
  R_Free(p); R_Free(Phat); R_Free(pt); 
  R_Free(Mhat);
  R_Free(km);
  R_Free(kmt);
  R_Free(mt);
  R_Free(cp);
  R_Free(res);
  R_Free(llik);
}
// end of mable_cd


/*///////////////////////////////////////////////////////////////*/
/*  Coordinate Descent optimization for MADE with degree m       */
/*   selected between [m,M] by the method of change-point        */
/*  x = n-vector or n x d matrix of sample data in [0,1]^d       */
/*  m d-vector of model degree(s)                                */
/*  alpha an K=prod(m+1)-vector so that the coefficients         */
/*       p_j=\alpha_j^2/\Vert\alpha\Vert^2,  j=0,1,\dots,K-1}.   */
/*///////////////////////////////////////////////////////////////*/
 
void made_cd(int *m, int *M, int *n, int *d, double *phat, double *x, double *Fn,  
      int *maxit, double *eps, double *level, double *pval, double *Dn, 
      double *lr, int *chpts, int *progress){
  int i, j, k, *km, K, maxK=1, *cp, max_k=0;
  int cp0=0, *Mhat, Khat=0, dd=d[0], *mt, *kmt; 
  double *p, *pt, *Phat, minD, pct=0.0, ttl, *ad, *res;
  //double tini=1e-4, sump;//, Dt, tini=1e-4;

  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  ad = R_Calloc(1, double);

  km = R_Calloc(dd+1, int);
  Mhat = R_Calloc(dd, int);
  kmt = R_Calloc(dd+1, int);
  mt = R_Calloc(dd, int);
  km[0] = 1;
  // Rprintf("km[%d]=%d\n",0,km[0]);
  for(i=1; i<=dd; i++){
    maxK*=M[i-1]+1;
    km[i]=km[i-1]*(m[i-1]+1);
    max_k += M[i-1]-m[i-1];
    // Rprintf("km[%d]=%d\n",i,km[i]);
  }
  //max_k-=(dd+1);
//  Rprintf("maxK=%d\n",maxK);
  Phat = R_Calloc(maxK, double);  
  pt = R_Calloc(maxK, double);  
  p = R_Calloc(maxK, double);  
  K=km[dd]; 
//  Rprintf("K=%d\n",K);
  //Rprintf("progress=%d\n",*progress);
  if(*progress!=0)
    Rprintf("\n Mable fit of multivariate data. This may take several minutes.\n\n");
  
  Rprintf("made_cd: level=%g\n",*level);
  ttl= (double) max_k;
  for(j=0;j<K;j++) p[j]=1.0;
  made_m_cd(m, n, &dd, p, x, Fn, maxit, eps, ad);
  Dn[0]=-log(ad[0]);
  minD=ad[0];
  for(j=0;j<K;j++){
    Phat[j]=p[j];
    phat[j]=p[j];
    pt[j]=p[j];
  }
  k = 1; 
  while(k<=max_k && pval[k]>*level){
      for(i=0; i<dd; i++){
      kmt[0]=1;
         for(j=0;j<dd;j++) {
            mt[j]=m[j];
            kmt[j+1]=kmt[j]*(m[j]+1);
         }
          pm2pmpe_k(p, pt, dd, mt, kmt, i);
          K=kmt[dd];
          for(j=0; j<K; j++) pt[j]=sqrt(pt[j]);
          made_m_cd(mt, n, &dd, pt, x, Fn, maxit, eps, ad);
          //Rprintf("k=%d, max_kdeg=%d, ad = %g\n", k, max_k, ad[0]);
          if(i==0 || ad[0]<=minD){
              Dn[k]=-log(ad[0]);
              minD=ad[0];
              Khat=K;
              for(j=0; j<dd; j++) Mhat[j]=mt[j];
              for(j=0; j<K; j++) Phat[j]=pt[j];
          }
      }
    
      for(i=0; i<dd; i++) m[i]=Mhat[i];
      for(j=0; j<Khat; j++) p[j]=Phat[j];
      if(k<=4){
        pval[k]=1.0;
        chpts[k]=0;
      }
      else{            
        cp[0]=k;
        chpt_exp(Dn, lr, res, cp);
        pval[k]=res[0];
        chpts[k]=cp[0];
        if(chpts[k]!=cp0){
            cp0=chpts[k];
            for(i=0; i<dd; i++) M[i]=m[i];
            for(j=0;j<Khat;j++) phat[j]=p[j];             
        }
      } 
      //Rprintf("\n chpts[%d]=%d, pval[%d]=%g,  Dn[%d]=%g\n",k, chpts[k], k, pval[k], k, Dn[k]);
      pct += 1.0;
      if(*progress!=0) ProgressBar(fmin2(1.0,pct/ttl),"");
      R_CheckUserInterrupt();
      k+=1;
  }
  d[0]=k-1; 
  R_Free(p); R_Free(Phat); R_Free(pt); 
  R_Free(Mhat);
  R_Free(km);
  R_Free(kmt);
  R_Free(mt);
  R_Free(cp);
  R_Free(res);
  R_Free(ad);
}
// end of made_cd


/* t=(t1,...,td) is in [0,1]^d */
void mable_mvdf(int *d, int *m, int *km, int *n, double *t, double *p, 
        double *mvdf, int *density){
  int i, j, K;
  double *tmp;
  K=km[*d];
  tmp = R_Calloc((*n)*K, double);
  if(*density==0) MV_pBeta(t, m, *n, *d, km, tmp);
  if(*density==1) MV_dBeta(t, m, *n, *d, km, tmp);
  for(i=0;i<*n;i++){
    mvdf[i]=0.0;
    for(j=0;j<K;j++) mvdf[i]+=p[j]*tmp[i+(*n)*j];
  }
  R_Free(tmp);
}
// end of mable-multivar.c
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*    Maximum Approximate Bernstein likelihood Estimate   */
/*              of Copula Density Function                */
/*                                                        */
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/*                                                            */
/*  Reference:                                                */
/*  Guan, Z.,(???) Bernstein Polynomial Model for             */
/*     Nonparametric Estimate of Copula Density,              */
/*////////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/* Multivariate Bernstein base polynomials */
/* Multivar_Beta: Returns a K x n matrix, K=(m1+1)...(md+1)   */
/* km[0]=1, km[1]=m1+1,km[2]=(m1+1)(m2+1),...,km[d]=(m1+1)...(md+1)=K,*/
/*  Column-major layout  */
/* beta(x1j, i1+1, m1+1-i1)...beta(xdj, id+1, md+1-id), 0<=ik<=mk, 1<=k<=d,  
         the (i1+km[1]*i2+...+km[d-1]*id,j)th element */
/*////////////////////////////////////////////////////////////*/
/* Multivariate Bernstein base polynomials                    */
/* Multivar_Beta: Returns a K x n matrix, K=(m1+1)...(md+1)   */
/* km[0]=1, km[1]=m1+1,km[2]=(m1+1)(m2+1),...,km[d]=(m1+1)...(md+1)=K,*/
/*  Column-major layout  */
/* beta(x1j, i1+1, m1+1-i1)...beta(xdj, id+1, md+1-id), 0<=ik<=mk, 1<=k<=d,  
         the (i1+km[1]*i2+...+km[d-1]*id,j)th element */
void dBeta_copula(double *x, int *m, int n, int d, int *km, double *PdBta, 
        double *SdBta){
  int i, j, jj, k, r, K, it;
  double tmp;
  K=km[d];
  for(j=0; j<n; j++){
    for(it=0;it<K;it++){
      PdBta[it+K*j]=1.0;
      SdBta[it+K*j]=0.0;
      r = it;
      for(k=d-1; k>0; k--){
        jj = r%km[k];
        i = (r-jj)/km[k];
        tmp=dbeta(x[j+n*k], i+1, m[k]+1-i, FALSE);
        PdBta[it+K*j]*=tmp;
        SdBta[it+K*j]+=tmp;
        r = jj;
        //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
      }
      tmp=dbeta(x[j], r+1, m[0]+1-r, FALSE); 
      PdBta[it+K*j]*=tmp; 
      SdBta[it+K*j]+=tmp; 
    }
  }
}
void dBeta_copula_one_obs(double *x, int *m, int j, int n, int d, int *km, 
        double *PdBta, double *SdBta){
  int i, jj, k, r, K, it;
  double tmp;
  K=km[d];
  
  for(it=0;it<K;it++){
    PdBta[it]=1.0;
    SdBta[it]=0.0;
    r = it;
    for(k=d-1; k>0; k--){
      jj = r%km[k];
      i = (r-jj)/km[k];
      tmp=dbeta(x[j+n*k], i+1, m[k]+1-i, FALSE);
      PdBta[it]*=tmp;
      SdBta[it]+=tmp;
      r = jj;
      //Rprintf("it=%d, k=%d, i=%d\n",  it, k, i);
    }
    tmp=dbeta(x[j], r+1, m[0]+1-r, FALSE); 
    PdBta[it]*=tmp; 
    SdBta[it]+=tmp; 
  }
}

/*/////////////////////////////////////////////////////*/
/* A slower version but takes less memory              */
/* EM Method for mixture of high dimimensional         */
/* beta(x1, i1+1, m1+1-i1),...,beta(xd, id+1, md+1-id),*/ 
/*               0<=ik<=mk, 1<=k<=d                    */
/*   with uniform marginal constraints                 */
/*/////////////////////////////////////////////////////*/
void em_copula_hd(double *x, double *p, int *m, int *km,
      int n, int d, int K, int maxit, double eps, 
      double *llik, int progress, int *conv){
  int i, j, it;
  double del, llik_nu, *pb, *sb, *pt, fp=0.0, ttl, sump;
  ttl=(double) maxit;
  conv[0] = 0;
  pb = R_Calloc(K, double);
  sb = R_Calloc(K, double);
  pt = R_Calloc(K, double);
  //llik[0] = -n*log(n);
  del = 10.0;
  it = 0;
  while(del>eps && it<maxit){
    llik_nu=0.0;
    for(i=0; i<K; i++) pt[i]=0.0;
    for(j=0; j<n; j++){
      fp = 0.0;
      dBeta_copula_one_obs(x, m, j, n, d, km, pb, sb);
      for(i=0; i<K; i++) {
        pb[i] = p[i]*pb[i];
        fp += pb[i];
      }
      for(i=0; i<K; i++) pt[i] += pb[i]/fp;
      //llik_nu += log(fp);
    }
    for(i=0; i<K; i++){
      p[i] = pt[i]/(double)n;
    }
    
    // uniform marginal constraints
    for(i=0; i<K; i++) pt[i]=0.0;
    for(j=0; j<n; j++){
      fp = 0.0;
      dBeta_copula_one_obs(x, m, j, n, d, km, pb, sb);
      for(i=0; i<K; i++) {
        pb[i] = p[i]*pb[i];
        fp += pb[i];
      }
      for(i=0; i<K; i++) pt[i] += p[i]*sb[i]/fp;
      //llik_nu += log(fp);
    }
    sump=0.0;
    for(i=0; i<K; i++){
      p[i] = pt[i]/(double) (d*n);
      sump+=p[i];
    }
    for(i=0; i<K; i++) p[i]/=sump;    
    // end uniform marginal constraints    
    for(j=0; j<n; j++){
      fp = 0.0;
      dBeta_copula_one_obs(x, m, j, n, d, km, pb, sb);
      for(i=0; i<K; i++) {
        pb[i] = p[i]*pb[i];
        fp += pb[i];
      }
      //for(i=0; i<K; i++) pt[i] += pb[i]/fp;
      llik_nu += log(fp);
    }

    del = fabs(llik[0]-llik_nu);
    it += 1;
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
  R_Free(pb);
  R_Free(sb);
  R_Free(pt);
}
/* end function em_copula_mixbeta */

/*/////////////////////////////////////////////////////*/
/* EM Method for mixture of low dimimensional          */
/* beta(x1, i1+1, m1+1-i1),...,beta(xd, id+1, md+1-id),*/ 
/*               0<=ik<=mk, 1<=k<=d                    */
/*   with uniform marginal constraints                 */
/*/////////////////////////////////////////////////////*/
void em_copula_ld(double *p, double *pBta, double *sBta, int *m, 
      int *km, int n, int d, int maxit, double eps, 
      double *llik, int progress, int *conv){
  int i, j, it,  K=km[d];
  //int i, j, it, jj, jt, k, l, r, itmp, K=km[d];
  double del, llik_nu, *fp;//*pBeta,  *pnu;
  double tmp, ttl, sump;
  ttl=(double) maxit;
  conv[0] = 0;
  fp = R_Calloc(n, double);
  llik[0] = loglik_bern_multivar(p, K, n, pBta);
  del = 10.0;
  it = 1;
  while(del>eps && it<maxit){
    for(j=0; j<n; j++){
      fp[j] = 0.0;
      for(i=0; i<K; i++){
        fp[j] += p[i]*pBta[i+K*j];
      }
    }
    del=0.0;
    // em iteration
    sump=0.0;
    for(i=0; i<K; i++){
      tmp = 0.0; 
      for(j=0; j<n; j++) tmp += pBta[i+K*j]/fp[j];
      del+=fabs(p[i]-p[i]*tmp/(double)n);
      p[i] = p[i]*tmp/(double)n;
      sump+=p[i];
    }
    for(i=0; i<K; i++) p[i]/=sump;    
    // End em iteration

    // uniform marginal constraints
    for(j=0; j<n; j++){
      fp[j] = 0.0;
      for(i=0; i<K; i++){
        fp[j] += p[i]*pBta[i+K*j];
      }
    }
    sump=0.0;
    for(i=0; i<K; i++){
      tmp = 0.0; 
      for(j=0; j<n; j++) tmp += sBta[i+K*j]/fp[j];
      del+=fabs(p[i]-p[i]*tmp/(double)(d*n));
      p[i] = p[i]*tmp/(double)(d*n);
      sump+=p[i];
    }
    //Rprintf("loop: sump = %g\n", sump);
    for(i=0; i<K; i++) p[i]/=sump;    
    // end uniform marginal constraints

    llik_nu = loglik_bern_multivar(p, K, n, pBta);
    del += fabs(llik[0]-llik_nu);
    it += 1;
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
 // sump=0.0;
  //  for(j=0;j<K;j++) sump+=p[j]; 
   //Rprintf("\n C-ld: sum(p)=%g\n", sump);
  //R_Free(pBeta);
  R_Free(fp);
  //Free(pnu);
}
/* end function em_copula_beta_mix */
void em_copula(double *x, double *p, int *m, int *km, int K, int n, int d, int maxit, double eps, 
      double *llik, int progress, int *conv, int hd){
    if(hd==0){
      double *pBta, *sBta;
      pBta = R_Calloc(n*K, double);// very BIG memory
      sBta = R_Calloc(n*K, double);// very BIG memory
      dBeta_copula(x, m, n, d, km, pBta,sBta);
      em_copula_ld(p, pBta, sBta, m, km, n, d, maxit, eps, llik, progress, conv);
      R_Free(pBta);
      R_Free(sBta);
    }
    else{
      em_copula_hd(x, p, m, km, n, d, K, maxit, eps, llik, progress, conv);
    }  
}
/*############################################################*/
// MABLE of copula density with uniform marginal constraints
//  and model degrees either given or preselected based on marginal data 
//  or to be selected by the method of change-point
/*############################################################*/
void mable_copula(int *M, int *n, int *d, int *search, double *phat, int *mhat, 
      double *x, int *maxit,  double *eps, double *level, double *pval, double *lk, 
      double *lr, int *chpts, int *progress, int *conv,  int *hd){
  int i, j, k, prgrs, *km, *m, K, maxM, *cp, max_deg;
  int cp0=0, *Mhat, nphat=0, *a, *b, Khat; 
  double *p, *llik, pct=0.0, ttl, *res, *Phat, pv=1.0;

  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);

  a = R_Calloc(*d, int);
  b = R_Calloc(*d-1, int);
  km = R_Calloc(*d+1, int);
   m = R_Calloc(*d, int);
  Mhat = R_Calloc(*d, int);
  llik = R_Calloc(1, double);
  km[0] = 1;
  maxM=0;
  // Rprintf("km[%d]=%d\n",0,km[0]);
  for(i=1; i<=*d; i++){
    km[i]=km[i-1]*(M[i-1]+1);
    maxM = imax2(maxM, M[i-1]);
    // Rprintf("km[%d]=%d\n",i,km[i]);
  }
  K=km[*d]; 
  Khat=K;
  //K0=K;
  //Rprintf("progrss=%d\n",*progress);
  if(*progress!=0)
    Rprintf("\n Mable fit of multivariate data. This may take several minutes.\n\n");
  if(*search==0){
    if(*progress==1) prgrs=1;
    else prgrs=0;
    //for(i=0;i<K;i++) phat[i]=(phat[i]+tini/(double) K)/(1.0+tini);
    for(i=0;i<K;i++) phat[i]=1.0/(double) K;
    em_copula(x, phat, M, km, K, *n, *d, *maxit, *eps, lk, prgrs, conv, *hd);
    for(i=0; i<*d; i++) mhat[i]=M[i];
    
    //for(j=0;j<K;j++) sump+=phat[j]; 
    //Rprintf("\n C: sum(p)=%g\n", sump);
  }    
  else{
    k = 2; 
    //pval[0]=1.0;
    prgrs=0;
    max_deg=0;
    for(i=0; i<*d; i++){
      //Rprintf("M[%d]=%d,\n",i,M[i]);
      max_deg+=M[i];
    }
    //Rprintf("max_deg=%d, K=%d\n",max_deg, K);
    //Rprintf("level=%g\n",*level);
    ttl= (double) max_deg;
    p = R_Calloc(K, double); // ! BIG memory
    Phat = R_Calloc(K, double); // ! BIG memory
    //lk[0]=-1000000.0;
    //llik[0]=-10000000.0;
    while(k<max_deg && pv>*level){        
      for(i=2; i<=*d; i++){
        for(j=0;j<i;j++) a[j]=j+1;
        while(a[0]<=*d-i+1){
            //for(j=0; j<i; j++)  Rprintf("\n a[%d]=%d, ", j, a[j]);
            //Rprintf("\n");
            for(j=0;j<i-1;j++) b[j]=j+1;
            while(b[0]<=k-i+1){
                //for(j=0; j<i-1; j++)  Rprintf("\n b[%d]=%d, ", j, b[j]);
                //Rprintf("\n");
                for(j=0;j<*d;j++) m[j]=0;
                m[a[0]-1]=b[0];
                for(j=1;j<i-1;j++) m[a[j]-1]=b[j]-b[j-1];
                m[a[i-1]-1]=k-b[i-2];
                //for(j=0; j<*d; j++)  Rprintf("\n m[%d]=%d, ", j, m[j]);
                //Rprintf("\n");
                for(j=0; j<*d; j++) km[j+1]=km[j]*(m[j]+1);             
                K=km[*d];
                for(j=0;j<K;j++) p[j]=1.0/(double) K;
                em_copula(x, p, m, km, K, *n, *d, *maxit, *eps, llik, prgrs, conv, *hd);
                //Rprintf("\n k=%d,   llik = %g,  lk[%d] = %g\n", k, llik[0], k-2, lk[k-2]);
                if(llik[0]>lk[k-2]){
                    lk[k-2]=llik[0];
                    Khat=K;
                    for(j=0; j<*d; j++) Mhat[j]=m[j];
                    for(j=0; j<K; j++) Phat[j]=p[j];
                }
                if(b[0]<k-i+1) next_combo(b, i-1, k-1);
                else break;
            }
            if(a[0]<*d-i+1) next_combo(a, i, *d);
            else break;
        }
      }
      if(k>=5){
        cp[0]=k-2;
        chpt_exp(lk, lr, res, cp);
        pval[k-2]=res[0];
        chpts[k-2]=cp[0];
      }
      else{            
        pval[k-2]=1.0;
        chpts[k-2]=0;
      } 
      if(chpts[k-2]!=cp0){
        cp0=chpts[k-2];
        for(i=0; i<*d; i++) mhat[i]=Mhat[i];
        nphat=Khat;
        for(j=0;j<nphat;j++) phat[j]=Phat[j];             
      }
      pv=pval[k-2];
      Rprintf("\n chpts[%d]=%d,  lk[%d]=%g, pv = %g\n",k-2, chpts[k-2], k-2, lk[k-2], pv);
      pct += 1.0;
      if(*progress!=0) ProgressBar(fmin2(1.0,pct/ttl),"");
      R_CheckUserInterrupt();
      k+=1;
    }
    pct = ttl;
    if(*progress!=0) ProgressBar(fmin2(1.0,pct/ttl),"");
    //for(j=0;j<nphat;j++) sump+=phat[j]; 
    //Rprintf("\n C: sum(p)=%g\n", sump);
    R_Free(p); R_Free(Phat); 
    for(i=0; i<*d; i++) M[i]=mhat[i];
    d[0]=k-1;      
  }      
  R_Free(a); R_Free(b); R_Free(Mhat); R_Free(m); R_Free(km);
  R_Free(llik); R_Free(cp); R_Free(res);
}
// end of mable-copula.c
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
    SEXP args, res, tmp;
    int i, j, k, m, d, ii, jj, pc=0;
    double te, *alpha;
    MableDRStruct MDS = (MableDRStruct) ex;
    m = MDS->m; i = MDS->i; j = MDS->j; k = MDS->k; 
    d = MDS->d; alpha = MDS->alpha;
    PROTECT(args = allocVector(REALSXP, n));
    for(ii = 0; ii < n; ii++) REAL(args)[ii] = x[ii];

    PROTECT(tmp = lang2(MDS->f, args));
    PROTECT(res = eval(tmp, MDS->env));

    if(length(res) != n*(d+1))
	error("evaluation of regression function(s) gave a result of wrong length");
    if(TYPEOF(res) == INTSXP) {
	     PROTECT(res = coerceVector(res, REALSXP));
       pc = 1;
	     //res = coerceVector(res, REALSXP);
       //UNPROTECT(1); /* uprotect the old res */
       //PROTECT(res);
    } else if(TYPEOF(res) != REALSXP)
	error("evaluation of regression function(s) gave a result of wrong type");

    for(ii=0;ii<n;ii++) { 
        te  = 0.0;
        for(jj=0;jj<=d;jj++) te+=REAL(res)[ii+n*jj]*alpha[jj];
        //x[ii] = REAL(res)[ii+n*j]*REAL(res)[ii+n*k]*dbeta(x[ii],i+1,m-i+1,FALSE)*exp(te);
        x[ii] = REAL(res)[ii+n*j]*REAL(res)[ii+n*k]*(m+1)*dbinom_raw(i,m,x[ii],1-x[ii], FALSE)*exp(te);
	    if(!R_FINITE(x[ii]))
	       error("non-finite r(x) value");
    }
    UNPROTECT(3+pc);
    //UNPROTECT(3);
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
    beta_alpha = R_Calloc(n*(m+1), double);
    y = R_Calloc(n, double);

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
    R_Free(beta_alpha); R_Free(y); 
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
    Pm = R_Calloc(N, double); 
    Pm_alfa = R_Calloc(N, double);
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
    R_Free(Pm); R_Free(Pm_alfa);
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
    dwt = R_Calloc(md, double);       
    ddwt =  R_Calloc(mdd, double);
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
    R_Free(dwt); R_Free(ddwt);   
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
    Pm_alfa = R_Calloc(Nm, double);
    dPm_alfa = R_Calloc(Nmd, double);
    ddPm_alfa = R_Calloc(Nmd*d1, double);
    dwt = R_Calloc(md, double);
    ddwt =R_Calloc(md*d1, double);
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
    R_Free(Pm_alfa); R_Free(dPm_alfa); R_Free(ddPm_alfa); R_Free(dwt); R_Free(ddwt);
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
    double *beta_x, *beta_y, *fm, *Tk, *H, *Jac, llik_new;//, sump;  
    double del_em, del_nt, lam, tmp;   
    n = nx+ny;
    beta_x = R_Calloc(nx*(m+1), double);
    beta_y = R_Calloc(ny*(m+1), double);
    H = R_Calloc(d+1, double);
    Tk = R_Calloc(m+1, double);
    fm = R_Calloc(n, double);
    Jac = R_Calloc((d+1)*(d+1), double);
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
    //sump=.0;
    //for(j=0; j<=m; j++) sump+=p[j]*wt[j];
        llik_new = loglik_bern(alpha, p, ry, beta_x, beta_y, m, nx, ny, d);
        del_em = fabs(llik_new-llik[0]);
        llik[0] = llik_new;
        it_em++;
    }
    R_Free(beta_x); R_Free(beta_y); R_Free(fm); R_Free(H); R_Free(Jac); R_Free(Tk); 
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
    double del_em, del_nt, *tmp;//, sump; 
    n = nx+ny;
    beta_x = R_Calloc(nx*(m+1), double);
    beta_y = R_Calloc(ny*(m+1), double);
    H = R_Calloc(d+1, double);
    Tk = R_Calloc(m+1, double);
    fm = R_Calloc(n, double);
    Jac = R_Calloc((d+1)*(d+1), double);
    tmp = R_Calloc(d+1, double);
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
    //sump=.0;
    //for(j=0; j<=m; j++) sump+=p[j];
    //sump=.0;
    //for(j=0; j<=m; j++) sump+=p[j]*wt[j];
    R_Free(beta_x); R_Free(beta_y); R_Free(fm); R_Free(Jac); R_Free(Tk); R_Free(H); 
    R_Free(tmp);  
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
  alpha_hat = R_Calloc((k+1)*(d+1), double);
  alpha = R_Calloc(d+1, double);
  lp=M[0]*(k+1)+(k+1)*(k+2)/2;
  phat = R_Calloc(lp, double);
  pval = R_Calloc(k+1, double);
  chpts = R_Calloc(k+1, int);
  p= R_Calloc(M[1]+1, double);
  wt = R_Calloc(M[1]+1, double);
  wt_ahat = R_Calloc(lp, double);
  llik = R_Calloc(1, double);
  lk = R_Calloc(k+1, double);
  lr = R_Calloc(k+1, double);
  z = R_Calloc(n0+n1, double);
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);

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
  R_Free(alpha_hat); R_Free(alpha); R_Free(phat); R_Free(pval);    
  R_Free(chpts); R_Free(p); R_Free(wt); R_Free(wt_ahat);  
  R_Free(llik); R_Free(lk); R_Free(lr); R_Free(z); R_Free(cp); R_Free(res);

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
  phat = R_Calloc(lp, double);
  pval = R_Calloc(k+1, double);
  chpts = R_Calloc(k+1, int);
  tmp = 1.0+(double)(k*(k-1));
  wt = R_Calloc(M[1]+1, double);
  wt_ahat = R_Calloc(lp, double);
  p= R_Calloc(M[1]+1, double);
  llik = R_Calloc(1, double);
  lk = R_Calloc(k+1, double);
  lr = R_Calloc(k+1, double);
  z = R_Calloc(n0+n1, double);
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);

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
  
  R_Free(lk); R_Free(lr); R_Free(wt); R_Free(wt_ahat); R_Free(p); R_Free(phat); 
  R_Free(chpts); R_Free(pval); R_Free(res); R_Free(llik); R_Free(cp); R_Free(z);
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
  double del_em, del_nt, *tmp, *h, *Sig, *Er, *dFm_s, *ddFm_s; 
  double *Tk, *Pm_s, *dPm_s, *ddPm_s, *Pi_theta, *T_zero;
  T_zero = R_Calloc(m1,double); Tk = R_Calloc(m1,double);
  Pi_theta = R_Calloc(Nm,double); tmp = R_Calloc(d1,double);
  Pm = R_Calloc(Nm,double); Pm_s = R_Calloc(Nm,double);
  dPm_s = R_Calloc(Nmd,double); ddPm_s = R_Calloc(Nmd*d1,double);
  Fm = R_Calloc(N,double); Fm_s = R_Calloc(N,double);
  H = R_Calloc(d1,double); Jac = R_Calloc(d1*d1,double);

  h = R_Calloc(d1*d1,double); Sig = R_Calloc(d1*d1,double);  
  dFm_s = R_Calloc(N*d1,double); ddFm_s = R_Calloc(N*d1*d1,double);
  Er = R_Calloc(d1,double);
  
  cpBeta(t, m, N, Pm); 
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
  
  
  R_Free(H); R_Free(Jac); R_Free(Tk); 
  R_Free(Er); R_Free(h); R_Free(Sig); R_Free(dFm_s); R_Free(ddFm_s);
  R_Free(tmp); R_Free(T_zero); R_Free(Pi_theta); R_Free(Pm); R_Free(Pm_s);
  R_Free(dPm_s); R_Free(ddPm_s); R_Free(Fm); R_Free(Fm_s);
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
  Tk = R_Calloc(m1,double);
  Pm = R_Calloc(Nm,double); Pm_s = R_Calloc(Nm,double);
  Fm = R_Calloc(N,double); Fm_s = R_Calloc(N,double);

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
  R_Free(Tk); R_Free(Pm); R_Free(Pm_s); R_Free(Fm); R_Free(Fm_s);
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
  phat = R_Calloc(lp, double);
  pval = R_Calloc(k+1, double);
  chpts = R_Calloc(k+1, int);
  tmp = 1.0+(double)(k*(k+1));
  //Rprintf("  dim:  %d\n", (k+1)*(d+1));
  lk = R_Calloc(k+1, double);
  lr = R_Calloc(k+1, double);
  p = R_Calloc(M[1]+1, double);
  wt = R_Calloc(M[1]+1, double);
  llik = R_Calloc(1, double);
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);

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

  R_Free(llik); R_Free(wt); R_Free(phat); R_Free(p); R_Free(lk); R_Free(lr);
  R_Free(cp); R_Free(res); R_Free(cp); R_Free(res);
  setAttrib(ans, R_NamesSymbol, ansnames);
  UNPROTECT(2);
  return ans;
}

/*////////////////////////////////////////////////////*/
/*  Model degree m selection by change-point method   */
/*    with (alpha-hat, phat) based on grouped data    */
/*////////////////////////////////////////////////////*/
SEXP C_mable_dr_group(SEXP args){
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
  phat = R_Calloc(lp, double);
  pval = R_Calloc(k+1, double);
  chpts = R_Calloc(k+1, int);
  tmp = 1.0+(double)(k*(k-1));
  //Rprintf("  dim:  %d\n", (k+1)*(d+1));
  lk = R_Calloc(k+1, double);
  lr = R_Calloc(k+1, double);
  alpha_hat = R_Calloc((k+1)*(d+1), double);
  alpha = R_Calloc(d+1, double);
  se = R_Calloc(d+1, double);
  se_hat = R_Calloc((k+1)*(d+1), double);
  p= R_Calloc(M[1]+1, double);
  wt = R_Calloc(M[1]+1, double);
  llik = R_Calloc(1, double);
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
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

  R_Free(llik); R_Free(wt); R_Free(p); R_Free(alpha_hat); R_Free(alpha); R_Free(phat); 
  R_Free(lk); R_Free(lr); R_Free(chpts); R_Free(pval); R_Free(cp); R_Free(res); 
  R_Free(se); R_Free(se_hat);
  setAttrib(ans, R_NamesSymbol, ansnames);
  UNPROTECT(2);
  return ans;
}


/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*  Maximum Approximate Bernstein likelihood Estimation   */
/*  in Accelerated Failure Time Regression model based    */
/*               on Interval Censored data                */
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
          double *BSz, double *BSz2, double *tau){
  int i,j, n=n0+n1;
  double llkhd, fz, dSz;
  llkhd = 0.0;
  for(i=0; i<n0; i++){
    fz = 0.0; 
    for(j=0; j<=m; j++){
      fz += p[j]*BSz2[i+n*j];
    }
    llkhd += gx[i]+log(fz);
  }
  for(i=n0; i<n; i++){
    dSz=0.0;
    for(j=0; j <= m; j++){
      dSz += p[j]*(BSz[i+n*j]-BSz2[i+n*j]); 
    }
    llkhd += log(dSz);
  }
  llkhd -= n0*log(tau[0]);
  return llkhd;
}
/*/////////////////////////////////////////////////////////////////*/
/* Derivatives of loglikelihood ell(gamma, p) wrt gamma, aft model */
/*/////////////////////////////////////////////////////////////////*/
// 
void logblik_aft_derv(double *gama, double *p, int d, int m, 
    double *x, double *x0, double *tau, double *gx, double *z, double *z2, 
    int n0, int n1, double *ell, double *dell, double *ddell){
  int i,j,k, n=n0+n1, nmp2=n*(m+2);
  double tmp1=0.0, tmp2=0.0, A, B, C;
  double *BSz, *BSz2, *bz, *bz2;   
  bz = R_Calloc(nmp2, double);  
  bz2 = R_Calloc(nmp2, double); 
  BSz = R_Calloc(nmp2, double);  
  BSz2 = R_Calloc(nmp2, double); 
  ell[0]=0.0;
  for(i=0; i<d; i++){ 
    dell[i]=0.0; 
    for(j=0; j<d; j++) ddell[i+d*j]=0.0;}
  for(k=0; k<n0; k++) ell[0] += gx[k];
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
      dell[i] += (1.0+tmp1)*(x[k+n*i]-x0[i]); 
      for(j=0; j<d; j++)
        ddell[i+d*j] -= (tmp1*tmp1-tmp2)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
    }
  }
  ell[0] -= n0*log(tau[0]); //????
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
      dell[i]-=A*(x[k+n*i]-x0[i])/tmp1;           
      for(j=0;j<d;j++){
        ddell[i+d*j]-=(A/tmp1)*(A/tmp1+1.0)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
        ddell[i+d*j]-=(tmp2/tmp1)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
      }
    }
  }
  //for(i=0; i<d; i++)
     //for(j=0;j<d;j++) Rprintf("\n ddell[%d,%d]=%f", i,j, ddell[i+d*j]);

  R_Free(bz); R_Free(bz2);
  R_Free(BSz); R_Free(BSz2);
}

/* Extended version of struct opt_struct */
typedef struct mable_aft_struct
{
    int m, n0, n1, known_tau;
    double *x0, *x, *p, *y, *y2;
    double *tau;    
} mable_aft_struct, *MableAFTStruct;

// minus log-likelihood of AFT model
static double deviance_aft(int npar, double *par, void *ex)
{
  int i, j, m, n0, n1, n;
  double val=0.0, vz, *z, *z2, *gx, *BSz, *BSz2;
  MableAFTStruct AS = (MableAFTStruct) ex;
  m = AS->m; n0 = AS->n0; n1 = AS->n1; n = n0+n1;
  
  z = R_Calloc(n, double);  
  z2 = R_Calloc(n, double);  
  gx = R_Calloc(n, double);  
  BSz = R_Calloc(n*(m+2), double);  
  BSz2 = R_Calloc(n*(m+2), double);  
  egxmx0(par, npar, AS->x, n, gx, AS->x0);
  if(AS->known_tau!=1){
    AS->tau[0]=AS->tau[1];
    for(i=0;i<n;i++){
      z[i] = AS->y[i]*gx[i];
      z2[i] = AS->y2[i]*gx[i];
      AS->tau[0]=fmax(AS->tau[0], z[i]);
      if(AS->y2[i]<=AS->tau[1]) AS->tau[0]=fmax(AS->tau[0], z2[i]);
    }
    AS->tau[0]+=1.0/(double) n;
    for(i=0;i<n;i++){
      z[i] = z[i]/AS->tau[0];
      z2[i] = z2[i]/AS->tau[0];
      gx[i] = log(gx[i]);
    }
  }
  else{
    for(i=0;i<n;i++){
      z[i] = AS->y[i]*gx[i];
      z2[i] = AS->y2[i]*gx[i];
      gx[i] = log(gx[i]);
    }
  }
  Bdata(z, m, 0, n, BSz);
  Bdata(z2, m, n0, n1, BSz2);


  for(i=0; i<n0; i++){
    vz = 0.0; 
    for(j=0; j<=m; j++){
      vz += AS->p[j]*BSz2[i+n*j];
    }
    val -= gx[i]+log(vz);
  }
  //Rprintf("\n val=%f,  \n",val);
  for(i=n0; i<n; i++){
    vz=0.0;
    for(j=0; j <= m; j++){
      vz += AS->p[j]*(BSz[i+n*j]-BSz2[i+n*j]); 
    }
    val -= log(vz);
  }
  //Rprintf("\n val=%f,  \n",val);
  val += n0*log(AS->tau[0]);
  //Rprintf("\n val=%f,  \n",val);
  R_Free(z); R_Free(z2); R_Free(gx); R_Free(BSz); R_Free(BSz2);
  return val;
}


// derivative of minus log-likelihood of AFT model
static void D_deviance_aft(int npar, double *par, double *df, void *ex)
{
  int i, j, k, m, n0, n1, n, nmp2;
  double tmp=0.0, A, B;//, C;
  double *z, *z2, *gx, *BSz, *BSz2, *bz, *bz2;   
  MableAFTStruct AS = (MableAFTStruct) ex;
  m = AS->m; n0 = AS->n0; n1 = AS->n1; n = n0+n1;
  nmp2=n*(m+2);
  z = R_Calloc(n, double);  
  z2 = R_Calloc(n, double);  
  gx = R_Calloc(n, double);  
  bz = R_Calloc(nmp2, double);  
  bz2 = R_Calloc(nmp2, double); 
  BSz = R_Calloc(nmp2, double);  
  BSz2 = R_Calloc(nmp2, double);  
  egxmx0(par, npar, AS->x, n, gx, AS->x0);
  if(AS->known_tau!=1){
    AS->tau[0]=AS->tau[1];
    for(i=0;i<n;i++){
      z[i] = AS->y[i]*gx[i];
      z2[i] = AS->y2[i]*gx[i];
      AS->tau[0]=fmax(AS->tau[0], z[i]);
      if(AS->y2[i]<=AS->tau[1]) AS->tau[0]=fmax(AS->tau[0], z2[i]);
    }
    AS->tau[0]+=1.0/(double) n;
    for(i=0;i<n;i++){
      z[i] = z[i]/AS->tau[0];
      z2[i] = z2[i]/AS->tau[0];
      gx[i] = log(gx[i]);
    }
  }
  else{
    for(i=0;i<n;i++){
      z[i] = AS->y[i]*gx[i];
      z2[i] = AS->y2[i]*gx[i];
      gx[i] = log(gx[i]);
    }
  }

  Bdata(z, m, 0, n, BSz);//1-B(z),1-B(z2) 
  Bdata(z2, m, 0, n, BSz2);
  Bdata(z, m, n, 0, bz); // beta(z), beta(z2)
  Bdata(z2, m, n, 0, bz2);
  for(i=0; i<npar; i++) df[i]=0.0; 

  for(k=0; k<n0; k++){
    A = 0.0;
    B = 0.0;
    //C = 0.0;
    for(j=0;j<m;j++){
      A += AS->p[j]*bz[k+n*j];
      B += AS->p[j]*(j*bz[k+n*j]-(j+1)*bz[k+n*(j+1)]);
      //C += AS->p[j]*(j*j*bz[k+n*j]-(j+1)*(2*j+1)*bz[k+n*(j+1)]+(j+1)*(j+2)*bz[k+n*(j+2)]);
    }  
    A += AS->p[m]*bz[k+n*m]; 
    B += AS->p[m]*m*bz[k+n*m];
    //C += AS->p[m]*m*m*bz[k+n*m];       
    tmp = B/A;
    for(i=0; i<npar; i++){
      df[i] -= (1.0+tmp)*(AS->x[k+n*i]-AS->x0[i]); 
    }
  }
  for(k=n0; k<n; k++){
    A = 0.0;
    B = 0.0;
    //C = 0.0;
    tmp = 0.0;
    for(j=0;j<=m;j++){
      A += AS->p[j]*(z[k]*bz[k+n*j]-z2[k]*bz2[k+n*j]);
      tmp += AS->p[j]*(BSz[k+n*j]-BSz2[k+n*j]);
      //B += AS->p[j]*(j*bz[k+n*j]-(j+1)*bz[k+n*(j+1)]);
      //C += AS->p[j]*(j*bz2[k+n*j]-(j+1)*bz2[k+n*(j+1)]);
    }
    for(i=0; i<npar; i++){
      df[i]+=A*(AS->x[k+n*i]-AS->x0[i])/tmp;           
    }
  }

  R_Free(z); R_Free(z2); R_Free(gx); R_Free(BSz); R_Free(BSz2); R_Free(bz); R_Free(bz2);
}


/*////////////////////////////////////////////////////*/
/*  Find maximizer gamma of ell(gamma, p) by Newton   */
/*       method for a gvien p for AFT model           */
/*////////////////////////////////////////////////////*/
void gofp_aft_nt(double *gama, int d, double *p, int m, double *y, double *y2, 
      double *x, double *x0, double *tau, double *gx, double *z, double *z2, 
      int n0, int n1, double *ell, double *dell, double *ddell, double eps, 
      int maxit, int prog, int known_tau, int *conv){
  int i,j, it=0, n=n0+n1;
  double delta=0.0, *tmp;  
  tmp = R_Calloc(d, double);
//  Rprintf("known_tau=%d\n",known_tau);
  logblik_aft_derv(gama, p, d, m, x, x0, tau, gx, z, z2, n0, n1, ell, dell, ddell);
  for(i=0;i<d;i++) delta+=fabs(dell[i]); 
  while(it<maxit && delta>eps){
    minverse(ddell, d);  
    for(i=0;i<d;i++){
      tmp[i] = 0.0;
      for(j=0;j<d;j++) tmp[i] += ddell[i+d*j]*dell[j];
    }
    delta = 0.0;
    for(i=0;i<d;i++){
      gama[i] -= tmp[i];
      delta += fabs(tmp[i]);
    }
    egxmx0(gama, d, x, n, gx, x0);
    if(known_tau!=1){
      tau[0]=tau[1];
      for(i=0;i<n;i++){
        z[i] = y[i]*gx[i];
        z2[i] = y2[i]*gx[i];
        tau[0]=fmax(tau[0], z[i]);
        if(y2[i]<=tau[1]) tau[0]=fmax(tau[0], z2[i]);
      }
      tau[0]+=1.0/(double) n;
      for(i=0;i<n;i++){
        z[i] = z[i]/tau[0];
        z2[i] = z2[i]/tau[0];
        gx[i] = log(gx[i]);
      }
    }
    else{
      for(i=0;i<n;i++){
        z[i] = y[i]*gx[i];
        z2[i] = y2[i]*gx[i];
        gx[i] = log(gx[i]);
      }
    }
    logblik_aft_derv(gama, p, d, m, x, x0, tau, gx, z, z2, n0, n1, ell, dell, ddell);
    for(i=0;i<d;i++) delta+=fabs(dell[i]);
    it++;
    R_CheckUserInterrupt();
  }
  *conv = (it<maxit) ? 0:1;
  if(prog==0) Rprintf("NT: m=%d, it=%d, del=%e, llik=%f\n",m,  it, delta, ell[0]);
  R_Free(tmp); 
}


/*////////////////////////////////////////////////////*/
/*   maximizer p of ell(gamma, p) for a gvien gamma   */
/*  gx:  gama*x.tilde, where gama is the given        */
/*       regression coefficient of the AFT model      */
/*  ss: step-size epsilon:  p = (1-ss)p+ss Psi(p),    */
/*        default ss=1 so that p = Psi(p)             */
/*////////////////////////////////////////////////////*/
void pofg_aft(double *p, int m, double *gx, int n0, int n1, 
    double *BSz, double *BSz2, double *tau, double *llik, 
    double eps, int maxit, int prog, int *conv, double *delta){
  int i, j, n=n0+n1, mp1=m+1, it=0;
  double  dlt=1.0, dSz;
  double *Tmp, *pnu, lik_nu; 
  Tmp=R_Calloc(mp1, double); 
  pnu=R_Calloc(mp1, double);
  conv[0]=0;
  while(dlt>eps && it<maxit){
    for(j=0;j<=m;j++) pnu[j]=0.0;
    // p = p *Psi(p) 
    lik_nu=0.0;
    for(i=0; i<n0;i++){
      dSz=0.0; 
      for(j=0;j<=m;j++){
        Tmp[j]=BSz2[i+n*j]*p[j];
        dSz+=Tmp[j];
      }
      for(j=0;j<=m;j++){
        pnu[j]+=Tmp[j]/dSz;
      }
      lik_nu += gx[i]+log(dSz);
    }
    for(i=n0; i<n;i++){
      dSz=0.0;  
      for(j=0;j<=m;j++){
        Tmp[j]=(BSz[i+n*j]-BSz2[i+n*j])*p[j];
        dSz+=Tmp[j];
      }
      for(j=0;j<=m;j++){
        pnu[j]+=Tmp[j]/dSz;
      }
      lik_nu += log(dSz);
    }  
    for(j=0;j<=m;j++) pnu[j] /= (double) n;
    //pnu<-(1-ss)*p+ss*pnu
    lik_nu -= n0*log(tau[0]);
    if(it>0) dlt=fabs(llik[0]-lik_nu);
    llik[0]=lik_nu;
    for(j=0;j<=m;j++) p[j]=pnu[j];
    it++;  
    R_CheckUserInterrupt();
  }
  if(prog==0) Rprintf("EM: m=%d, it=%d, del=%e, llik=%f\n",m,  it, dlt, llik[0]);
  conv[0]=(it<maxit)? 0:1;
  delta[0]=dlt;
   R_Free(pnu); R_Free(Tmp); 
}
/*////////////////////////////////////////////////////*/
/*  Maximum approx. Bernstein likelihood estimate of  */
/*   (gamma, p) with a fixed degree m for AFT model   */
/*////////////////////////////////////////////////////*/
void mable_aft_m(double *gama, double *p, int *dm, double *x, double *y, double *y2, 
       double *tau, int *N, double *x0, double *ell, double *ddell, double *EPS, 
       int *MAXIT, int *progress, int *conv, double *delta, int *known_tau, int *method){
  int i, n0=N[0], n1=N[1], n=n0+n1, d=dm[0], m=dm[1],  nmp2 =n*(m+2), it=0;
  int maxit=MAXIT[0], maxit_em=MAXIT[1], prog=1, nbt1=0; 
  double eps=EPS[0], eps_em=EPS[1], pct=0.0;
  double *z, *z2, *gx, *BSz, *BSz2, dlt=1.0;//*xt, tini, 
  double *new_ell, *dell, val = 0.0, abstol=-INFINITY, reltol=eps;
  int *mask, trace=0, fncount = 0, grcount = 0, nREPORT=10;
  int ifail = 0, conv_nt=0;
  MableAFTStruct AS;
  AS = (MableAFTStruct) R_alloc(1, sizeof(mable_aft_struct));
  //tini=0.00001;// tini is used to make sure p is in interior of S(m)
  new_ell = R_Calloc(1, double);  
  dell = R_Calloc(d, double);  
  z = R_Calloc(n, double);  
  z2 = R_Calloc(n, double);  
  gx = R_Calloc(n, double);  
  BSz = R_Calloc(nmp2, double);  
  BSz2 = R_Calloc(nmp2, double);  
  mask = R_Calloc(d, int);
  for (i = 0; i < d; i++) mask[i] = 1;
  egxmx0(gama, d, x, n, gx, x0);
  if(*known_tau!=1){
    tau[0]=tau[1];
    for(i=0;i<n;i++){
      z[i] = y[i]*gx[i];
      z2[i] = y2[i]*gx[i];
      tau[0]=fmax(tau[0], z[i]);
      if(y2[i]<=tau[1]) tau[0]=fmax(tau[0], z2[i]);   
    }
    tau[0]+=1.0/(double) n;
    for(i=0;i<n;i++){
      z[i] = z[i]/tau[0];
      z2[i] = z2[i]/tau[0];
      gx[i] = log(gx[i]);
    }
  }
  else{
    for(i=0;i<n;i++){
      z[i] = y[i]*gx[i];
      z2[i] = y2[i]*gx[i];
      gx[i] = log(gx[i]); //???
      if(y2[i]<=1 && z2[i]>1) nbt1++;
    }
    if(nbt1==n){
      Rprintf("\n");
      warning("May need to try another baseline 'x0' and/or a larger truncation time 'tau'.\n");}
  }
  Bdata(z, m, 0, n, BSz);
  Bdata(z2, m, n0, n1, BSz2);
  if(*progress==1){
    Rprintf("\n Mable fit of AFT model with a given degree ... \n"); 
    ProgressBar(pct,"");} 
  pofg_aft(p, m, gx, n0, n1, BSz, BSz2, tau, ell, eps_em, maxit_em, prog, conv, delta);
  //Rprintf("\n ell=%f, ell0=%f\n",ell[0], ell[1]);
  //ell[1]=ell[0];
  
  if(*method==1){
    logblik_aft_derv(gama, p, d, m, x, x0, tau, gx, z, z2, n0, n1, ell, dell, ddell);
    if(matrix_singular(ddell, d)==1){
      *method=0;
      //Rprintf("\n Singular Hessian matrix. Use Quasi-Newton Method.");
    }
  }

  AS->m=m; AS->n0=n0; AS->n1=n1; AS->known_tau=*known_tau;
  AS->x0=x0; AS->x=x; AS->y=y; AS->y2=y2; AS->p=p; AS->tau=tau;
  //Rprintf("\n gamma = (%f",gama[0]);
  //for(i=1;i<d;i++) Rprintf(" ,%f", gama[i]);
  //Rprintf(").\n");
  //Rprintf("\n Initial value = %f.\n", deviance_aft(d, gama, AS));
  while(it<maxit && (dlt>eps || ell[0]<ell[1])){
    if(*method==0){
      vmmin(d, gama, &val, deviance_aft, D_deviance_aft, maxit, trace, mask, abstol, 
	      reltol, nREPORT, (void *)AS, &fncount, &grcount, &ifail);
      new_ell[0]=-val;
      conv[1]=ifail;
      egxmx0(gama, d, x, n, gx, x0);
      if(*known_tau!=1){
        tau[0]=AS->tau[1];
        for(i=0;i<n;i++){
          z[i] = y[i]*gx[i];
          z2[i] = y2[i]*gx[i];
          tau[0]=fmax(tau[0], z[i]);
          if(y2[i]<=tau[1]) tau[0]=fmax(tau[0], z2[i]);
        }
        tau[0]+=1.0/(double) n;
        for(i=0;i<n;i++){
          z[i] = z[i]/tau[0];
          z2[i] = z2[i]/tau[0];
          gx[i] = log(gx[i]);
        }
      }
      else{
        for(i=0;i<n;i++){
          z[i] = y[i]*gx[i];
          z2[i] = y2[i]*gx[i];
          gx[i] = log(gx[i]);
        }
      }
    }
    else {//if(*method==1)
      gofp_aft_nt(gama, d, p, m, y, y2, x, x0, tau, gx, z, z2, n0, n1, 
            new_ell, dell, ddell, eps, maxit, prog, known_tau[0], &conv_nt); 
      conv[1]=conv_nt;
    }
    nbt1=0;
    Bdata(z, m, 0, n, BSz);
    Bdata(z2, m, n0, n1, BSz2);
    pofg_aft(p, m, gx, n0, n1, BSz, BSz2, tau, new_ell, eps_em, maxit_em, prog, conv, delta);
    dlt = fabs(new_ell[0]-ell[0]);
    ell[0]=new_ell[0];
// 
    logblik_aft_derv(gama, p, d, m, x, x0, tau, gx, z, z2, n0, n1, ell, dell, ddell);
    if(matrix_singular(ddell, d)==1) *method=0;
    else *method=1;
    
    AS->p=p; AS->tau=tau;
    
    //Rprintf("\n ell=%f, ell0=%f\n", ell[0], new_ell[0]);
    //Rprintf("\n pmp1=%d, x0=%f,  %f\n", *pmp1, x0[0], x0[1]);
    pct=fmin(it/(double)maxit, fmax(0.0, 1-fabs(dlt-eps)));
    if(*progress==1) ProgressBar(pct,""); 
    it++;
    R_CheckUserInterrupt();
    // Rprintf("         mable-m: it=%d, del=%f, ell=%f\n", it, del, ell[0]);
  }
  if(*progress==1){
    ProgressBar(1.0,""); 
    Rprintf("\n");}
  delta[0]=dlt;
  conv[0]*=2; 
  if(it<maxit) conv[0] += 1;
    //warning("\nThe maximum iterations were reached \nwith a delta = %f.\n", dlt);
  
  //Rprintf("mable-m: it=%d, delta=%f\n", it, dlt);
  // calculate Hessian matrix?? 
  // this is unnecessary because the model will be fit again after m is selected
  if(*method==0)
    logblik_aft_derv(gama, p, d, m, x, x0, tau, gx, z, z2, n0, n1, ell, dell, ddell);
  if(matrix_singular(ddell, d)==0) 
      minverse(ddell, d); //Sig=-n*ddell
  //else Rprintf("\n Singular Hessian matrix.\n");
  R_Free(new_ell); R_Free(dell); R_Free(z); R_Free(z2); 
  R_Free(gx); R_Free(BSz); R_Free(BSz2); R_Free(mask);
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
/*      dm: (d,m)                                            */
/*       x: d-dim covariate                                  */
/*      x0: baseline covariate value, default is 0           */
/*///////////////////////////////////////////////////////////*/
void mable_aft(int *M, double *gama, int *dm, double *p, double *x, double *y,    
      double *y2, double *tau, int *N, double *x0, double *lk, double *lr, 
      double *ddell, double *EPS, int *MAXIT, int *progress, double *pval, 
      int *chpts, double *level, int *conv, int *known_tau){
  int i, j, d=dm[0], k=M[1]-M[0], tmp=0, cp0=1, m1=1; //, cp1=1
  int m, mp1, lp, prg=1-*progress, *cp, method=1, *icon; 
  double pct=0.0, ttl, lnn=-1.0e20, pv0=1.0, delta=0.0;//, pv1=1.0;   
  double *res, *phat, *ghat, *ell, *lrcp, *dlt;//, sm=0.0; 
  m=M[0]; 
  lp=(k+1)*(2*m+k+2)/2;
  icon = R_Calloc(2, int); // (conv_em,conv_nt)
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  phat=R_Calloc(lp, double);
  ghat=R_Calloc(d*(k+1), double);
  ell=R_Calloc(2, double);
  lrcp=R_Calloc(k+1, double); // k or k+1 ??
  dlt=R_Calloc(1, double);
  if(*progress==1) {Rprintf("\n Mable fit of AFT model ... \n");
      ProgressBar(0.0,""); }
  ttl=(double) (k+2)*(k+1);
  mp1=m+1;
  dm[1]=m;
  ell[1]=lnn;
  icon[0]=0;
  icon[1]=0;
  // change to use "Nelder-Mead" method???? with maxit=500
  mable_aft_m(gama, p, dm, x, y, y2, tau, N, x0, ell, ddell, EPS, MAXIT, &prg, 
        icon, dlt, known_tau, &method);
  method=1; //??
  tmp=mp1;
  for(i=0;i<tmp;i++) phat[i]=p[i];
  for(i=0;i<d;i++) ghat[i]=gama[i];
  lk[0]=ell[0];
  ell[1]=ell[0];
  pval[0]=1.0;
  chpts[0]=0;
  delta=dlt[0];
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
    //Rprintf("\n m0=%d, m1=%d, m=%d, k=%d, lp=%d, tmp=%d\n",M[0], M[1], m, k, lp, tmp+mp1);
    mable_aft_m(gama, p, dm, x, y, y2, tau, N, x0, ell, ddell, EPS, MAXIT, &prg, 
      icon, dlt, known_tau, &method);
    method=1; //??
    for(j=0;j<mp1;j++) phat[j+tmp]=p[j];
    tmp+=mp1;
    //if((i+1)*d>(k+1)*d) error("\n (i+1)*d>(k+1)*d\n");
    for(j=0;j<d; j++) ghat[j+i*d]=gama[j];
    lk[i]=ell[0]; 
    ell[1]=ell[0]; 
    //Rprintf("\n lk[%d]=%f, lp=%d, temp=%d\n",i, lk[i], lp, tmp);
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
    if(pval[i]<pv0){
      pv0=pval[i];
      cp0=chpts[i];
      delta=dlt[0];
      m1 = m;
      for(j=0;j<i; j++) lrcp[j]=lr[j];
    }
    //Rprintf("\n chpts[%d]=%d, pval[%d]=%f\n",i, chpts[i], i, pval[i]);
 
    R_CheckUserInterrupt();
    pct +=2*(i+1); 
    if(*progress==1) ProgressBar(fmin(1.0,pct/ttl),""); 
    i++;
  }
  if(*progress==1){
    ProgressBar(1.00,"");
    Rprintf("\n");}
  //if(icon[0]>0) Rprintf("\nAt least one EM iteration did not converge.\n"); 
  //if(icon[1]>0) Rprintf("\nAt least one (quasi) Newton iteration did not converge.\n"); 
  if(m==M[1]){
    Rprintf("\nMaximum degree reached.\n"); 
    if(pv0>.2) conv[0]=1; 
    else conv[0]=0;
    //if(pv0<=0.1) 
    Rprintf("A degree with smallest p-value of the change-point %f is returned.\n", pv0); 
    //else error("Search for a model degree failed with smallest p-value of the change-point %f.\n", pv0);
  }
  else M[1]=m1;
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  m=cp0+M[0];
  dm[1]=m;
  //if(tmp+m+1>lp) error("\n tmp+m+1>lp\n");
  for(j=0;j<=m;j++) p[j]=phat[tmp+j];
  for(j=0; j<d; j++) gama[j]=ghat[d*cp0+j];
  // print gamma
  //Rprintf("\n gamma = (%.3f",gama[0]);
  //for(i=1;i<d;i++) Rprintf(" ,%.3f", gama[i]);
  //Rprintf(").");
  // print p
  //Rprintf("\n p = (%.4f",p[0]);
  //for(i=1;i<=m;i++) Rprintf(" ,%.4f", p[i]);
  //Rprintf(").");
  // Add an argument method = 0 (quasi-newton) or 1 (newton) or other methods???
  // if ddell is singular method=0, or method = 1 in the last mable fit below
  //logblik_aft_derv(gama, p, d, m, x, x0, tau, gx, z, z2, n0, n1, ell, dell, ddell);
  //if(matrix_singular(ddell, d)==0) 
  icon[0]=0;
  icon[1]=0;
  method=1;
  mable_aft_m(gama, p, dm, x, y, y2, tau, N, x0, ell, ddell, EPS, MAXIT, &prg, 
        icon, dlt, known_tau, &method);
  // print gamma
  //Rprintf("\n gamma = (%.3f",gama[0]);
  //for(i=1;i<d;i++) Rprintf(" ,%.3f", gama[i]);
  //Rprintf(").");
  // print p
  //Rprintf("\n p = (%.4f",p[0]);
  //for(i=1;i<=m;i++) Rprintf(" ,%.4f", p[i]);
  //Rprintf(").");
  /* decimal convergence code :            */
  /* conv=(cp>0)+2*(em>0 || nt>0 || it=maxit in mablem) */
  if(icon[0]>0 || icon[1]>0) conv[0]+=2;
  //if(icon[1]>0) conv[0]+=4;
  for(j=0;j<M[1]-M[0]; j++) lr[j]=lrcp[j];
  *level=delta;
  dm[0]=M[1]-M[0];
  R_Free(icon); R_Free(cp); R_Free(phat); R_Free(ghat);
  R_Free(res); R_Free(ell); R_Free(lrcp); R_Free(dlt); 
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
void mable_aft_gamma(int *M, double *gama, int *dm, double *x, double *y, 
    double *y2, int *N, double *x0, double *lk, double *lr, double *p, 
    double *ddell, double *eps, int *maxit, int *progress, double *pval, 
    int *chpts, double *level, int *conv, double *delta, double *tau, int *known_tau){        
  int i,j, d=dm[0], k=M[1]-M[0], *cp, cp0=1, cp1=1, n0=N[0], n1=N[1], n=n0+n1, tmp, itmp=0;
  int m=M[1], mp1=m+1,  mp2=m+2, lp=(k+1)*M[0]+(k+1)*(k+2)/2, nbt1=0; 
  double tini=.0001, pct=0.0, ttl, *z, *z2, *res, pv0=1.0, pv1=1.0;
  double *ell, *dell, *gx, *BSz, *BSz2, *phat; 
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  phat=R_Calloc(lp, double);
  ell = R_Calloc(1, double);
  dell = R_Calloc(d, double);
  BSz = R_Calloc(n*mp2, double);  
  BSz2 = R_Calloc(n*mp2, double); 
  z = R_Calloc(n, double);
  z2 = R_Calloc(n, double);
  gx = R_Calloc(n, double);
  if(*progress==1) {Rprintf("\n Mable fit of AFT model with given regression coefficients ... \n");
      ProgressBar(0.0,""); }
  ttl = (double)((k+2)*(k+1));
  egxmx0(gama, d, x, n, gx, x0);
  if(*known_tau!=1){
    tau[0]=tau[1]; 
    for(i=0;i<n;i++){
      z[i] = y[i]*gx[i];
      z2[i] = y2[i]*gx[i];
      tau[0]=fmax(tau[0],z[i]);
      if(y2[i]<=tau[1]) tau[0]=fmax(tau[0],z2[i]);
    }
    tau[0]+=1.0/(double) n;
    for(i=0;i<n;i++) {
      z[i] = z[i]/tau[0];
      z2[i] = z2[i]/tau[0];
      gx[i] = log(gx[i]);
    }
  }
  else{
    for(i=0;i<n;i++){
      z[i] = y[i]*gx[i];
      z2[i] = y2[i]*gx[i];
      gx[i] = log(gx[i]);
      // add check to see any z, z2 are bigger than 1
      if (y2[i]<=1 && z2[i]>1)  nbt1++;
    }
    if(nbt1==n){
      Rprintf("\n");
      warning("May need to try another baseline 'x0' and/or a larger truncation time 'tau'.\n");
    }
  }
  m=M[0]; 
  mp1=m+1;
  mp2=m+2;
  Bdata(z, m, 0, n, BSz);
  Bdata(z2, m, n0, n1, BSz2);
  //for(i=0;i<=m;i++) p[i]=1.0/(double) mp1; 
  pofg_aft(p, m, gx, n0, n1, BSz, BSz2, tau, ell, *eps, *maxit, *progress, conv, delta);
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
    pofg_aft(p, m, gx, n0, n1, BSz, BSz2, tau, ell, *eps, *maxit, *progress, conv, delta);
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
      //warning("\nThe maximum candidate degree has been reached \n
      //with a p-value of the change-point %f.\n", res[0]);
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
  R_Free(phat);  R_Free(ell);  R_Free(dell);
  R_Free(BSz); R_Free(BSz2); R_Free(z); R_Free(z2); 
  R_Free(gx);  R_Free(cp); R_Free(res);
}
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*  Maximum Approximate Bernstein likelihood Estimation   */
/*     in Proportional Hazards Regression model based     */
/*               on Interval Censored data                */
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
    //Rprintf("lk1: lk=%f\n", llkhd);
    for(i=n0; i<n; i++){
        Sy=0.0;
        Sy2=0.0;
        for(j=0; j<=mp1; j++){
            Sy += p[j]*BSy[i+n*j];
            Sy2 += p[j]*BSy2[i+n*j]; 
        }
        //Rprintf("Sy: Sy=%f\n", Sy);
        llkhd += log(R_pow(Sy, egx[i])-R_pow(Sy2, egx[i]));
    }
    //Rprintf("lk2: lk=%f\n", llkhd);
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
        ell[0] += egxt;
        egxt= exp(egxt);
        ell[0] += log(Sy2[k])+(egxt-1.0)*log(Sy[k]);
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
        //tmp=0.0;
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
    tmp = R_Calloc(d, double);
    Sy = R_Calloc(n, double);
    Sy2 = R_Calloc(n, double);
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
    R_Free(tmp); 
    R_Free(Sy); 
    R_Free(Sy2);
}
/*////////////////////////////////////////////////////*/
/*   Initializing p for fm(.|x1; p), x1=x0(gama1)     */
/*       using fm(.|x0; p0), where x0=x0(gama0)       */
/*            dgx0 = gama1*(x1-x0)                    */
/*////////////////////////////////////////////////////*/
void initialize_p(double *p, int m, double dgx0){
    int i, j, mp1=m+1;
    double pi0=0.0, sum_p=0.0, edgx0, *tmp, *Tmp;
    tmp=R_Calloc(mp1, double);
    Tmp=R_Calloc(mp1, double);
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
        //Tmp[i]+=p[mp1];
    }
    for(i=0; i<=m;i++){ 
        p[i]=edgx0*R_pow(Tmp[i], edgx0-1.0)*tmp[i];
        sum_p+=p[i];
    }
    //Rprintf("Init: sum_p=%f\n",sum_p);
    for(i=0; i<=m; i++){
    p[i]=pi0*p[i]/sum_p;
    //Rprintf("Init: p[i]=%f\n",p[i]);
    }
    p[mp1]=1-pi0;
    R_Free(tmp); 
    R_Free(Tmp); 
}

/*////////////////////////////////////////////////////*/
/* maximizer p of ell(gamma, p) for a gvien gamma     */
/*  egx: exp(gama*x.tilde), where gama is the given   */
/*       regression coefficient                       */
/*////////////////////////////////////////////////////*/
void pofg_ph(double *p, int m, double *egx, int n0, int n1, double *BSy, double *BSy2, 
        double *llik, double eps, int maxit, int prog, int *conv, double *delta){
    int i, j, n=n0+n1, mp1=m+1, mp2=m+2, it=0;
    double sum_egx=0.0,  del=1.0, Sp, Sp2;
    double *Tmp, *Tmp2, *pnu, tmp1, tmp2, llik_nu;
    Tmp=R_Calloc(mp2, double);
    Tmp2=R_Calloc(mp2, double);
    pnu=R_Calloc(mp2, double);
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
    R_Free(Tmp); 
    R_Free(Tmp2); 
    R_Free(pnu);//Free(Sp);  Free(Sp2);
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
  phat = R_Calloc(lp, double);
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  ell = R_Calloc(1, double);
  dell = R_Calloc(d, double);
  BSy = R_Calloc(n*mp2, double);  
  BSy2 = R_Calloc(n*mp2, double); 
  Sy = R_Calloc(n, double);
  Sy2 = R_Calloc(n, double);
  egx = R_Calloc(n, double);
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
  //for(i=0;i<=m;i++) p[i]=*pi0/(double) mp1; 
  //p[mp1]=1.0-*pi0;
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
  R_Free(phat);  
  R_Free(cp);  
  R_Free(res);  
  R_Free(ell);  
  R_Free(dell);
  R_Free(BSy); 
  R_Free(BSy2); 
  R_Free(Sy); 
  R_Free(Sy2); 
  R_Free(egx);  
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
  ell1 = R_Calloc(1, double);  
  dell = R_Calloc(d, double);  
  egx = R_Calloc(n, double);  
  BSy = R_Calloc(n*mp2, double);  
  BSy2 = R_Calloc(n*mp2, double);  
  gnu = R_Calloc(d, double);  
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
  R_Free(BSy); 
  R_Free(BSy2); 
  R_Free(gnu); 
  R_Free(ell1); 
  R_Free(dell); 
  R_Free(egx); 
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
  int m, *cp, tmp, lp, m1=1, cp0=1;//, cp1=1;  
  double *ghat, *phat, *res, *ell, pct, ttl, *lrcp, pv0=1.0;//, pv1=1.0; 
  lp=M[0]*(k+1)+(k+1)*(k+4)/2;
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  phat=R_Calloc(lp, double);
  ghat=R_Calloc(d*(k+1), double);
  ell=R_Calloc(1, double);
  lrcp=R_Calloc(k, double);
  //egx=Calloc(n, double);
  if(*progress==1) {Rprintf("\n Mable fit of Cox PH regression model ... \n");
      ProgressBar(0.0,""); }
  ttl=(double)(k+2)*(k+1);
  m=M[0]; 
  //for(i=0;i<=m;i++) p[i]=*pi0/(double)(m+1);
  //p[m+1] = 1-*pi0;
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
    if(pval[i]<pv0){
      pv0=pval[i];
      cp0=chpts[i];
      m1 = m;
      for(j=0;j<i; j++) lrcp[j]=lr[j];
    }
//    if(chpts[i]>chpts[i-1]){
//      cp1=chpts[i];
//    }
//    if(cp0<cp1) pv1=pval[i];
//    else pv0=pval[i];
//    if(pv1<pv0){
//      cp0=cp1;
//      pv0=pv1;
//    }
//    else pv0=pval[i];
    R_CheckUserInterrupt();
    pct +=2*(i+1)/ttl;
    if(*progress==1) ProgressBar(fmin(1.0,pct),""); 
    i++;
  }
  if(*progress==1){
    ProgressBar(1.0,"");
    Rprintf("\n");}
  if(m==M[1]){
    conv[0]+=1; 
    Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value of the change-point %f is returned.\n", pv0);}
  //else conv[0]=0;
  M[1]=m1;
  tmp=cp0*(M[0]*2+(cp0+3))/2;
  dm[1]=cp0+M[0];
  m=dm[1];
  for(j=0;j<=m+1;j++) p[j]=phat[tmp+j];
  for(j=0; j<dm[0]; j++) gama[j]=ghat[dm[0]*cp0+j];
  mable_ph_m(gama, p, dm, x, y, y2, N, x0,  ell, ddell, EPS, MAXIT, &prg, conv, res);
  for(j=0;j<m1-M[0]; j++) lr[j]=lrcp[j];
  if(*progress==1) Rprintf("\n");
  R_Free(phat); R_Free(ghat); R_Free(ell); 
  //Free(egx);
  R_Free(cp); R_Free(res); R_Free(lrcp);
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
  Tmp=R_Calloc(n*mp2, double);
  Tmp2=R_Calloc(n*mp2, double);
  pnu=R_Calloc(mp2, double);
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
  R_Free(Tmp);  R_Free(Tmp2);  R_Free(pnu);//R_Free(Sp);  Free(Sp2); Free(egx); 
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
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  egx=R_Calloc(n, double);
  phat=R_Calloc(lp, double);
  ell = R_Calloc(1, double);
  BSy = R_Calloc(n*mp2, double);  
  BSy2 = R_Calloc(n*mp2, double); 
  Sy = R_Calloc(n, double);
  Sy2 = R_Calloc(n, double);
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
  R_Free(cp);  
  R_Free(res);  
  R_Free(egx);  
  R_Free(phat);  
  R_Free(ell);  
  R_Free(BSy); 
  R_Free(BSy2); 
  R_Free(Sy); 
  R_Free(Sy2); 
}
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////*/
/*                                                        */
/*                    C Program for                       */
/*  Maximum Approximate Bernstein likelihood Estimation   */
/*   in Proportional Odds Rates Regression model based    */
/*                  Interval Censored data                */
/*////////////////////////////////////////////////////////*/
/*////////////////////////////////////////////////////////////*/
/*                                                            */
/*  Reference:                                                */
/*   Zhong Guan, Maximum Approximate Bernstein Likelihood     */
/*         Estimation in Proportional Odds Rates Model        */
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
/*  MABLE for eta-PO Model:                                            */
/*        [1-S^eta(t|x)]/S^eta(t|x)                                    */
/*        ------------------------- = exp[gamma'x]                     */
/*        [1-S^eta(t|0)]/S^eta(t|0)                                    */
/*   where S(t|x) is the survival function given covariate x.          */
/*  Maximum Approximate Bernstein Likelihood Estimation of survival    */
/*     survival function S(t|x) and regression coefficients gamma      */
/*       based on interval censored data (y=(y1,y2), x, delta)         */
/* (y1,y2): censoring interval containing event time t                 */
/*       x: d-vector of covariate values                               */
/*   delta: censoring indicator, 0 uncensored, 1: interval censored    */
/* The data are arranged so that the 1st n0 obs are uncensored (y1=y2) */
/*  (delta=0) and the rest n1=n-n0 obs are interval censored(delta=1)  */
/*/////////////////////////////////////////////////////////////////////*/
/*/////////////////////////////////////////////////////////////////*/
/*   Log-Likelihood ell(gamma, p) for po model, where              */
/*   gamma=(gamma1,...,gamma_d), p=(p0, p1, ..., pm),              */
/*     egx: exp(gamma*x.tilde)                                     */
/*     BSy: 1-B(y1)                                                */
/*    BSy2: (beta(y1),1-B(y2)),                                    */
/*/////////////////////////////////////////////////////////////////*/
//
double log_blik_po(int m, double *egx, int n0, int n1, 
          double *Sy, double *Sy2, double eta){
  int i, n=n0+n1;
  double llkhd, Stx, Stx2, eta1=1.0/eta;
  llkhd = 0.0;
  for(i=0; i<n0; i++){
    llkhd += log(egx[i]*Sy2[i])-(eta1+1.0)*log(egx[i]+(1.0-egx[i])*R_pow(Sy[i],eta));
  }
  //Rprintf("lk1: lk=%f\n", llkhd);
  for(i=n0; i<n; i++){
    Stx = R_pow(egx[i]+(1.0-egx[i])*R_pow(Sy[i],eta), eta1);
    Stx2 = R_pow(egx[i]+(1.0-egx[i])*R_pow(Sy2[i],eta), eta1);
    llkhd += log(Sy[i]/Stx-Sy2[i]/Stx2);
  }
  //Rprintf("lk2: lk=%f\n", llkhd);
  return llkhd;
}
/*/////////////////////////////////////////////////////////////////*/
/*      Derivatives of loglikelihood ell(gamma, p) wrt gamma       */
/*/////////////////////////////////////////////////////////////////*/
//
void logblik_po_derv(double *gama, int d, double *x, double *x0, double *egx, int n0, int n1,
      double *Sy, double *Sy2, double *ell, double *dell, double *ddell, double eta){
  int i,j,k, n=n0+n1;
  double tmp0=0.0, tmp=0.0, tmp2=0.0;
  double Stx, Stx2, dStx, Stxe, Stxe2, eta1=1.0/eta;
  ell[0]=0.0;
  for(i=0; i<d; i++){ dell[i]=0.0; 
    for(j=0; j<d; j++) ddell[i+d*j]=0.0;}
  for(k=0; k<n0; k++){
    //egxt=0.0;
    //for(i=0; i<d; i++) egxt += gama[i]*(x[k+n*i]-x0[i]);
    ell[0] += log(egx[k]); //egxt;
    //egxt= exp(egxt);
    Stx = R_pow(Sy[k],eta);
    //tmp = egxt+(1.0-egxt)*Stx;
    tmp = egx[k]+(1.0-egx[k])*Stx;
    ell[0] += log(Sy2[k])-(eta1+1.0)*log(tmp);
    Stx2 = Stx/tmp;
    for(i=0; i<d; i++){
      dell[i] += (1.0-(1.0+eta1)*(1.0-Stx2))*(x[k+n*i]-x0[i]); 
      for(j=0; j<d; j++)
        ddell[i+d*j] -= (1.0+eta1)*Stx2*(1.0-Stx2)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
    }
  }
  for(k=n0; k<n; k++){
    //egxt=0.0;
    //for(i=0; i<d; i++) egxt+= gama[i]*(x[k+n*i]-x0[i]);
    //egxt= exp(egxt);
    //Stx=Sy[k]/R_pow(egxt+(1.0-egxt)*R_pow(Sy[k],eta), eta1); 
    //Stx2=Sy2[k]/R_pow(egxt+(1.0-egxt)*R_pow(Sy2[k],eta), eta1);  
    Stx=Sy[k]/R_pow(egx[k]+(1.0-egx[k])*R_pow(Sy[k],eta), eta1); 
    Stx2=Sy2[k]/R_pow(egx[k]+(1.0-egx[k])*R_pow(Sy2[k],eta), eta1);  
    dStx =Stx-Stx2;
    ell[0] += log(dStx);
    Stxe = R_pow(Stx,eta);
    tmp = (1.0-(1.0+eta)*Stxe)*(1.0-Stxe)*Stx*eta1*eta1/dStx;
    if(Sy2[k]>0) Stxe2 = R_pow(Stx2,eta);      
    else Stxe2=0.0;
    tmp2 = (1.0-(1.0+eta)*Stxe2)*(1.0-Stxe2)*Stx2*eta1*eta1/dStx;
    tmp0 = ((Stxe-1.0)*Stx-(Stxe2-1.0)*Stx2)*eta1/dStx;
    for(i=0; i<d; i++){
      dell[i]+=(x[k+n*i]-x0[i])*tmp0;           
      for(j=0;j<d;j++){ 
        ddell[i+d*j]-=tmp0*tmp0*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
        ddell[i+d*j]+=(tmp-tmp2)*(x[k+n*i]-x0[i])*(x[k+n*j]-x0[j]);
      }
    }
  }
}
/*/////////////////////////////////////////////////////////////////*/
/*   Derivatives of loglikelihood ell(gamma, p) wrt (gamma, eta)   */
/*/////////////////////////////////////////////////////////////////*/
// gamma=theta[0:(d-1)], eta=theta[d]
void dllik_gamma_eta(double *theta, int d, double *x, double *x0, int n0, int n1,
      double *Sy, double *Sy2, double *ell, double *dell, double *ddell){
  int i,j,k, n=n0+n1, d1=d+1;
  double egxt, tmp=0.0, tmp2, SlnS, SlnS2, Sgg=0.0, Sgg2=0.0, xti, xtj;
  double Se, Se2, Sg, Sg2, Sge, Sge2, See, See2;
  double Stx, Stx2, dStx, Stxe, Stxe2, eta=theta[d], eta1=1.0/eta;
  ell[0]=0.0;
  for(i=0; i<=d; i++){ dell[i]=0.0; 
    for(j=0; j<=d; j++) ddell[i+d1*j]=0.0;}
  for(k=0; k<n0; k++){
    egxt=0.0;
    for(i=0; i<d; i++) egxt += theta[i]*(x[k+n*i]-x0[i]);
    ell[0] += egxt;
    egxt = exp(egxt);
    Stx = R_pow(Sy[k],eta);
    tmp = egxt+(1.0-egxt)*Stx;
    ell[0] += log(Sy2[k])-(eta1+1.0)*log(tmp);
    Stx2 = Stx/tmp;
    if(Sy[k]>0){
      SlnS=log(Stx)*Stx2;  SlnS2=log(Stx)*log(Stx)*Stx2;}
    else {SlnS=0;SlnS2=0;}
    for(i=0; i<d; i++){
      xti = x[k+n*i]-x0[i];
      dell[i] += (1.0-(1.0+eta1)*(1.0-Stx2))*xti; 
      ddell[i+d1*d]+=(eta1*eta1*((1-Stx2)+(1+eta)*egxt*SlnS/tmp))*xti;
      //ddell[d+d1*i] = ddell[i+d1*d];
      for(j=0; j<d; j++){
        xtj = x[k+n*j]-x0[j];
        ddell[i+d1*j] -= (1.0+eta1)*Stx2*(1.0-Stx2)*xti*xtj;
      }
    }
    dell[d] += eta1*eta1*(log(tmp)+(1.0+eta)*(egxt-1)*SlnS);///????
    ddell[d+d1*d]+=log(tmp)+(egxt-1)*SlnS;
    ddell[d+d1*d]=ddell[d+d1*d]-.5*(eta+1)*egxt*(egxt-1)*SlnS2/tmp;
  }
  ddell[d+d1*d]*= -2*eta1*eta1*eta1;
  //Rprintf("\nn0=%d, n1=%d\n", n0, n1);
  //Rprintf("eta=%f\n", eta);
  //Rprintf("ell=%f\n", ell[0]);
  //Rprintf("dell=%f,%f,%f\n", dell[0], dell[1], dell[d]);
  //Rprintf("ddell=%f,%f,%f,%f\n", ddell[0], ddell[1], ddell[2], ddell[d+d1*d]);
  for(k=n0; k<n; k++){
    egxt=0.0;
    for(i=0; i<d; i++) egxt+= theta[i]*(x[k+n*i]-x0[i]);
    egxt= exp(egxt);
    tmp=egxt+(1.0-egxt)*R_pow(Sy[k],eta);
    Stx=Sy[k]/R_pow(tmp, eta1); 
    tmp2=egxt+(1.0-egxt)*R_pow(Sy2[k],eta);
    Stx2=Sy2[k]/R_pow(tmp2, eta1);  
    dStx =Stx-Stx2;
    ell[0] += log(dStx);
    Stxe = R_pow(Stx,eta);
    Sg=eta1*(Stxe-1.0)*Stx;
    Sgg = -eta1*(1.0-(1.0+eta)*Stxe)*Sg;
    //Se = eta1*Stx*(log(Sy[k])-log(Stx)+(egxt-1)*log(Sy[k])*Stxe);
    Se = eta1*Stx*(eta1*log(tmp)+(egxt-1)*log(Sy[k])*Stxe);
    Sge=eta1*(eta1*Stx*(1-Stxe)-Se+Stxe*(Stx*log(Stx)+(eta+1)*Se));
    See=eta1*Se*(-2+log(Sy[k])*egxt/tmp-log(Stx));
    See+=eta1*Stx*Stxe*log(Sy[k])*log(Sy[k])*egxt*(egxt-1)/tmp;
    if(Sy2[k]>0){
      Stxe2 = R_pow(Stx2,eta); 
      //Se2 = eta1*Stx2*(log(Sy2[k])-log(Stx2)+(egxt-1)*log(Sy2[k])*Stxe2);
      Se2 = eta1*Stx2*(eta1*log(tmp2)+(egxt-1)*log(Sy2[k])*Stxe2);
      Sge2=eta1*(eta1*Stx2*(1-Stxe2)-Se2+Stxe2*(Stx2*log(Stx2)+(eta+1)*Se2));
      See2=eta1*Se2*(-2+log(Sy2[k])*egxt/tmp2-log(Stx2));
      See2+=eta1*Stx2*Stxe2*log(Sy2[k])*log(Sy2[k])*egxt*(egxt-1)/tmp2;
    }     
    else{
      Stxe2=0.0;
      Se2 = 0.0;
      Sge2=0.0;
      See2=0.0;
    }
    Sg2=eta1*(Stxe2-1.0)*Stx2;
    Sgg2 = -eta1*(1.0-(1.0+eta)*Stxe2)*Sg2;
    tmp = (Sg-Sg2)/dStx;
    for(i=0; i<d; i++){
      xti = x[k+n*i]-x0[i];
      dell[i]+=xti*tmp;           
      ddell[i+d1*d]+=xti*((Sge-Sge2)-(Se-Se2)*tmp)/dStx;
      //ddell[d+d1*i] = ddell[i+d1*d];
      for(j=0;j<d;j++){ 
        xtj = x[k+n*j]-x0[j];
        ddell[i+d1*j]+=((Sgg-Sgg2)/dStx-tmp*tmp)*xti*xtj;
      }
    }
    dell[d] += (Se-Se2)/dStx;
    ddell[d+d1*d]+=(See-See2)/dStx-(Se-Se2)*(Se-Se2)/(dStx*dStx);
  }
  for(i=0; i<d; i++) ddell[d+d1*i] = ddell[i+d1*d];
  //Rprintf("ell=%f\n", ell[0]);
  //Rprintf("dell=%f\n", dell[d]);
  //Rprintf("ddell=%f\n", ddell[d+d1*d]);
}

/*///////////////////////////////////////////////////////*/
/*     Maximize ell(gamma, eta, p) with gvien p with     */
/*  f(t|x0) being approximated by Bernstein polynomial   */
/*///////////////////////////////////////////////////////*/
// gamma=theta[0:(d-1)], eta=theta[d]

void geofp_po(double *theta, int d, double *p, int m, double *x, double *x0,   
       int n0, int n1, double *Sy, double *Sy2, double *ell, double *dell, 
       double *ddell, double eps, int maxit, int prog){
  int i, it=0, d1=d+1;
  double *del, *tmp;
  del=R_Calloc(1, double);
  tmp = R_Calloc(d1, double);
  //Rprintf("NT: eta=%f\n", theta[d]);
  dllik_gamma_eta(theta, d, x, x0, n0, n1, Sy, Sy2, ell, dell, ddell);
  del[0]=0.0;
  for(i=0;i<=d;i++) del[0]+=fabs(dell[i]); 
  while(it<maxit && del[0]>eps){
    //Rprintf("NT: m=%d, it=%d, del=%e, eta=%f\n",m,  it, del, theta[d]);
    newton_iter(ddell, d1, dell, theta, del); 
    theta[d] = fmax(.1,theta[d]);
    dllik_gamma_eta(theta, d, x, x0, n0, n1, Sy, Sy2, ell, dell, ddell);
    for(i=0;i<=d;i++) del[0]+=fabs(dell[i]);
    it++;
    R_CheckUserInterrupt();
  }
  if(prog==0) Rprintf("NT: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del[0], ell[0]);
  R_Free(del); 
  R_Free(tmp); 
}

/*//////////////////////////////////////////////////////////*/
/* Maximizer gamma of ell(gamma, eta, p) with gvien (eta,p) */
/* with f(t|x0) being approximated by Bernstein polynomial  */
/*//////////////////////////////////////////////////////////*/

void gofp_po(double *gama, int d, double *p, int m, double *x, double *x0,   
       double *egx, int n0, int n1, double *Sy, double *Sy2, double *ell, 
       double *dell, double *ddell, double eps, int maxit, int prog, double eta){
  int i, it=0, n=n0+n1;
  double *del, *tmp;//, *egx;
  del = R_Calloc(1, double);
  tmp = R_Calloc(d, double);
  //egx = R_Calloc(n, double); 
  //Rprintf("NT: gama=%f\n", gama[0]);
  //egx_x0(gama, d, x, n, egx, x0);//???? add egx to args
  //logblik_po_derv(gama, d, x, x0, egx, n0, n1, Sy, Sy2, ell, dell, ddell, eta);// remove???
  del[0]=0.0;
  for(i=0;i<d;i++) del[0]+=fabs(dell[i]); 
  while(it<maxit && del[0]>eps){
    //Rprintf("NT: m=%d, it=%d, del=%e, gama=%f\n",m,  it, del, gama[0]);
    newton_iter(ddell, d, dell, gama, del); 
    //egx_x0(gama, d, x, n, egx, x0);
    egxmx0(gama, d, x, n, egx, x0);
    logblik_po_derv(gama, d, x, x0, egx, n0, n1, Sy, Sy2, ell, dell, ddell, eta);
    for(i=0;i<d;i++) del[0]+=fabs(dell[i]);
    it++;
    R_CheckUserInterrupt();
  }
  if(prog==0) Rprintf("NT: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del[0], ell[0]);
  R_Free(del); 
  R_Free(tmp); 
  //R_Free(egx);
}


/* Quasi-Newton method for maximizing likelihood with a given p */

/* Extended version of struct opt_struct */
typedef struct mable_po_struct
{
    int m, n0, n1;
    double *x0, *x,  *Sy, *Sy2;//*p,
    double eta;    
} mable_po_struct, *MablePOStruct;

// deviance: minus log-likelihood of PO model to be called by vmmin()
// with par = gama, npar = d;  

static double deviance_po(int npar, double *par, void *ex){
  int i, n, n0;
  double *egx, val, Stx, Stx2, eta1, eta;
  MablePOStruct PS = (MablePOStruct) ex;
  n0 = PS->n0; n = n0+PS->n1; 
  eta = PS->eta; eta1=1.0/eta;
  //mp1=m+1; 
  egx = R_Calloc(n, double);  
  egxmx0(par, npar, PS->x, n, egx, PS->x0);
  val = 0.0;
  for(i=0; i<n0; i++){
    val -= log(egx[i]*PS->Sy2[i])-(eta1+1.0)*log(egx[i]+(1.0-egx[i])*R_pow(PS->Sy[i],eta));
  }
  //Rprintf("lk1: lk=%f\n", val);
  for(i=n0; i<n; i++){
    Stx = R_pow(egx[i]+(1.0-egx[i])*R_pow(PS->Sy[i],eta), eta1);
    Stx2 = R_pow(egx[i]+(1.0-egx[i])*R_pow(PS->Sy2[i],eta), eta1);
    val -= log(PS->Sy[i]/Stx-PS->Sy2[i]/Stx2);
  }
  Rprintf("lk2: lk=%f\n", val);
  R_Free(egx);  
  return val;
}

// derivative of deviance wrt gamma for PO model to be called by vmmin()
// with par = gamma, npar = d; 

static void D_deviance_po(int npar, double *par, double *df, void *ex)
{
  int i, k, n, n0;
  double egxt, tmp0=0.0, tmp=0.0;
  double Stx, Stx2, dStx, Stxe, Stxe2, eta1, eta;
  MablePOStruct PS = (MablePOStruct) ex;
  n0 = PS->n0; n = n0+ PS->n1;
  eta = PS->eta; eta1=1.0/eta;

  for(i=0; i<npar; i++)  df[i]=0.0; 
  for(k=0; k<n0; k++){
    egxt=0.0;
    for(i=0; i<npar; i++) egxt += par[i]*(PS->x[k+n*i]-PS->x0[i]);
    egxt= exp(egxt);
    Stx = R_pow(PS->Sy[k],eta);
    tmp = egxt+(1.0-egxt)*Stx;
    Stx2 = Stx/tmp;
    for(i=0; i<npar; i++){
      df[i] -= (1.0-(1.0+eta1)*(1.0-Stx2))*(PS->x[k+n*i]-PS->x0[i]); 
    }
  }
  for(k=n0; k<n; k++){
    egxt=0.0;
    for(i=0; i<npar; i++) egxt+= par[i]*(PS->x[k+n*i]-PS->x0[i]);
    egxt= exp(egxt);
    Stx=PS->Sy[k]/R_pow(egxt+(1.0-egxt)*R_pow(PS->Sy[k],eta), eta1); 
    Stx2=PS->Sy2[k]/R_pow(egxt+(1.0-egxt)*R_pow(PS->Sy2[k],eta), eta1);  
    dStx =Stx-Stx2;
    Stxe = R_pow(Stx,eta);
    //tmp = (1.0-(1.0+eta)*Stxe)*(1.0-Stxe)*Stx*eta1*eta1/dStx;
    if(PS->Sy2[k]>0) Stxe2 = R_pow(Stx2,eta);      
    else Stxe2=0.0;
    //tmp2 = (1.0-(1.0+eta)*Stxe2)*(1.0-Stxe2)*Stx2*eta1*eta1/dStx;
    tmp0 = ((Stxe-1.0)*Stx-(Stxe2-1.0)*Stx2)*eta1/dStx;
    for(i=0; i<npar; i++){
      df[i]-=(PS->x[k+n*i]-PS->x0[i])*tmp0;           
    }
  }
}

/* end of quasi-Newton method */


/*////////////////////////////////////////////////////*/
/*   Initializing p for fm(.|x1; p), x1=x0(gama1)     */
/*       using fm(.|x0; p0), where x0=x0(gama0)       */
/*         dgx0 = gama1*(x1-x0) for PO model          */
/*////////////////////////////////////////////////////*/
void initialize_p_po(double *p, int m, double dgx0, double eta){
  int i, j, mp1=m+1;
  double pi0=0.0, sum_p=0.0, edgx0, *tmp, *Tmp;
  tmp=R_Calloc(mp1, double);
  Tmp=R_Calloc(mp1, double);
  edgx0=exp(dgx0);
  pi0=1.0-R_pow(p[mp1], edgx0);
  //Rprintf("Init: pi0=%f\n",pi0);
  for(i=0; i<=m;i++){ 
    tmp[i]=0.0;
    Tmp[i]=0.0;
    for(j=0;j<mp1;j++){
      tmp[i]+=p[j]*dbeta(i/(double) m, j+1,m-j+1, FALSE);
      Tmp[i]+=p[j]*(1.0-pbeta(i/(double) m, j+1,m-j+1, TRUE, FALSE)); 
    }
    //Tmp[i]+=p[mp1];
  }
  for(i=0; i<=m;i++){ 
    p[i]=edgx0*tmp[i]/(R_pow(edgx0+(1.0-edgx0)*R_pow(Tmp[i], eta),1.0+1.0/eta));
    sum_p+=p[i];
  }
  //Rprintf("Init: sum_p=%f\n",sum_p);
  for(i=0; i<=m; i++){
    p[i]=pi0*p[i]/sum_p;
    //Rprintf("Init: p[i]=%f\n",p[i]);
  }
  p[mp1]=1-pi0;
  R_Free(tmp); 
  R_Free(Tmp); 
}

/*/////////////////////////////////////////////////////////////*/
/* maximizer p of ell(gamma, eta, p) for a gvien (gamma, eta)  */
/*     egx: exp(gama*x.tilde), where gama is the given         */
/*             regression coefficient                          */
/*/////////////////////////////////////////////////////////////*/
void pofg_po(double *p, int m, double *egx, int n0, int n1, 
      double *BSy, double *BSy2, double *llik, double eps, 
      int maxit, int prog, int *conv, double *delta, double eta){
  int i, j, n=n0+n1, mp1=m+1, mp2=m+2, it=0;
  double del=1.0, Sp, Sp2, lam=0.0, Spe, Spe2, eta1=1.0/eta;
  double *Bp, *Bp2, *pnu, tmp0, tmp, tmp2, llik_nu; 
  Bp=R_Calloc(mp2, double);
  Bp2=R_Calloc(mp2, double);
  pnu=R_Calloc(mp2, double);
  //Rprintf("eta=%f\n",eta);
  while(del>eps && it<maxit){
    for(j=0;j<mp2;j++) pnu[j]=0.0;
    // p = p *psi(p)/lambda(p)
    llik_nu=0.0;
    for(i=0; i<n0;i++){
        llik_nu+=log(egx[i]);
      Sp=0.0; Sp2=0.0; 
      for(j=0;j<mp2;j++){
        Bp[j]=BSy[i+n*j]*p[j];
        Sp+=Bp[j];
        Bp2[j]=BSy2[i+n*j]*p[j];
        Sp2+=Bp2[j];
      }
      Spe = R_pow(Sp, eta);
      Spe2 = R_pow(Sp, eta-1.0);
      llik_nu+=log(Sp2)-(1+eta1)*log(egx[i]+(1.0-egx[i])*Spe);
      for(j=0;j<=mp1;j++){
        pnu[j]+=(eta+1)*(egx[i]-1.0)*Spe2*Bp[j]/(egx[i]+(1.0-egx[i])*Spe);
        pnu[j]+=Bp2[j]/Sp2;
      }
      //pnu[mp1]+=(eta+1)*(egx[i]-1.0)*Spe2*Bp[mp1]/(egx[i]+(1.0-egx[i])*Spe);
      lam+=(egx[i]-1.0)*Spe/(egx[i]+(1-egx[i])*Spe);
    }
    lam = n0 + (eta+1)*lam;
    for(i=n0; i<n;i++){
      Sp=0.0; Sp2=0.0;
      for(j=0;j<mp2;j++){
        Bp[j]=BSy[i+n*j]*p[j];
        Sp+=Bp[j];
        Bp2[j]=BSy2[i+n*j]*p[j];
        Sp2+=Bp2[j];
      }
      Spe = R_pow(Sp, eta);
      Spe2 = R_pow(Sp2, eta);
      tmp=egx[i]+(1-egx[i])*Spe;
      tmp2=egx[i]+(1-egx[i])*Spe2;
      tmp0=Sp/R_pow(tmp, eta1)-Sp2/R_pow(tmp2, eta1);
      tmp=R_pow(tmp, 1.0+eta1);
      tmp2=R_pow(tmp2, 1.0+eta1);
      for(j=0;j<=mp1;j++){
        pnu[j]+=egx[i]*(Bp[j]/tmp-Bp2[j]/tmp2)/tmp0;
      }
      lam+=egx[i]*(Sp/tmp-Sp2/tmp2)/tmp0;
      llik_nu+=log(tmp0);
    }  
    for(j=0;j<=mp1;j++) pnu[j] /= lam;
    del=fabs(llik[0]-llik_nu);
    it++;  llik[0]=llik_nu;
    for(j=0;j<=mp1;j++) p[j]=pnu[j];
    lam=0.0;
    R_CheckUserInterrupt();
  }
  if(prog==0) Rprintf("EM: m=%d, it=%d, del=%e, llik=%f\n",m,  it, del, llik[0]);
  conv[0]=0;
  delta[0]=del;
  if(it==maxit){
    conv[0]+=1;
    //warning("\nThe maximum iterations were reached \nwith a delta = %f.\n", del);
  }
  R_Free(Bp); 
  R_Free(Bp2); 
  R_Free(pnu); 
}
/*////////////////////////////////////////////////////////////*/
/* Maximum Approximate Profile Likelihhod Estimation in       */
/*          PO model with a given (gamma, eta)                */
/* M: set of positive integers as candidate degrees of        */ 
/*       Bernstein poly model                                 */
/* gama: an efficient estimate of regression coefficient      */ 
/*       gamma, for data without covariate we set gama=0      */
/*    x: covariate centered at x0 satisfying                  */
/*             gama'x0=min{gama'xi, i=1,...,n}                */
/*////////////////////////////////////////////////////////////*/
void mable_po_gamma(int *M, double *gama, int *dm, double *x, 
      double *y, double *y2, int *N, double *x0, double *lk, double *lr, 
      double *p, double *ddell, double *eps, int *maxit, int *progress,
      double *level, double *pval, int *chpts, int *conv, double *delta, double *eta){
  int i, j, d=dm[0], k=M[1]-M[0], n0=N[0], n1=N[1], n=n0+n1;
  int *cp, lp, tmp, m=M[1],  mp1=m+1, mp2=m+2, itmp=0, cp0=1, cp1=1; 
  double tini=.000001, pct, ttl, *res, pv0=1.0, pv1=1.0;
  double *ell, *dell, *egx, *BSy, *BSy2, *phat;  
  lp=M[0]*(k+1)+(k+1)*(k+2)/2;
//  lp=M[0]*(k+1)+(k+1)*(k+4)/2;
  phat = R_Calloc(lp, double);
  cp = R_Calloc(1, int);
  res = R_Calloc(1, double);
  ell = R_Calloc(1, double);
  dell = R_Calloc(d, double);
  BSy = R_Calloc(n*mp2, double);  
  BSy2 = R_Calloc(n*mp2, double); 
  //Sy = R_Calloc(n, double);
  //Sy2 = R_Calloc(n, double);
  egx = R_Calloc(n, double);
  if(*progress==1) {
    Rprintf("\n Mable fit of PO model with given regression coefficients ... \n");
    ProgressBar(0.0,""); }
  ttl = (double)((k+2)*(k+1));
  egx_x0(gama, d, x, n, egx, x0);
  //egxmx0(gama, d, x, n, egx, x0);
  // add check to see if any egx is less than 1
  for(i=0;i<n;i++) 
    if(egx[i]<1) {
      Rprintf("\n");
      error("Try another baseline 'x0'.\n");}
  m=M[0]; 
  mp1=m+1;
  mp2=m+2;
  for(i=0;i<=m;i++) p[i]=1.0/(double) mp1; 
  p[mp1]=0.0;
  //p[mp1]=1.0-*pi0;
  if(m>0){
    Bdata(y, m, 0, n, BSy);
    Bdata(y2, m, n0, n1, BSy2);
    pofg_po(p, m, egx, n0, n1, BSy, BSy2, ell, *eps, *maxit, *progress, conv, delta, *eta);
    itmp+=conv[0];
    lk[0]=ell[0]; 
  }
  else{
    lk[0]=0;
    for(i=0;i<n0;i++) lk[0]+=log(egx[i])-(1.0-1.0/(*eta))*log(egx[i]+(1.0-egx[i])*R_pow(1.0-y[i],*eta));
    for(i=n0;i<n;i++) 
      lk[0]+=log((1.0-y[i])/R_pow(egx[i]+(1.0-egx[i])*R_pow(1.0-y[i],*eta),1.0/(*eta))-(1.0-y2[i])/R_pow(egx[i]+(1.0-egx[i])*R_pow(1.0-y2[i],*eta),1.0/(*eta)));
  }
  tmp=mp1;
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
    for(j=0;j<mp1;j++) p[j]=(p[j]+tini/(double) mp2)/(1.0+tini);
    pofg_po(p, m, egx, n0, n1, BSy, BSy2, ell, *eps, *maxit, *progress, conv, res, *eta);
    lk[i]=ell[0];
    for(j=0;j<mp1;j++) phat[j+tmp]=p[j];
    tmp += mp1;
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
      Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with the smallest p-value, %f, of the change-point is returned.\n", pv0);
    }
    delta[0]=res[0];
    delta[1]=pv0;
  }
  M[1]=m;
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  dm[1]=cp0+M[0];
  m=dm[1];
  for(j=0;j<=m;j++) p[j]=phat[tmp+j];
  R_Free(phat);  
  R_Free(cp);  
  R_Free(res);  
  R_Free(ell);  
  R_Free(dell);
  R_Free(BSy); 
  R_Free(BSy2); 
  //R_Free(Sy); 
  //R_Free(Sy2); 
  R_Free(egx);  
}

/*////////////////////////////////////////////////////*/
/*  Maximum approx. Bernstein likelihood estimate of  */
/*           (gamma, p) with a fixed degree m         */
/*////////////////////////////////////////////////////*/

void mable_po_m(double *gama, double *p, int *dm, double *x, double *y, double *y2, 
      int *N, double *x0, double *ell, double *ddell, double *EPS, int *MAXIT,
      int *progress, int *conv, double *delta, double *eta, int *eta_known, int *method){
  int i, j, n0=N[0], n1=N[1], n=n0+n1, d=dm[0], m=dm[1], mp2=m+2, it=0;
  int maxit=MAXIT[0], maxit_em=MAXIT[1], maxit_nt=MAXIT[2], prog=1; 
  double eps=EPS[0], eps_em=EPS[1], eps_nt=EPS[2];
  double tini, *BSy, *BSy2, *Sy, *Sy2, *gnu, del, pct=0.0; 
  double *ell1, *dell, *egx;
  double val = 0.0, abstol=-INFINITY, reltol=eps;
  int *mask, trace=0, fncount = 0, grcount = 0, nREPORT=10;
  int ifail = 0;
  MablePOStruct PS;
  PS = (MablePOStruct) R_alloc(1, sizeof(mable_po_struct));
  tini=0.00001;// tini is used to make sure p is in interior of S(m+1)
  ell1 = R_Calloc(1, double);  
  dell = R_Calloc(d+(*eta_known!=1), double);  
  egx = R_Calloc(n, double);  
  BSy = R_Calloc(n*mp2, double);  
  BSy2 = R_Calloc(n*mp2, double);  
  Sy = R_Calloc(n, double);
  Sy2 = R_Calloc(n, double);
  gnu = R_Calloc(d+(*eta_known!=1), double);  
  Bdata(y, m, 0, n, BSy);
  Bdata(y2, m, n0, n1, BSy2);
  mask = R_Calloc(d, int);
  for (i = 0; i < d; i++) mask[i] = 1;
  //Rprintf("\nx0=%f, gamma=%f\n",x0[0], gama[0]);
  egx_x0(gama, d, x, n, egx, x0);
  //egxmx0(gama, d, x, n, egx, x0);
  //Rprintf("\nx0=%f, gamma=%f\n",x0[0], gama[0]);
  // add check to see if any egx is less than 1
  for(i=0;i<n;i++) 
    if (egx[i]<1){
      Rprintf("\n");
      //Rprintf("egx[%d]=%f\n", i, egx[i]);
      error("Try another baseline 'x0'.\n");}
  for(j=0;j<d;j++)  gnu[j]=gama[j];
  if(*eta_known!=1) gnu[d]=eta[0]; // ?????
  //Rprintf("gnu=%f\n", gnu[0]);
  if(m>0) pofg_po(p, m, egx, n0, n1, BSy, BSy2, ell, eps_em, maxit_em, prog, conv, delta, eta[0]);
  //Rprintf("\n ?  x0=%f, gamma=%f\n",x0[0], gama[0]);
    fm_Sm(p, m, BSy, BSy2, n, Sy, Sy2);

  logblik_po_derv(gama, d, x, x0, egx, n0, n1, Sy, Sy2, ell, dell, ddell, eta[0]);
  if(*method==1){
    if(matrix_singular(ddell, d)==1){
      *method=0;
      //Rprintf("\n Singular Hessian matrix. Use Quasi-Newton Method.");
    }
  }


  PS->m=m; PS->n0=n0; PS->n1=n1; 
  PS->x0=x0; PS->x=x; //PS->p=p; 
  PS->eta=eta[0]; PS->Sy=Sy; PS->Sy2=Sy2;
  
  if(*method==1){
      if(*eta_known!=1)
        geofp_po(gnu, d, p, m, x, x0, n0, n1, Sy, Sy2, ell, dell, ddell, eps_nt, maxit_nt, prog);
      else
        gofp_po(gnu, d, p, m, x, x0, egx, n0, n1, Sy, Sy2, ell, dell, ddell, eps_nt, maxit_nt, prog, eta[0]);
  }
  else {
      vmmin(d, gnu, &val, deviance_po, D_deviance_po, maxit, trace, mask, abstol, 
	      reltol, nREPORT, (void *)PS, &fncount, &grcount, &ifail);
      ell[0]=-val;  
  }
  del=0.0;
  for(i=0;i<d;i++){
    del+=fabs(gnu[i]-gama[i]);
    gama[i]=gnu[i];
  }
  if(*eta_known!=1){ 
    del+=fabs(gnu[d]-eta[0]);
    eta[0]=gnu[d];}
  if(m==0) del=0.0;
  //Rprintf("\n ??  x0=%f, gamma=%f\n",x0[0], gama[0]);
  if(*progress==1){
    Rprintf("\n Mable fit of PO model with a given degree ... \n"); 
    ProgressBar(pct,""); }
  while(it<maxit && del>eps){
    egx_x0(gama, d, x, n, egx, x0);
    //egxmx0(gama, d, x, n, egx, x0);
    //Rprintf("\n ???  x0=%f, gamma=%f\n",x0[0], gama[0]);
    // add check to see if any egx is less than 1
    for(i=0;i<n;i++) 
        if (egx[i]<1) {
            Rprintf("\n");
            //Rprintf("egx[%d]=%f\n", i, egx[i]);
            error("Try another baseline 'x0'.\n");}
    for(i=0;i<=m; i++) p[i]=(p[i]+tini/(double) (m+1))/(1.0+tini); 
    pofg_po(p, m, egx, n0, n1, BSy, BSy2, ell1, eps_em, maxit_em, prog, conv, delta, eta[0]);
    //PS->p=p;
    fm_Sm(p, m, BSy, BSy2, n, Sy, Sy2);
    PS->x0=x0;   
    PS->eta=eta[0]; PS->Sy=Sy; PS->Sy2=Sy2;
///////////////
    logblik_po_derv(gama, d, x, x0, egx, n0, n1, Sy, Sy2, ell, dell, ddell, eta[0]);
    if(matrix_singular(ddell, d)==1) *method=0;
    else *method=1;
/////////////////
    if(*method==1){
        if(*eta_known!=1)
          geofp_po(gnu, d, p, m, x, x0, n0, n1, Sy, Sy2, ell1, dell, ddell, eps_nt, maxit_nt, prog);
        else
          gofp_po(gnu, d, p, m, x, x0, egx, n0, n1, Sy, Sy2, ell1, dell, ddell, eps_nt, maxit_nt, prog, eta[0]);
    }
    else{
      vmmin(d, gnu, &val, deviance_po, D_deviance_po, maxit, trace, mask, abstol, 
	      reltol, nREPORT, (void *)PS, &fncount, &grcount, &ifail);
      ell1[0]=-val;  
    }
    del=0.0;
    for(i=0;i<d;i++){
        del+=fabs(gnu[i]-gama[i]);
        gama[i]=gnu[i];
    }
    if(*eta_known!=1){ 
      del+=fabs(gnu[d]-eta[0]);
      eta[0]=gnu[d];}
    del+=fabs(ell1[0]-ell[0]);
    ell[0]=ell1[0];
//    logblik_po_derv(gama, d, x, x0, egx, n0, n1, Sy, Sy2, ell, dell, ddell, eta[0]);
//    if(matrix_singular(ddell, d)==1) *method=0;
//    else *method=1;
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
    if(del>10*eps) conv[0]+=1;
    //Rprintf("\nThe maximum iterations were reached \nwith a delta = %f.\n", del);
  }
  //Rprintf("conv=%d, it=%d, maxit=%d.\n",conv[0], it, maxit);
  //Rprintf("mable-m: it=%d, del=%f\n", it, del);
  minverse(ddell, d); //Sig=-n*ddell
  R_Free(BSy); 
  R_Free(BSy2); 
  R_Free(Sy); 
  R_Free(Sy2); 
  R_Free(gnu); 
  R_Free(ell1); 
  R_Free(dell); 
  R_Free(egx); 
}

/*////////////////////////////////////////////////////////*/
/*  MABLE of (gamma, p) and an optimal degree m           */
/*                                                        */
/*    M: set of positive integers as candidate degrees    */
/*       of Bernstein poly model                          */
/*   gama: initial value of regression coefficient gamma  */
/*   phat: (m1+2)-vector, first mt+2 are mable of p       */
/*         obtained with gamma=gama_hat                   */
/*     x: covariate                                       */
/*    x0: baseline x s.t. gama'x0=min{gama'xi, i=1,...,n} */
/*////////////////////////////////////////////////////////*/
void mable_po(int *M, double *gama, int *dm, double *p, double *x,   
      double *y, double *y2, int *N, double *x0, double *lk, double *lr, 
      double *ddell, double *EPS, int *MAXIT, int *progress, double *level,
      double *pval, int *chpts, int *conv, double *eta, int *eta_known){
  int d=dm[0], i, j, l, k=M[1]-M[0], prg=1-*progress;
  int m, *cp, tmp, lp, method=1, cp0=0;//, cp1=0;  
  double *ehat, *ghat, *phat, *hess, *pv, *del, *ell, pct, ttl;//, pv0=1.0, pv1=1.0; 
  lp=M[0]*(k+1)+(k+1)*(k+2)/2; //M[0]*(k+1)+(k+1)*(k+4)/2;
  cp = R_Calloc(1, int);
  pv = R_Calloc(1, double);
  del = R_Calloc(1, double);
  phat=R_Calloc(lp, double);
  ghat=R_Calloc(d*(k+1), double);
  ehat=R_Calloc(k+1, double);
  hess=R_Calloc(d*d*(k+1),double); // for hessian matrix ddell
  ell=R_Calloc(1, double);
  //egx=R_Calloc(n, double);
  //Rprintf("\nx0=%f\n",x0[0]);
  //Rprintf("gamma=%f\n", gama[0]);
  //Rprintf("x[1:3]=%f,%f,%f\n",x[0], x[1], x[2]);
  if(*progress==1) {Rprintf("\n Mable fit of PO regression model ... \n");
      ProgressBar(0.0,""); }
  ttl=(double)(k+2)*(k+1);
  m=M[0]; 
  for(i=0;i<=m;i++) p[i]=1.0/(double)(m+1);
  p[m+1] = 0.0;//1-*pi0;
  dm[1]=m;
  mable_po_m(gama, p, dm, x, y, y2, N, x0, ell, ddell, EPS, MAXIT, &prg, conv, del, eta, eta_known, &method);
  for(i=0;i<d;i++) {
    ghat[i]=gama[i];
    for(l=0;l<d;l++) hess[i+d*l]=ddell[i+d*l];
  }
  ehat[0]=eta[0];
  lk[0]=ell[0];
  tmp=m+1;
  for(i=0;i<tmp;i++) phat[i]=p[i];
  pval[0]=1.0;
  chpts[0]=0;
  pct = 2/ttl;
  if(*progress==1) ProgressBar(pct,""); 
  i=1;
  while(i<=k && pval[i-1]>*level){
    p[m+2] = 0.0;//p[m+1];
    p[m+1] = (m+1)*p[m]/(double)(m+2);
    for(j=m; j>=1; j--) p[j] = (p[j-1]*j+p[j]*(m+1-j))/(double)(m+2);
    p[0] = (m+1)*p[0]/(double)(m+2);
    m=M[0]+i;
    dm[1]=m;
    for(j=0;j<=m;j++) p[j]=(p[j]+.000001/(double)(m+1))/1.000001;
    mable_po_m(gama, p, dm, x, y, y2, N, x0,  ell, ddell, EPS, MAXIT, &prg, conv, del, eta, eta_known, &method);
    lk[i]=ell[0];  
    for(j=0;j<=m;j++) phat[j+tmp]=p[j];
    tmp += m+1;
    for(j=0; j<d; j++) {
        ghat[d*i+j]=gama[j];
        for(l=0;l<d;l++) hess[i*d*d+j+d*l]=ddell[j+d*l];
    }
    ehat[i]=eta[0];
    //Rprintf("lk[%d]=%f\n",i, lk[i]);
    if(i>=3){
      cp[0]=i;
      chpt_exp(lk, lr, pv, cp);
      pval[i]=pv[0];
      chpts[i]=cp[0];
    }
    else{            
      pval[i]=1.0;
      chpts[i]=0;
    }
//    if(chpts[i]>chpts[i-1]){
//      cp1=chpts[i];
//    }
//    if(cp0<cp1) pv1=pval[i];
//    else pv0=pval[i];
//    if(pv1<pv0){
//      cp0=cp1;
//      pv0=pv1;
//    }
//    else pv0=pval[i];
    R_CheckUserInterrupt();
    pct +=2*(i+1)/ttl;
    if(*progress==1) ProgressBar(fmin2(1.0,pct),""); 
    i++;
  }
//  Rprintf("\n ok1\n");
  if(*progress==1){
    ProgressBar(1.0,"");
    Rprintf("\n");}
  if(m==M[1]){
    //if(pval[i-1]>0.2) conv[0]+=1; 
    Rprintf("\nThe maximum candidate degree has been reached. \nA model degree with p-value of the change-point %f is returned.\n", pval[i-1]);}
  M[1]=m;
  cp0=chpts[i-1];
  tmp=cp0*(M[0]*2+(cp0+1))/2;
  dm[1]=cp0+M[0];
  //tmp=cp[0]*(M[0]*2+(cp[0]+3))/2;
  //dm[1]=cp[0]+M[0];
  for(j=0;j<=dm[1];j++) p[j]=phat[j+tmp];
  m=dm[1];
//  Rprintf("ok2\n");
//  Rprintf("conv=%d\n",conv[0]);
  for(j=0; j<d; j++) {
    gama[j]=ghat[d*cp0+j];
    for(l=0;l<d;l++) ddell[j+l*d]=hess[cp0*d*d+j+l*d];
  }
//  if(matrix_singular(ddell, d)==1) {
//      method=1;
//      mable_po_m(gama, p, dm, x, y, y2, N, x0,  ell, ddell, EPS, MAXIT, &prg, conv, del, eta, eta_known, &method);
//  }
//Rprintf("ok3\n");
  eta[0]=ehat[cp0];
  pval[i]=del[0];
  if(*progress==1) Rprintf("\n");
  R_Free(phat);
  R_Free(ghat);
  R_Free(ehat);
  R_Free(hess);
  R_Free(ell); 
  //R_Free(egx);
  R_Free(cp);
  R_Free(pv);
  R_Free(del);
}

/*///////////////////////////////////////*/
/*  Parametric PO with Weibull Baseline  */
/*///////////////////////////////////////*/
/*/////////////////////////////////////////////////////////////////*/
/*   loglikelihood ell(theta, sigma, kappa) and  derivatives       */
/*/////////////////////////////////////////////////////////////////*/
// gamma=theta[0:(d-1)], eta=theta[d], sigma=theta[d+1], kappa=theta[d+2]
void dlik_weibull(double *theta, int d, double *x, double *y, double *y2, 
      int n0, int n1, double *ell, double *dell, double *ddell){
  int i,j,k, n=n0+n1, d1=d+1, d2=d+2, d3=d+3;
  double sig=theta[d+1], ka=theta[d+2], eta=theta[d], eta1=1.0/eta, ysigk, ysigk2;
  double egx, tmp=0.0, xi, xj, tmp1, tmp2, tmp3, lnysig;
  double Sg, Sg2, Sgg, Sgg2, Se, Se2, Sge, Sge2, See, See2;
  double Sgs, Sgs2, Ses, Ses2, Sgk, Sgk2;
  double Ss, Ss2, Sk, Sk2, Sss, Sss2, Ssk, Ssk2, Sek, Sek2, Skk, Skk2;
  double Sy, Sy2, Stx, Stx2, dStx, Stxe, Stxe2;
  ell[0]=0.0;
  for(i=0; i<=d2; i++){ dell[i]=0.0; 
    for(j=0; j<=d2; j++) ddell[i+d3*j]=0.0;}
  for(k=0; k<n0; k++){
    Sy=R_pow(y[k]/sig, ka);
    Sy2=exp(-Sy);//S0
    lnysig=log(y[k]/sig);
    egx=0.0;
    for(i=0; i<d; i++) egx += theta[i]*x[k+n*i];
    ell[0] += egx+(ka-1)*lnysig-Sy;
    egx = exp(egx);
    Stx = R_pow(Sy2,eta);
    tmp = egx+(1.0-egx)*Stx;
    ell[0] += -(eta1+1.0)*log(tmp);
    Stx2 = Stx/tmp; //Stx^eta
    tmp1=egx/tmp; //1-(1-egx)*Stx^eta
    for(i=0; i<d; i++){
      xi = x[k+n*i];
      dell[i] += (1.0-(1.0+eta1)*(1.0-Stx2))*xi; 
      ddell[i+d3*d]+=(eta1*eta1*(1-Stx2)-(1+eta1)*egx*Sy*Stx2/tmp)*xi;
      for(j=0; j<d; j++){
        xj = x[k+n*j];
        ddell[i+d3*j] -= (1.0+eta1)*Stx2*(1.0-Stx2)*xi*xj;
      }
      ddell[i+d3*d1]+=(eta+1)*(ka/sig)*Stx2*Sy*tmp1*xi;
      ddell[i+d3*d2]+=(eta+1)*Stx2*Sy*tmp1*lnysig*xi;
    }
    dell[d] += eta1*eta1*log(tmp)-(1.0+eta1)*(egx-1)*Stx2*Sy;//???
    tmp2=Sy*(1-(1+eta)*(1-egx)*Stx2);
    dell[d1]+= (tmp2-1)*ka/sig;
    dell[d2]+= (1-tmp2)*lnysig+1/ka;
    ddell[d+d3*d]+=log(tmp)+(egx-1)*Stx2*log(Stx);
    ddell[d+d3*d]-=.5*(eta+1)*egx*(egx-1)*Stx2*log(Stx)*log(Stx)/tmp;
    ddell[d+d3*d1]+=(egx-1)*(ka/sig)*Sy*(1-(eta+1)*Sy*tmp1)*Stx2;
    ddell[d+d3*d2]-=(egx-1)*Sy*lnysig*(1-(eta+1)*Sy*tmp1)*Stx2;
    tmp3 = eta*(eta+1)*(1-egx)*Sy*Sy*Stx2*tmp1;
    ddell[d1+d3*d1]+= -(ka+1)*(tmp2-1)*ka/sig/sig;
    ddell[d1+d3*d1]+= -R_pow(ka/sig,2)*(1+tmp3);
    ddell[d1+d3*d2]+= (1+1/ka)*(tmp2-1)*ka/sig;
    ddell[d1+d3*d2]+= (ka/sig)*log(y[k]/sig)*(1+tmp3);
    ddell[d2+d3*d2]+= (1+lnysig)*((1-tmp2)*lnysig+1/ka)-(ka+1)/ka/ka-lnysig*(lnysig+1/ka);
    ddell[d2+d3*d2]+= -R_pow(lnysig,2)*tmp3;
  }
  ell[0]+=n0*(log(ka)-log(sig));
  ddell[d+d3*d]*= -2*eta1*eta1*eta1;
  //for(i=0;i<=d2;i++) Rprintf("theta[%d]=%f\n", i, theta[i]);

  //Rprintf("ell=%f\n", ell[0]);
  //Rprintf("dell=%f\n", dell[0]);
  //Rprintf("ddell=%f\n", ddell[d+d1*d]);
  //Rprintf("ddell=%f\n", ddell[d+d1*d]);
  for(k=n0; k<n; k++){
    egx=0.0;
    for(i=0; i<d; i++) egx += theta[i]*x[k+n*i];
    egx= exp(egx);
    ysigk=R_pow(y[k]/sig, ka);
    Sy=exp(-ysigk);
    tmp=egx+(1.0-egx)*R_pow(Sy,eta);
    Stx=Sy/R_pow(tmp, eta1); 
    Stxe = R_pow(Stx,eta);
    tmp1=1-(1.0-egx)*Stxe;
    if(y[k]>0){
      Sg=(Stxe-1.0)*Stx/eta;
      Sgg = -(1.0-(1.0+eta)*Stxe)*Sg*eta1;
      //Sgg = (1.0-(1.0+eta)*Stxe)*(1.0-Stxe)*Stx*eta1*eta1;
      Se = eta1*Stx*(eta1*log(tmp)+(1-egx)*Stxe*ysigk);
      Sge=eta1*(-Sg-Se+Stxe*(Stx*log(Stx)+(eta+1)*Se));
      //See=-eta1*Se*(2+ysigk*tmp1+log(Stx));
      //See+=eta1*Stx*Stxe*ysigk*ysigk*(egx-1)*tmp1;
      See=-eta1*Se*(2+ysigk*egx/tmp+log(Stx));
      See+=eta1*Stx*Stxe*ysigk*ysigk*egx*(egx-1)/tmp;
    
      Ss=(ka/sig)*ysigk*Stx*egx/tmp; 
      tmp2=-eta1*((1+(tmp1-1)*(1+eta))*ysigk+log(Stx)+1);
      lnysig=log(y[k]/sig);
      //Sk= -Stx*tmp1*ysigk*lnysig;
      Sk = -(sig/ka)*Ss*lnysig;
      Sgk = -eta1*(1-(1.0+eta)*Stxe)*Sk;
      Ssk =(1/ka)*Ss+(ka/sig)*(ysigk*(1+(tmp1-1)*(1+eta))-1)*Sk;
      Skk =-(Sk+sig*Ssk*lnysig)/ka;
      Sgs = -eta1*(1-(1.0+eta)*Stxe)*Ss;
      //Ses =Ss*tmp2+eta1*Stx*tmp1*(ka/sig)*ysigk;
      Ses =(Se*(1+(tmp1-1)*(1+eta))+Stx*Stxe*log(Stx))*(ka/sig)*ysigk;
      Sss =((ka/sig)*ysigk*(1+(tmp1-1)*(1+eta))-(ka+1)/sig)*Ss;
      //Sek =Sk*tmp2-eta1*Stx*tmp1*ysigk*lnysig*(egx!=1);
      Sek =-(sig/ka)*Ses*lnysig;
    }
    else{
      Sg=0; Sgg=0; Se=0, Sge=0; See=0; Ss=0; Sk=0; Sgk = 0; 
      Sek =0; Ssk =0; Skk =0; Sgs=0; Ses=0; Sss=0;
    }
    //Rprintf("Stx=%f\n", Stx);
    if(y2[k]>0){
      ysigk2=R_pow(y2[k]/sig, ka);
      Sy2=exp(-ysigk2);
      tmp=egx+(1.0-egx)*R_pow(Sy2,eta);
      Stx2=Sy2/R_pow(tmp, eta1);  
      Stxe2 = R_pow(Stx2,eta); 
      tmp1=1-(1-egx)*Stxe2;
      Sg2=(Stxe2-1.0)*Stx2/eta;
      Sgg2 = -(1.0-(1.0+eta)*Stxe2)*Sg2*eta1;
      //Sgg2 = (1.0-(1.0+eta)*Stxe2)*(1.0-Stxe2)*Stx2*eta1*eta1;
      Se2 = eta1*Stx2*(eta1*log(tmp)+(1-egx)*Stxe2*ysigk2);
      Sge2=eta1*(-Sg2-Se2+Stxe2*(Stx2*log(Stx2)+(eta+1)*Se2));
      //See2=-eta1*Se2*(2+ysigk2*tmp1+log(Stx2));
      //See2+=eta1*Stx2*Stxe2*ysigk2*ysigk2*(egx-1)*tmp1;
      See2=-eta1*Se2*(2+ysigk2*egx/tmp+log(Stx2));
      See2+=eta1*Stx2*Stxe2*ysigk2*ysigk2*egx*(egx-1)/tmp;
      
      Ss2=(ka/sig)*ysigk2*Stx2*egx/tmp;
      tmp2=-eta1*((1+(tmp1-1)*(1+eta))*ysigk2+log(Stx2)+1);
      lnysig=log(y2[k]/sig);
      Sk2 = -(sig/ka)*Ss2*lnysig;
      Sgk2 = -eta1*(1-(1.0+eta)*Stxe2)*Sk2;
      Ssk2 =(1/ka)*Ss2+(ka/sig)*(ysigk2*(1+(tmp1-1)*(1+eta))-1)*Sk2;
      Skk2 =-(Sk2+sig*Ssk2*lnysig)/ka;
      Sgs2 = -eta1*(1-(1.0+eta)*Stxe2)*Ss2;
      //Ses2 =Ss2*tmp2+eta1*Stx2*tmp1*(ka/sig)*ysigk2*(egx!=1);
      Ses2 =(Se2*(1+(tmp1-1)*(1+eta))+Stx2*Stxe2*log(Stx2))*(ka/sig)*ysigk2;
      Sss2 =((ka/sig)*ysigk2*(1+(tmp1-1)*(1+eta))-(ka+1)/sig)*Ss2;
      //Sek2 =Sk2*tmp2-eta1*Stx2*tmp1*ysigk2*lnysig*(egx!=1);
      Sek2 =-(sig/ka)*Ses2*lnysig;
    }
    else{  
      Stx2=0; Stxe2 =0; Sg2=0; Sgg2 =0; Se2 =0; Sge2=0; See2=0; Sk2=0; Sgs2 =0;
      Sgk2 =0; Ses2 =0; Sek2 =0; Sss2 =0; Ssk2 =0; Skk2 =0; Ss2=0;
    }
    dStx =Stx-Stx2;  
    //Rprintf("Stx=%f, Stx2=%f, dStx=%f\n", Stx, Stx2, dStx); 
    ell[0] += log(dStx);
    tmp = (Sg-Sg2)/dStx;
    for(i=0; i<d; i++){
      xi = x[k+n*i];
      dell[i]+=xi*tmp;           
      ddell[i+d3*d]+=xi*((Sge-Sge2)-(Se-Se2)*tmp)/dStx;
      ddell[i+d3*d1]+=xi*((Sgs-Sgs2)-(Ss-Ss2)*tmp)/dStx;
      ddell[i+d3*d2]+=xi*((Sgk-Sgk2)-(Sk-Sk2)*tmp)/dStx;
      for(j=0;j<d;j++){ 
        xj = x[k+n*j];
        ddell[i+d3*j]+=((Sgg-Sgg2)/dStx-tmp*tmp)*xi*xj;        
      }
    }
    dell[d] += (Se-Se2)/dStx;
    dell[d1]+= (Ss-Ss2)/dStx;
    dell[d2]+= (Sk-Sk2)/dStx;
    ddell[d+d3*d]+=(See-See2)/dStx-(Se-Se2)*(Se-Se2)/(dStx*dStx);
    ddell[d+d3*d1]+=(Ses-Ses2)/dStx-(Se-Se2)*(Ss-Ss2)/(dStx*dStx);
    ddell[d+d3*d2]+=(Sek-Sek2)/dStx-(Se-Se2)*(Sk-Sk2)/(dStx*dStx);
    ddell[d1+d3*d1]+=(Sss-Sss2)/dStx-(Ss-Ss2)*(Ss-Ss2)/(dStx*dStx);
    ddell[d1+d3*d2]+=(Ssk-Ssk2)/dStx-(Ss-Ss2)*(Sk-Sk2)/(dStx*dStx);
    ddell[d2+d3*d2]+=(Skk-Skk2)/dStx-(Sk-Sk2)*(Sk-Sk2)/(dStx*dStx);
  }
  for(i=0; i<d; i++){
    ddell[d+d3*i] = ddell[i+d3*d];
    ddell[d1+d3*i]= ddell[i+d3*d1];
    ddell[d2+d3*i]= ddell[i+d3*d2];
  }
  ddell[d1+d3*d]= ddell[d+d3*d1];
  ddell[d2+d3*d]= ddell[d+d3*d2];
  ddell[d2+d3*d1]=ddell[d1+d3*d2];
  
  
  //Rprintf("ell=%f\n", ell[0]);
  //for(i=0;i<=d2;i++) {
  //  Rprintf("dell[%d]=%f\n", i, dell[i]);
  //  for(j=0;j<=d2;j++) Rprintf("ddell[%d,%d]=%f\n",i,j, ddell[i+d3*j]); 
  //}
  //R_Free(Sy); R_Free(Sy2);
}
// gamma=theta[0:(d-1)], sigma=theta[d], kappa=theta[d+1]
// eta is known
void dlik_weibull_eta(double *theta, int d, double *x, double *y, double *y2, 
      int n0, int n1, double *ell, double *dell, double *ddell, double eta){
  int i,j,k, n=n0+n1, d1=d+1, d2=d+2;
  double sig=theta[d], ka=theta[d1], eta1=1.0/eta, ysigk, ysigk2;
  double egx, tmp=0.0, xi, xj, tmp1, tmp2, tmp3, lnysig;
  double Sg, Sg2, Sgg, Sgg2, Sgs, Sgs2, Sgk, Sgk2;
  double Ss, Ss2, Sk, Sk2, Sss, Sss2, Ssk, Ssk2, Skk, Skk2;
  double Sy, Sy2, Stx, Stx2, dStx, Stxe, Stxe2;
  ell[0]=0.0;
  for(i=0; i<d2; i++){ dell[i]=0.0; 
    for(j=0; j<d2; j++) ddell[i+d2*j]=0.0;}
  for(k=0; k<n0; k++){
    Sy=R_pow(y[k]/sig, ka);
    Sy2=exp(-Sy);//S0
    lnysig=log(y[k]/sig);
    egx=0.0;
    for(i=0; i<d; i++) egx += theta[i]*x[k+n*i];
    ell[0] += egx+(ka-1)*lnysig-Sy;
    egx = exp(egx);
    Stx = R_pow(Sy2,eta);
    tmp = egx+(1.0-egx)*Stx;
    ell[0] += -(eta1+1.0)*log(tmp);
    Stx2 = Stx/tmp; //Stx^eta
    tmp1=egx/tmp; //1-(1-egx)*Stx^eta
    for(i=0; i<d; i++){
      xi = x[k+n*i];
      dell[i] += (1.0-(1.0+eta1)*(1.0-Stx2))*xi; 
      for(j=0; j<d; j++){
        xj = x[k+n*j];
        ddell[i+d2*j] -= (1.0+eta1)*Stx2*(1.0-Stx2)*xi*xj;
      }
      ddell[i+d2*d]+=(eta+1)*(ka/sig)*Stx2*Sy*tmp1*xi;
      ddell[i+d2*d1]+=(eta+1)*Stx2*Sy*tmp1*lnysig*xi;
    }
    tmp2=Sy*(1-(1+eta)*(1-egx)*Stx2);
    dell[d]+= (tmp2-1)*ka/sig;
    dell[d1]+= (1-tmp2)*lnysig+1/ka;
    tmp3 = eta*(eta+1)*(1-egx)*Sy*Sy*Stx2*tmp1;
    ddell[d+d2*d]+= -(ka+1)*(tmp2-1)*ka/sig/sig;
    ddell[d+d2*d]+= -R_pow(ka/sig,2)*(1+tmp3);
    ddell[d+d2*d1]+= (1+1/ka)*(tmp2-1)*ka/sig;
    ddell[d+d2*d1]+= (ka/sig)*log(y[k]/sig)*(1+tmp3);
    ddell[d1+d2*d1]+= (1+lnysig)*((1-tmp2)*lnysig+1/ka)-(ka+1)/ka/ka-lnysig*(lnysig+1/ka);
    ddell[d1+d2*d1]+= -R_pow(lnysig,2)*tmp3;
  }
  ell[0]+=n0*(log(ka)-log(sig));
  //for(i=0;i<d2;i++) Rprintf("theta[%d]=%f\n", i, theta[i]);

  //Rprintf("ell=%f\n", ell[0]);
  //Rprintf("dell=%f\n", dell[0]);
  //Rprintf("ddell=%f\n", ddell[d+d1*d]);
  //Rprintf("ddell=%f\n", ddell[d+d1*d]);
  for(k=n0; k<n; k++){
    egx=0.0;
    for(i=0; i<d; i++) egx += theta[i]*x[k+n*i];
    egx= exp(egx);
    ysigk=R_pow(y[k]/sig, ka);
    Sy=exp(-ysigk);
    tmp=egx+(1.0-egx)*R_pow(Sy,eta);
    Stx=Sy/R_pow(tmp, eta1); 
    Stxe = R_pow(Stx,eta);
    tmp1=1-(1.0-egx)*Stxe;
    if(y[k]>0){
      Sg=(Stxe-1.0)*Stx/eta;
      Sgg = -(1.0-(1.0+eta)*Stxe)*Sg*eta1;
    
      Ss=(ka/sig)*ysigk*Stx*egx/tmp; 
      tmp2=-eta1*((1+(tmp1-1)*(1+eta))*ysigk+log(Stx)+1);
      lnysig=log(y[k]/sig);
      //Sk= -Stx*tmp1*ysigk*lnysig;
      Sk = -(sig/ka)*Ss*lnysig;
      Sgk = -eta1*(1-(1.0+eta)*Stxe)*Sk;
      Ssk =(1/ka)*Ss+(ka/sig)*(ysigk*(1+(tmp1-1)*(1+eta))-1)*Sk;
      Skk =-(Sk+sig*Ssk*lnysig)/ka;
      Sgs = -eta1*(1-(1.0+eta)*Stxe)*Ss;
      Sss =((ka/sig)*ysigk*(1+(tmp1-1)*(1+eta))-(ka+1)/sig)*Ss;
    }
    else{
      Sg=0; Sgg=0;  Ss=0; Sk=0; Sgk = 0; 
      Ssk =0; Skk =0; Sgs=0; Sss=0;
    }
    //Rprintf("Stx=%f\n", Stx);
    if(y2[k]>0){
      ysigk2=R_pow(y2[k]/sig, ka);
      Sy2=exp(-ysigk2);
      tmp=egx+(1.0-egx)*R_pow(Sy2,eta);
      Stx2=Sy2/R_pow(tmp, eta1);  
      Stxe2 = R_pow(Stx2,eta); 
      tmp1=1-(1-egx)*Stxe2;
      Sg2=(Stxe2-1.0)*Stx2/eta;
      Sgg2 = -(1.0-(1.0+eta)*Stxe2)*Sg2*eta1;
      
      Ss2=(ka/sig)*ysigk2*Stx2*egx/tmp;
      tmp2=-eta1*((1+(tmp1-1)*(1+eta))*ysigk2+log(Stx2)+1);
      lnysig=log(y2[k]/sig);
      Sk2 = -(sig/ka)*Ss2*lnysig;
      Sgk2 = -eta1*(1-(1.0+eta)*Stxe2)*Sk2;
      Ssk2 =(1/ka)*Ss2+(ka/sig)*(ysigk2*(1+(tmp1-1)*(1+eta))-1)*Sk2;
      Skk2 =-(Sk2+sig*Ssk2*lnysig)/ka;
      Sgs2 = -eta1*(1-(1.0+eta)*Stxe2)*Ss2;
      Sss2 =((ka/sig)*ysigk2*(1+(tmp1-1)*(1+eta))-(ka+1)/sig)*Ss2;
    }
    else{  
      Stx2=0; Stxe2 =0; Sg2=0; Sgg2 =0; Sk2=0; Sgs2 =0;
      Sgk2 =0; Sss2 =0; Ssk2 =0; Skk2 =0; Ss2=0;
    }
    dStx =Stx-Stx2;  
    //Rprintf("Stx=%f, Stx2=%f, dStx=%f\n", Stx, Stx2, dStx); 
    ell[0] += log(dStx);
    tmp = (Sg-Sg2)/dStx;
    for(i=0; i<d; i++){
      xi = x[k+n*i];
      dell[i]+=xi*tmp;           
      ddell[i+d2*d]+=xi*((Sgs-Sgs2)-(Ss-Ss2)*tmp)/dStx;
      ddell[i+d2*d1]+=xi*((Sgk-Sgk2)-(Sk-Sk2)*tmp)/dStx;
      for(j=0;j<d;j++){ 
        xj = x[k+n*j];
        ddell[i+d2*j]+=((Sgg-Sgg2)/dStx-tmp*tmp)*xi*xj;        
      }
    }
    dell[d]+= (Ss-Ss2)/dStx;
    dell[d1]+= (Sk-Sk2)/dStx;
    ddell[d+d2*d]+=(Sss-Sss2)/dStx-(Ss-Ss2)*(Ss-Ss2)/(dStx*dStx);
    ddell[d+d2*d1]+=(Ssk-Ssk2)/dStx-(Ss-Ss2)*(Sk-Sk2)/(dStx*dStx);
    ddell[d1+d2*d1]+=(Skk-Skk2)/dStx-(Sk-Sk2)*(Sk-Sk2)/(dStx*dStx);
  }
  for(i=0; i<d; i++){
    ddell[d+d2*i]= ddell[i+d2*d];
    ddell[d1+d2*i]= ddell[i+d2*d1];
  }
  ddell[d1+d2*d]=ddell[d+d2*d1];
  
  
  //Rprintf("ell=%f\n", ell[0]);
  //for(i=0;i<d2;i++) {
  //  Rprintf("dell[%d]=%f\n", i, dell[i]);
  //  for(j=0;j<d2;j++) Rprintf("ddell[%d,%d]=%f\n",i,j, ddell[i+d2*j]); 
  //}
  //R_Free(Sy); R_Free(Sy2);
}

/*///////////////////////////////*/
/*  Newton method for GPO model  */
/*    with Weibull baseline      */
/*///////////////////////////////*/
// if eta is unknown, theta=(gamma,eta,sigma,kappa)
// if eta is known, theta=(gamma,sigma,kappa,eta)
void weib_gpo(double *theta, int *d, double *x, int *n0, int *n1, double *y,     
       double *y2, double *lk, double *ddell, double *eps, 
       int *maxit, int *prog, int *conv, double *del, int *eta_known){
  int i, it=0, d2=*d+2, np=d2+(*eta_known!=1);
  double *dell, eta;
  dell = R_Calloc(np, double);
  if(*eta_known!=1){
    dlik_weibull(theta, *d, x, y, y2, *n0, *n1, lk, dell, ddell);
    del[0]=0.0;
    for(i=0;i<=d2;i++) del[0]+=fabs(dell[i]); 
    //Rprintf("NT: eta=%f\n", theta[*d]);
    //Rprintf("NT: it=%d, del=%f, llik=%f\n", it, del, lk[0]);
    while(it<*maxit && del[0]>*eps){
      //Rprintf("NT: it=%d, del=%e, eta=%f\n", it, del[0], theta[*d]);
      //minverse(ddell, d3);  
      newton_iter(ddell, np, dell, theta, del); 
      theta[*d] = fmax(0.01, theta[*d]);
      dlik_weibull(theta, *d, x, y, y2, *n0, *n1, lk, dell, ddell);
      for(i=0;i<=d2;i++) del[0]+=fabs(dell[i]);
      it++;
      R_CheckUserInterrupt();
    }
  }
  else{
    eta=theta[d2];
    dlik_weibull_eta(theta, *d, x, y, y2, *n0, *n1, lk, dell, ddell, eta);
    del[0]=0.0;
    for(i=0;i<d2;i++) del[0]+=fabs(dell[i]); 
    //Rprintf("NT: eta=%f\n", theta[*d]);
    //Rprintf("NT: it=%d, del=%f, llik=%f\n", it, del, lk[0]);
    while(it<*maxit && del[0]>*eps){
      //Rprintf("NT: it=%d, del=%e, eta=%f\n", it, del[0], theta[*d]);
      newton_iter(ddell, np, dell, theta, del); 
      dlik_weibull_eta(theta, *d, x, y, y2, *n0, *n1, lk, dell, ddell, eta);
      for(i=0;i<d2;i++) del[0]+=fabs(dell[i]);
      it++;
      R_CheckUserInterrupt();
    }
  }  
  if(it<*maxit) conv[0]=0;
  else conv[0]=1;
  if(*prog==1) Rprintf("NT: it=%d, del=%f, llik=%f\n", it, del[0], lk[0]);
  R_Free(dell);
  //R_Free(A); R_Free(b);
}

