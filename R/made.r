#################################################################
#         Minimum Approximate Distance Estimation
#################################################################


#' Exponential change-point
#' @param x  a vector of nondecreasing values of log-likelihood or -log of distance
#' @return a list of exponential change-point, its p-value, and the likelihood ratios 
#'           of the exponential change-point model
#' @export
chpt.exp=function(x){
    lr=NULL
    k = length(x)-1
    maxLR=0.0
    lnk=log(k)
    llnk=log(lnk)
    lr0 = k*log((x[k+1]-x[1])/k)
    lr[k]=0.0
    for(i in 1:(k-1)){
        lr[i] = lr0-i*log((x[i+1]-x[1])/i)-(k-i)*log((x[k+1]-x[i+1])/(k-i))
        if(lr[i]>maxLR){
            chpt = i+1
            maxLR=lr[i]
        }
    }
    # p-value of change-point
    pv=1.0-exp(-2*lnk*lnk*sqrt(llnk*pi)*exp(-2*sqrt(maxLR*llnk)))
    list(lr=lr, chpt=chpt, pvalue=pv)
}

# Calculate pt of degree mt=m+e_t from p of degree m
pm2pmpe_k<-function(p, m, t){
  d<-length(m)
  e<-diag(1,d)
  I<-rep(0,d)
  mt<-m+e[,t]
  km<-c(1,cumprod(m+1))
  kmt<-c(1,cumprod(mt+1))
  K<-km[d+1];
  K1<-kmt[d+1];
  p1<-rep(0,K1)
#  it<-K-1
#  while(it>=0){
  it<-0
  while(it<K){
    r <- it
    k<-d-1
    while(k > 0){
      jj<-r%%km[k+1]
      I[k+1]<-(r-jj)/km[k+1]
      r<-jj
      k<-k-1
    }
    I[1]<-r
    jt1<-1+sum((I+e[,t])*kmt[1:d])
    jt2<-1+sum(I*kmt[1:d])
    it<-it+1
    #if(d<=2) cat("it=",it,"jt1=",jt1,"jt2=",jt2,"I=",I,"\n")
    p1[jt1]<-p1[jt1]+(I[t]+1)*p[it]/(mt[t]+1)
    p1[jt2]<-p1[jt2]+(mt[t]-I[t])*p[it]/(mt[t]+1)
    # p1[jt1]<-p1[jt1]+(I[t]+1)*p[it+1]/(mt[t]+1)
    # p1[jt2]<-p1[jt2]+(mt[t]-I[t])*p[it+1]/(mt[t]+1)
    # it<-it-1
  }
  p1
}

#' Multivariate empirical cumulative distribution evaluated at sample data
#' @param x an \code{n x d} matrix of data values, 
#'     rows are \code{n} observations of \eqn{X=(X_1,...,X_d)}
#' @return a vector of \code{n} values 
#' @export
mvecdf<-function(x){
  n<-nrow(x)
  d<-ncol(x)
  ans<-rep(0,n)
  for(i in 1:n)
    #ans[i]<-mean(apply(x, 1, function(y) all(y<=x[i,]))) 
    for(j in 1:n) ans[i]<-ans[i]+(all(x[j,]<=x[i,]))
  ans/n
}

#' Component Beta cumulative distribution functions of the Bernstein polynomial model
#' @param x  \code{n x d} matrix, rows are n observations of X=(X1,...,Xd)
#' @param m  vector of d nonnegative integers \code{m=(m[1],...,m[d])}.
#' @return an \code{n x K} matrix, \code{K=(m[1]+1)...(m[d]+1)}.
#' @importFrom stats pbeta
#' @export
mvpbeta<-function(x, m){
  n<-nrow(x)
  d<-ncol(x)
  km<-c(1,cumprod(m+1))
  K<-km[d+1]
  dm<-matrix(1, nrow=n, ncol=K)
  it<-0
  while(it<K){
    r <- it
    for(k in (d-1):1){
      jj <- r%%km[k+1]
      i <- (r-jj)/km[k+1]
      dm[,it+1]<-dm[,it+1]*pbeta(x[,k+1], i+1, m[k+1]-i+1)
      r <- jj
    }
    dm[,it+1]<-dm[,it+1]*pbeta(x[,1], r+1, m[1]-r+1)
    it<-it+1
  }
  dm
}


#' Minimum Approximate Distance Estimate of univariate Density Function 
#'  with given model degree(s)
# 
#' @param x an \code{n x d} matrix or \code{data.frame} of multivariate sample of size \code{n}
#' @param m a positive integer or a vector of \code{d} positive integers specify  
#'    the given model degrees for the joint density.
#' @param p initial guess of \code{p}
#' @param interval a vector of two endpoints or a \code{2 x d} matrix, each column containing 
#'    the endpoints of support/truncation interval for each marginal density.
#'    If missing, the i-th column is assigned as \code{c(min(x[,i]), max(x[,i]))}.
#' @param method method for finding minimum distance estimate. "em": EM like method;
#     "qp": the quadratic programming algorithm using packages \code{quadprog} or \code{LowRankQP}.
#' @param maxit the maximum iterations
#' @param eps the criterion for convergence
#' @details 
#'   A \eqn{d}-variate cdf \eqn{F} on a hyperrectangle \eqn{[a, b]
#'   =[a_1, b_1] \times \cdots \times [a_d, b_d]} can be approximated 
#'   by a mixture of \eqn{d}-variate beta cdfs on \eqn{[a, b]}, 
#'   \eqn{\beta_{mj}(x) = \prod_{i=1}^dB_{m_i,j_i}[(x_i-a_i)/(b_i-a_i)]},
#'   with proportion \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}. 
#'   With a given model degree \code{m}, the parameters \code{p}, the mixing
#'   proportions of the beta distribution, are calculated as the minimizer of the
#'   approximate \eqn{L_2} distance between the empirical distribution and 
#'   the Bernstein polynomial model. The quadratic programming with linear constraints
#'   is used to solve the problem.
#' @return An invisible \code{mable} object with components
#' \itemize{
#'   \item \code{m} the given model degree(s) 
#'   \item \code{p} the estimated vector of mixture proportions 
#'       with the given optimal degree(s) \code{m}
#'   \item \code{interval} support/truncation interval \code{[a, b]}
#'   \item \code{D}  the minimum distance at degree \code{m}
#'  }
#' @importFrom quadprog solve.QP
#' @importFrom LowRankQP LowRankQP
#' @export
madem.density<-function(x, m, p=rep(1,prod(m+1))/prod(m+1), interval=NULL, 
        method=c("qp","em"), maxit=10000, eps=1e-7){
  data.name<-deparse(substitute(x))
  n<-NROW(x)
  d<-NCOL(x)
  xNames<-names(x)
  method <- match.arg(method)
  if(is.null(xNames)) 
    for(i in 1:d) xNames[i]<-paste("x",i,sep='')
  x<-as.matrix(x)   
  if(is.null(interval))  
    interval<-apply(x, 2, range) 
  else if(is.matrix(interval)){
    if(nrow(interval)!=2 || ncol(interval)!=d)
      stop("Invalid 'interval'.\n")
  }
  else{
    if(length(interval)!=2) 
      stop("Invalid 'interval'.\n")
    else if(any(apply(x, 2, min)<interval[1])||any(apply(x, 2, max)>interval[2]))
      stop("Invalid 'interval'.\n")
    else interval<-matrix(rep(interval, d), nrow=2)
  }
  x<-t((t(x)-interval[1,])/apply(interval,2,diff))
  if(any(x<0) || any(x>1)) stop("All data values must be contained in 'interval'.")
  y<-mvecdf(x)
  ssy<-sum(y^2)
  if(all(m==0)){ # uniform 
    p=1
    Dm=mean((apply(x,1,prod)-y)^2)
  }
  else{
    K=prod(m+1) 
    B=mvpbeta(x,m)
    if(method=="em"){
        C<-B-y
        CC<-t(C)%*%C
        del<-1; it<-0
        while(it<maxit && del>eps){
            Dp<-as.numeric(t(p)%*%CC%*%p)
            pt<-p
            p<-p*CC%*%p^2/Dp
            p<-p^2
            del<-sum(abs(p-pt))
            it<-it+1
        } 
        Dm<-Dp 
    }
    else{
        # equality constraints
        bvec=1 # sum(p)=1
        Amat=rep(1, K) # sum(p)=1
        # add nonnegative constraints
        Amat=rbind(Amat, diag(1,K))
        bvec=c(bvec,rep(0, K))
        Q=diag(0, K)
        for(i in 1:n) Q=Q+B[i,]%*%t(B[i,])
        cvec=t(B)%*%y
        
        #cat("dim A=",dim(Amat),"\n")
        tryCatch({res<-solve.QP(Dmat=Q, dvec=cvec, Amat=t(Amat), bvec=bvec, meq=1); 
                p=res$solution; Dm=(2*res$value+ssy)/n},
            error=function(e){
              if(geterrmessage()
                  =="matrix D in quadratic function is not positive definite!"){
    #            suppressMessages(require(LowRankQP))
                res<-LowRankQP(Vmat=Q, dvec=-cvec, Amat=matrix(Amat[1,], ncol=K), 
                                bvec=bvec[1], uvec=rep(1, K), method="CHOL")
                p<<-res$alpha
                Dm<<-(t(res$alpha)%*%Q%*%res$alpha-2*t(cvec)%*%res$alpha+ssy)/n
              }
              else stop(geterrmessage())
            }
        )
        #p=res$solution; D=(2*res$value+d)/n
    }
  }
  p<-zapsmall(p, digits = 6)
  ans=list(m=m, p=p, interval=interval, D=Dm)
  ans$xNames<-xNames
  if(d>1) ans$data.type<-"mvar"
  class(ans)<-"mable"
  invisible(ans)
}
# Choosing optimal degree
#' Minimum Approximate Distance Estimate of Density Function with an optimal model degree
#' @param x an \code{n x d} matrix or \code{data.frame} of multivariate sample of size \code{n}
#' @param M0 a positive integer or a vector of \code{d} positive integers specify
#'    starting candidate degrees for searching optimal degrees.  
#' @param M a positive integer or a vector of \code{d} positive integers specify  
#'    the maximum candidate or the given model degrees for the joint density.
#' @param search logical, whether to search optimal degrees between \code{M0} and \code{M} 
#'    or not but use \code{M} as the given model degrees for the joint density.
#' @param interval a vector of two endpoints or a \code{2 x d} matrix, each column containing 
#'    the endpoints of support/truncation interval for each marginal density.
#'    If missing, the i-th column is assigned as \code{c(min(x[,i]), max(x[,i]))}.
#' @param mar.deg logical, if TRUE, the optimal degrees are selected based  
#'    on marginal data, otherwise, the optimal degrees are chosen the joint data. See details.   
#' @param method method for finding minimum distance estimate. "em": EM like method;
#     "qp": the quadratic programming algorithm using packages \code{quadprog} or \code{LowRankQP}.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.
#' @param progress if TRUE a text progressbar is displayed
#' @details 
#'   A \eqn{d}-variate cdf \eqn{F} on a hyperrectangle \eqn{[a, b]
#'   =[a_1, b_1] \times \cdots \times [a_d, b_d]} can be approximated 
#'   by a mixture of \eqn{d}-variate beta cdfs on \eqn{[a, b]}, 
#'   \eqn{\beta_{mj}(x) = \prod_{i=1}^dB_{m_i,j_i}[(x_i-a_i)/(b_i-a_i)]},
#'   with proportion \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}. 
#'   With a given model degree \code{m}, the parameters \code{p}, the mixing
#'   proportions of the beta distribution, are calculated as the minimizer of the
#'   approximate \eqn{L_2} distance between the empirical distribution and 
#'   the Bernstein polynomial model. The quadratic programming with linear constraints
#'   is used to solve the problem.
#'   If \code{search=TRUE} then the model degrees are chosen using a method of change-point based on 
#'   the marginal data if \code{mar.deg=TRUE} or the joint data if \code{mar.deg=FALSE}. 
#'   If \code{search=FALSE}, then the model degree is specified by \eqn{M}.
#' @return An invisible \code{mable} object with components
#' \itemize{
#'   \item \code{m} the given model degree(s) 
#'   \item \code{p} the estimated vector of mixture proportions 
#'       with the given optimal degree(s) \code{m}
#'   \item \code{interval} support/truncation interval \code{[a, b]}
#'   \item \code{D}  the minimum distance at degree \code{m}
#'  \item \code{convergence} An integer code. 0 indicates successful completion(the EM iteration is   
#'    convergent). 1 indicates that the iteration limit \code{maxit} had been reached in the EM iteration;
#'  }
#' @export
made.density<-function(x, M0=1L, M, search=TRUE, interval=NULL, mar.deg=TRUE, 
      method=c("qp","em"), controls = mable.ctrl(), progress=TRUE){
  data.name<-deparse(substitute(x))
  n<-NROW(x)
  d<-NCOL(x)
  xNames<-names(x)
  #method <- match.arg(method)
  if(is.null(xNames)) 
    for(i in 1:d) xNames[i]<-paste("x",i,sep='')
  x<-as.matrix(x)   
  if(is.null(interval))  
    interval<-apply(x, 2, range) 
  else if(is.matrix(interval)){
    if(nrow(interval)!=2 || ncol(interval)!=d)
      stop("Invalid 'interval'.\n")
  }
  else{
    if(length(interval)!=2) 
      stop("Invalid 'interval'.\n")
    else if(any(apply(x, 2, min)<interval[1])||any(apply(x, 2, max)>interval[2]))
      stop("Invalid 'interval'.\n")
    else interval<-matrix(rep(interval, d), nrow=2)
  }
  x<-t((t(x)-interval[1,])/apply(interval,2,diff))
  if(any(x<0) || any(x>1)) stop("All data values must be contained in 'interval'.")
  if(d==1){
    xbar<-mean(x); s2<-var(x); 
    m0<-max(1,ceiling(xbar*(1-xbar)/s2-3)-2)
    if(is.vector(x)) x<-matrix(x,ncol=d)
  }
  else{
    xbar<-colMeans(x)
    s2<-colSums(t(t(x)-xbar)^2)/(n-1)
    m0<-pmax(1,ceiling(xbar*(1-xbar)/s2-3)-2)
  }
  if(missing(M) || length(M)==0 || any(M<=0)) stop("'M' is missing or nonpositvie.\n")
  else if(length(M)<d) M<-rep(M,d)
  M<-M[1:d]
  if(min(M)<5 && search) warning("'M' are too small for searching optimal degrees.\n")
  if(length(M0)==1) M0<-rep(M0,d)
  else if(length(M0)!=d){
    if(search){
      if(min(M0)<max(M)) 
        warning("Length of 'M0' does not match 'ncol(x)'. Use least 'M0'.\n")
      else stop("'M0' are too big for searching optimal degrees.\n")
      M0<-rep(max(M0),d)}
  }
  cat("m0=",m0,"M=",M,"\n")
  convergence<-0
  Ord<-function(i){
    if(any(1:3==i%%10) && !any(11:13==i)) switch(i%%10, "st", "nd", "rd")
    else "th"}
  if(search && mar.deg && d>1){
    M0<-pmax(m0,M0)
    for(i in 1:d){
      message("mable fit of the ",i, Ord(i), " marginal data.\n")
      res<-mable(x[,i], M=c(M0[i],M[i]), c(0,1), controls=controls, progress = TRUE)
      pval<-res$pval
      message("m[",i,"]=", res$m," with p-value ",  pval[length(pval)],"\n")
      M[i]<-res$m
    }
    search<-FALSE
  }
  if(search){
    cat("M0=",M0,"M=",M,"length p=",length(p),"\n")
      level<-controls$sig.level
      Dn<-NULL
      pval<-NULL
      chpts<-NULL
      p<-NULL
      pval[1]<-1
      chpts[1]<-1
      cp<-1
      e<-diag(1,d)
      m<-M[1,]
      min.pd<-sum(m)
      max.pd<-sum(M[2,])
      cat("m=",m, "\n")
      k<-1
      pval[k]<-1
      max.k<-max.pd-min.pd
      #cat("K=",prod(m+1),"length.alpha=",length(alpha), "\n")
      tmp<-madem.density(x, m, p=rep(1,prod(m+1))/prod(m+1), interval=rbind(rep(0,d), rep(1,d)),
            method=method, maxit=controls$maxit, eps=controls$eps)
      min.D<-tmp$D
      Dn[k]<-min.D
      pmt<-tmp$p
      k<-2
      while(k<=max.k && pval[k-1]>level){
        for(i in 1:d){
          tmp<-madem.density(x, m=m+e[,i], p=p, interval=rbind(rep(0,d), rep(1,d)), 
                method=method, maxit=controls$maxit, eps=controls$eps)
          if(i==1 || tmp$D<min.D){
            min.D<-tmp$D
            Dn[k]<-min.D
            mt<-m+e[,i]
            pmt<-tmp$p
            t<-i
          }
          Dn[k]<-min(Dn[k-1],Dn[k]) #???
        }
        if(k<=4){
          pval[k]<-1
          chpts[k]<-1
          lr<-rep(0,k)
        }
        else{
          res.cp<-chpt.exp(-log(Dn))
          pval[k]<-res.cp$pvalue
          chpts[k]<-res.cp$chpt
          if(chpts[k]!=cp){
            cp<-chpts[k]
            lr<-res.cp$lr
            mhat<-mt
            phat<-pmt
          }
        }
        #cat("k=",k,"pval=",pval[k], "Dn=",Dn[k], "\n")
        m<-mt
        p<-pmt
        k<-k+1
        #cat("m=",m, "\n")
      }
      p<-round(phat,7)
      nz<-which(p!=0)
      out<-list(m=mhat,p=p/sum(p), nz=nz, M=M, D=Dn, interval=interval, lr=lr, 
                chpts=chpts, pval=pval)
      if(d>1) out$data.type<-"mvar"
      class(out)<-"mable"
  }
  else{
      out<-madem.density(x, m=M, interval=rbind(rep(0,d), rep(1,d)), 
                method=method, maxit=controls$maxit, eps=controls$eps)
  }
  
  message("\n MABLE for ", d,"-dimensional data:")
  if(search){
    message("Model degrees m = (", out$m[1], append=FALSE)
    for(i in 2:d) message(", ",out$m[i], append=FALSE)
    if(mar.deg)
      message(") are selected based on marginal data m=", out$m, "\n")
    else{
      message(") are selected between M0 and M, inclusive, where \n", append=FALSE)
      message("M0 = ", M0, "\nM = ",M,"\n")
    }
    out$M0<-M0
    out$M<-M
  }
  else{
    message("Model degrees are specified: M=", M)
  }
  invisible(out)
}



#############################################################
#' Minimum Approximate Distance Estimate
#'  of Multivariate Density Function
#' @param x an \code{n x d} matrix or \code{data.frame} of multivariate sample of size \code{n}
#' @param M0 a positive integer or a vector of \code{d} positive integers specify
#'    starting candidate degrees for searching optimal degrees.  
#' @param M a positive integer or a vector of \code{d} positive integers specify  
#'    the maximum candidate or the given model degrees for the joint density.
#' @param search logical, whether to search optimal degrees between \code{M0} and \code{M} 
#'    or not but use \code{M} as the given model degrees for the joint density.
#' @param interval a vector of two endpoints or a \code{2 x d} matrix, each column containing 
#'    the endpoints of support/truncation interval for each marginal density.
#'    If missing, the i-th column is assigned as \code{c(min(x[,i]), max(x[,i]))}.
#' @param mar.deg logical, if TRUE, the optimal degrees are selected based  
#'    on marginal data, otherwise, the optimal degrees are chosen the joint data. See details.   
#' @param method method for finding minimum distance estimate. "cd": coordinate-descent;
#     "quadprog": the quadratic programming algorithm using packages \code{quadprog} or \code{LowRankQP}.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.
#' @param progress if TRUE a text progressbar is displayed
#' @details 
#'   A \eqn{d}-variate density \eqn{f} on a hyperrectangle \eqn{[a, b]
#'   =[a_1, b_1] \times \cdots \times [a_d, b_d]} can be approximated 
#'   by a mixture of \eqn{d}-variate beta densities on \eqn{[a, b]}, 
#'   \eqn{\beta_{mj}(x) = \prod_{i=1}^d\beta_{m_i,j_i}[(x_i-a_i)/(b_i-a_i)]/(b_i-a_i)},
#'   with proportion \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}. 
#'   If \code{search=TRUE} then the model degrees are chosen using a method of change-point based on 
#'   the marginal data if \code{mar.deg=TRUE} or the joint data if \code{mar.deg=FALSE}. 
#'   If \code{search=FALSE}, then the model degree is specified by \eqn{M}.
#'   For large data and multimodal density, the search for the model degrees is 
#'   very time-consuming. In this case, it is suggested that use \code{method="cd"} 
#'   and select the degrees based on marginal data using \code{\link{mable}} or 
#'   \code{\link{optimable}}.
#' @return  A list with components
#' \itemize{
#'  \item \code{m} a vector of the selected optimal degrees by the method of change-point
#'  \item \code{p} a vector of the mixture proportions \eqn{p(j_1, \ldots, j_d)}, arranged in the 
#'   column-major order of \eqn{j = (j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}.
#'  \item \code{minD}  the minimum distance at an optimal degree \code{m}
#'  \item \code{pval}  the p-values of change-points for choosing the optimal degrees for the 
#'    marginal densities
#'  \item \code{M} the vector \code{(m1, m2, ... , md)}, where \code{mi} is the largest candidate 
#'    degree when the search stoped for the \code{i}-th marginal density
#'  \item \code{interval} support hyperrectangle \eqn{[a, b]=[a_1, b_1] \times \cdots \times [a_d, b_d]}
#'  \item \code{convergence} An integer code. 0 indicates successful completion(the EM iteration is   
#'    convergent). 1 indicates that the iteration limit \code{maxit} had been reached in the EM iteration;
#' }
#' @importFrom quadprog solve.QP
#' @importFrom LowRankQP LowRankQP
#' @examples
#' ## Old Faithful Data
#' \donttest{  
#'  library(mable)
#'  a<-c(0, 40); b<-c(7, 110)
#'  ans<- made.mvar(faithful, M = c(46,19), search =FALSE, method="quadprog", 
#'  interval = rbind(a,b), progress=FALSE)
#'  plot(ans, which="density") 
#'  plot(ans, which="cumulative")
#' }
#' @seealso \code{\link{mable}}, \code{\link{optimable}}
#' @keywords distribution nonparametric multivariate 
#' @concept multivariate Bernstein polynomial model
#' @concept density estimation
#' @author Zhong Guan <zguan@iu.edu>
#' @references 
#' Guan, Z. (2016) Efficient and robust density estimation using Bernstein type polynomials. 
#'    \emph{Journal of Nonparametric Statistics}, 28(2):250-271.
#'
#' Wang, T. and Guan, Z.,(2019) Bernstein Polynomial Model for Nonparametric Multivariate Density,    
#'    \emph{Statistics}, Vol. 53, no. 2, 321-338  
#' @export
made.mvar<-function(x, M0=1L, M, search=TRUE, interval=NULL, mar.deg=TRUE, 
      method=c("cd","quadprog"), controls = mable.ctrl(), progress=TRUE){
  data.name<-deparse(substitute(x))
  n<-NROW(x)
  d<-NCOL(x)
  xNames<-names(x)
  method <- match.arg(method)
  if(is.null(xNames)) 
    for(i in 1:d) xNames[i]<-paste("x",i,sep='')
  x<-as.matrix(x)  
  if(is.null(interval))  
    interval<-apply(x, 2, range) 
  else if(is.matrix(interval)){
    if(nrow(interval)!=2 || ncol(interval)!=d)
      stop("Invalid 'interval'.\n")
  }
  else{
    if(length(interval)!=2) 
      stop("Invalid 'interval'.\n")
    else if(any(apply(x, 2, min)<interval[1])||any(apply(x, 2, max)>interval[2]))
      stop("Invalid 'interval'.\n")
    else interval<-matrix(rep(interval, d), nrow=2)
  }
  x<-t((t(x)-interval[1,])/apply(interval,2,diff))
  if(any(x<0) || any(x>1)) stop("All data values must be contained in 'interval'.")
  if(d==1){
    xbar<-mean(x); s2<-var(x); 
    m0<-max(1,ceiling(xbar*(1-xbar)/s2-3)-2)
    if(is.vector(x)) x<-matrix(x,ncol=d)
  }
  else{
    xbar<-colMeans(x)
    s2<-colSums(t(t(x)-xbar)^2)/(n-1)
    m0<-pmax(1,ceiling(xbar*(1-xbar)/s2-3)-2)
  }
  if(missing(M) || length(M)==0 || any(M<=0)) stop("'M' is missing or nonpositvie.\n")
  else if(length(M)<d) M<-rep(M,d)
  M<-M[1:d]
  if(min(M)<5 && search) warning("'M' are too small for searching optimal degrees.\n")
  if(length(M0)==1) M0<-rep(M0,d)
  else if(length(M0)!=d){
    if(search){
      if(min(M0)<max(M)) 
        warning("Length of 'M0' does not match 'ncol(x)'. Use least 'M0'.\n")
      else stop("'M0' are too big for searching optimal degrees.\n")
      M0<-rep(max(M0),d)}
  }
  #cat("m0=",m0,"M=",M,"\n")
  convergence<-0
  Ord<-function(i){
    if(any(1:3==i%%10) && !any(11:13==i)) switch(i%%10, "st", "nd", "rd")
    else "th"}
  if(search && mar.deg && d>1){
    M0<-pmax(m0,M0)
    for(i in 1:d){
      message("mable fit of the ",i, Ord(i), " marginal data.\n")
      res<-mable(x[,i], M=c(M0[i],M[i]), c(0,1), controls=controls, progress = TRUE)
      pval<-res$pval
      message("m[",i,"]=", res$m," with p-value ",  pval[length(pval)],"\n")
      M[i]<-res$m
    }
    search<-FALSE
  }

  if(method=="cd"){
    out<-list(interval=interval, xNames=xNames)  
    Fn<-mvecdf(x) 
      if(search){          
          pd<-sum(M) # polynomial degree
          #cat("Maximum polynomial degree=",pd,"\n")
          Dn<-rep(0, pd+1)
          p<-rep(0, ceiling((pd/d+1)^d))
          #cat("M0=",M0,"M=",M,"length p=",length(p),"\n")
          pval<-rep(1, pd+1)
          lr<-rep(0, pd+1)
          chpts<-rep(0, pd+1)
          ## Call C made_cd
          res<-.C("made_cd",
            as.integer(M0), as.integer(M), as.integer(n), as.integer(d), 
            as.double(p), as.double(x), as.double(Fn), as.integer(controls$maxit), 
            as.double(controls$eps), as.double(controls$sig.level),
             as.double(pval), as.double(Dn), as.double(lr), 
            as.integer(chpts), as.logical(progress))
          m<-res[[2]]; K<-prod(m+1); 
          out$p<-res[[5]][1:K];  
          k<-res[[4]];  
          out$khat<-k
          D<-res[[12]][1:(k+1)]
          out$lr<-res[[13]][1:(k+1)]
          out$chpts<-res[[14]][1:(k+1)]
          out$minD<-D[out$chpts[k]+1]
      }
      else{
          m<-M; K<-prod(m+1); 
          p<-rep(1, K)/K
          Dn<-0
          #cat("d=",d, "n=",n, "sum p=", sum(p),  "sum Fn=", sum(Fn), "\n")
          ## Call C made_m_cd
          res<-.C("made_m_cd",
            as.integer(m), as.integer(n), as.integer(d), as.double(p), as.double(x),
            as.double(Fn), as.integer(controls$maxit), as.double(controls$eps), as.double(Dn))
          out$p<-res[[4]][1:K]; 
          D<-res[[9]][1]              
          out$minD<-D
      }
      out$m=m;
      out$D<-D
      if(d>1) out$data.type<-"mvar"
      class(out)<-"mable"
  }
  else{
      if(search)  
        out<-made.density(x, M0, M, search, interval=rbind(rep(0,d),rep(1,d)), 
                mar.deg, controls, progress)
      else
        out<-madem.density(x, M, interval=rbind(rep(0,d),rep(1,d)))
  }
  
  message("\n MABLE for ", d,"-dimensional data:")
  if(search){
    message("Model degrees m = (", out$m[1], append=FALSE)
    for(i in 2:d) message(", ",out$m[i], append=FALSE)
    if(mar.deg)
      message(") are selected based on marginal data m=", out$m, "\n")
    else{
      message(") are selected between M0 and M, inclusive, where \n", append=FALSE)
      message("M0 = ", M0, "\nM = ",M,"\n")
    }
    out$M0<-M0
    out$M<-M
  }
  else{
    message("Model degrees are specified: M=", M)
  }
  invisible(out)
}


#' Matrix of the uniform marginal constraints
#' @param m  vector of d nonnegative integers \code{m=(m[1],...,m[d])}.
#' @details the matrix of the uniform marginal constraints \eqn{A} is used to form the 
#'     linear equality constraints on parameter \eqn{p}: \eqn{Ap=1/(m+1)}.
#' @return an \code{|m| x K} matrix, \code{|m|=m[1]+...+m[d])}, \code{K=(m[1]+1)...(m[d]+1)}.
#' @export
umc.mat<-function(m){
  d<-length(m)
  sm<-c(0,cumsum(m))
  km<-c(1,cumprod(m+1))
  K <- km[d+1] 
  umc<-matrix(0, nrow=sum(m), ncol=K)
  it<-0
  while(it<K){
      r <- it 
      for(k in (d-1):1){
          jj <- r%%km[k+1] 
          i <- (r-jj)/km[k+1] 
          if(i>0) umc[sm[k+1]+i, it+1] <- 1
          r <- jj 
      }
      if(r>0) umc[r, it+1] <- 1
      it<-it+1
  }
  umc
}
#umc.mat(c(2,3))



#' Minimum Approximate Distance Estimate of Copula with given model degrees 
# 
#' @param u an \code{n x d} matrix of (pseudo) observations.
#' @param m \code{d}-vector of model degrees
#' @details With given model degrees \code{m}, the parameters \code{p}, the mixing
#'   proportions of the beta distribution, are calculated as the minimizer of the
#'   approximate \eqn{L_2} distance between the empirical distribution and 
#'   the Bernstein polynomial model. 
#' @return An invisible \code{mable} object with components
#' \itemize{
#'   \item \code{m} the given degree
#'   \item \code{p} the estimated vector of mixture proportions 
#'      \eqn{p = (p_0, \ldots, p_m)}
#'       with the given degree \code{m}
#'   \item \code{D}  the minimum distance at degree \code{m}
#'  }
#' @importFrom quadprog solve.QP
#' @importFrom LowRankQP LowRankQP
#' @export
madem.copula<-function(u, m){
#  require(quadprog)
  d<-ncol(u)
  n=nrow(u)
  eps<-.Machine$double.eps^.5
  if(any(0-u>eps) || any(u-1>eps))
  stop("all 'u' must be in [0,1].")

  y=mvecdf(u)
  ssy<-sum(y^2)
  if(all(m==0)){ # uniform 
    p=1
    Dm=mean((apply(u,1,prod)-y)^2)
  }
  else{
    sm<-sum(m)
    K=prod(m+1)
    # equality constraints
    bvec=1 # sum(p)=1
    bvec=c(bvec, rep(1/(m+1),m)) # add UMC
    Amat=rep(1, K) # sum(p)=1
    Amat=rbind(Amat, umc.mat(m)) # add UMC
    # add nonnegative constraints
    Amat=rbind(Amat, diag(1,K))
    bvec=c(bvec,rep(0, K))
    Q=diag(0, K)
    B=mvpbeta(u,m)
    for(i in 1:n) Q=Q+B[i,]%*%t(B[i,])
    cvec=t(B)%*%y
    
    #cat("dim A=",dim(Amat),"\n")
    tryCatch({res<-solve.QP(Dmat=Q, dvec=cvec, Amat=t(Amat), bvec=bvec, meq=1+sm); 
            p=res$solution; Dm=(2*res$value+ssy)/n},
        error=function(e){
          if(geterrmessage()
              =="matrix D in quadratic function is not positive definite!"){
#            suppressMessages(require(LowRankQP))
            res<-LowRankQP(Vmat=Q, dvec=-cvec, Amat=matrix(Amat[1+(0:sm),], ncol=K), 
                            bvec=bvec[1+(0:sm)], uvec=rep(1, K), method="CHOL")
            p<<-res$alpha
            Dm<<-(t(res$alpha)%*%Q%*%res$alpha-2*t(cvec)%*%res$alpha+ssy)/n
          }
          else stop(geterrmessage())
        }
    )
    #p=res$solution; D=(2*res$value+d)/n
  }
  p<-zapsmall(p, digits = 6)
  ans=list(m=m, p=p/sum(p), D=Dm, interval=matrix(rep(0:1,d), nrow=2))
  ans$xNames<-"u"
  ans$data.type<-"copula"
  ans$margin<-NULL
  class(ans)<-"mable"
  invisible(ans)
}


#' Minimum Approximate Distance Estimate of Copula Density
#' @param x an \code{n x d} matrix of data values
#' @param unif.mar marginals are all uniform (\code{x} contain pseudo observations) 
#'         or not.   
#' @param M d-vector of preselected or maximum model degrees
#' @param search logical, whether to search optimal degrees between \code{0} and \code{M} 
#'    or not but use \code{M} as the given model degrees for the joint density.
#' @param interval a 2 by d matrix specifying the support/truncate interval of \code{x}, 
#'      if \code{unif.mar=TRUE} then \code{interval} is the unit hypercube
#' @param pseudo.obs When \code{unif.mar=FALSE}, use \code{"empirical"} distribution to 
#'    create pseudo observations, or use \code{"mable"} of marginal cdfs to create 
#'    pseudo observations
#' @param sig.level significance level for p-value of change-point
#' @details With given model degrees \code{m}, the parameters \code{p}, the mixing
#'   proportions of the beta distribution, are calculated as the minimizer of the
#'   approximate \eqn{L_2} distance between the empirical distribution and 
#'   the Bernstein polynomial model. The optimal model degrees \code{m} are chosen by
#'   a change-point method. The quadratic programming with linear constraints is 
#'   used to solve the problem.
#' @return An invisible \code{mable} object with components
#' \itemize{
#'   \item \code{m} the given degree
#'   \item \code{p} the estimated vector of mixture proportions 
#'      \eqn{p = (p_0, \ldots, p_m)}
#'       with the given degree \code{m}
#'   \item \code{D}  the minimum distance at degree \code{m}
#'  }
#' @export
made.copula=function(x, unif.mar= FALSE, M=30, search=TRUE, interval=NULL,  
               pseudo.obs=c("empirical","mable"), sig.level=0.01){
  pseudo.obs<-match.arg(pseudo.obs)
  data.name<-deparse(substitute(x))
  n=nrow(x)
  d=ncol(x)
  xNames<-names(x)
  if(is.null(xNames)) 
    for(i in 1:d) xNames[i]<-paste("x",i,sep='')
  if(missing(M) || length(M)==0 || any(M<=0)) 
    stop("'M' is missing or nonpositvie.\n")
  else if(length(M)<d) M<-rep(M,d)
  M<-M[1:d]
  if(min(M)<5 && search) 
    warning("'M' are too small for searching optimal degrees.\n")
  eps<-.Machine$double.eps^.5
  if(unif.mar) {
    interval=matrix(rep(0:1,d), nrow=2)
    if(any(-x>eps) || any(x-1>eps))
      stop("Data 'x' are not all in the unit hypercube.")
    u<-x
  }
  else{
    invalid.interval<-FALSE
    if(is.matrix(interval)){
      if(any(dim(interval)!=c(2,d))) 
        invalid.interval<-TRUE
      else if(any(t(interval[1,]-t(x))>eps) || any(t(t(x)-interval[2,])>eps))
        invalid.interval<-TRUE
    }
    else{
      if(length(interval)!=2) 
        invalid.interval<-TRUE
      else if(any(t(interval[1,]-t(x))>eps) || any(t(t(x)-interval[2,])>eps))
        invalid.interval<-TRUE
      else 
        interval<-matrix(rep(interval, d), nrow=2)
    }
    if(is.null(interval) || invalid.interval){  
      interval<-apply(x, 2, extendrange)
      message("'interval is missing or invalid, set as 'extended range'.\n")
    } 
    #for(i in 1:d) {
    #  x[,i]<-(x[,i]-interval[1,i])/diff(interval[,i])
    #  interval[,i]=0:1}
    #x=(x-rep(1,n)%*%t(interval[1,]))/(rep(1,n)%*%t(interval[2,]-interval[1,]))
    margin<-list(NULL)  
    for(k in 1:d){
        margin[[k]]<-mable(x[,k], M=c(0,M[k]), interval=interval[,k], progress =FALSE)
    }
    u<-NULL
    if(pseudo.obs=="empirical")   
      for(k in 1:d) u<-cbind(u,ecdf(x[,k])(x[,k])*n/(n+1)) 
    else{
      for(k in 1:d){
        u<-cbind(u,pmixbeta(x[,k], p=margin[[k]]$p, 
                            interval=margin[[k]]$interval))
      }
    }  
  }
  if(search){
    cp<-1
    E<-diag(1,d)
    max.pd<-sum(M)
    pval<-rep(1, max.pd+1)
    chpts<-rep(0, max.pd+1)
    D<-rep(Inf, max.pd+1)
    pval[1]=pv=1.0
    mk<-mhat<-rep(0,d)
    m<-mk
    res<-madem.copula(u,mk)
    D[1]<-Dhat<-res$D
    phat<-res$p
    #nphat<-1
    k<-1
    while(k<max.pd && pv>sig.level){
      k<-k+1
      for(i in 1:d){
        mk<-m+E[i,]
        km<-c(1,cumprod(mk+1))
        #K<-km[d+1]
        res<-madem.copula(u,mk)
        if(res$D<D[k]){
          D[k]<-res$D
          mhat<-mk
          p<-res$p
        }
      }
      m<-mhat
      phat<-p
      if(k<4){
        pval[k]=1
        chpts[k]=1
      }
      else{
        res=chpt.exp(-log(D[1:k]))
        pval[k]=res$pvalue
        chpts[k]=res$chpt
        pv=pval[k]
      }
      if(chpts[k]!=cp){
        cp=chpts[k]
        mhat=m
        #nphat=K
        phat=p; 
        Dhat<-D[k]
      }
      cat("\r MADE Copula Finished:", format(round(k/max.pd,3)*100, nsmall=1), " % ")
    }
    cat("\r MADE Copula Finished:", format(100, nsmall=1), " % \n")
    cat("m=",mhat,"\n")
    phat<-zapsmall(phat, digits = 6)
    ans<-list(m=mhat, p=phat/sum(phat), D=D[1:k], Dhat=Dhat, 
              interval=matrix(rep(0:1,d), nrow=2), 
              chpts=chpts[1:k], pval=pval[1:k])
  }
  else {
    res<-madem.copula(u, m=M)
    ans<-list(m=m, p=res$p, interval=interval, D=res$D, dim=d, data.type="mvar")#,
      #M=M[,otk], Delta=Dm[otk], pval=pval[otk], cp=cp[otk])
  }
  if(!unif.mar) ans$margin<-margin
  else ans$margin<-NULL
  ans$xNames<-xNames
  ans$data.type<-"copula"
  class(ans)="mable"
  invisible(ans)
}



