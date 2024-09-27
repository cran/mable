#############################################################
## Maximum Approximate Bernstein Likelihood Estimate
##  of Multivariate Density Function
##
## References:
##  Guan, Z. (2016) Efficient and robust density estimation using Bernstein type polynomials. \emph{Journal of Nonparametric Statistics}, 28(2):250-271.
##  Wang, T. and Guan, Z.,(2019) Bernstein Polynomial Model    
##     for Nonparametric Multivariate Density, Statistics,     
##              Vol. 53, no. 2, 321-338                        
#############################################################
# setwd("C:\\Users\\zguan\\Documents\\papers\\proportional odds ratio\\parametrization using bernstein polynomials\\multivariate-bernstein-poly-model")
# setwd("C:\\Documents and Settings\\zguan\\My Documents\\papers\\proportional odds ratio\\parametrization using bernstein polynomials\\multivariate-bernstein-poly-model")
# dyn.load("mable-multivar")
## Call C functions by .C
# R --arch x64 CMD SHLIB mable-multivar.c
# R CMD SHLIB mable-multivar.c
#############################################################
#' Maximum Approximate Bernstein Likelihood Estimate
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
#' @param method method for finding maximum likelihood estimate. "cd": coordinate-descent;
#     "em": the EM like algorithm; "lmem": a slower version EM like algorithm which requires 
#'    less memory for data that are high dimensional/large sample.
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
#'  \item \code{mloglik}  the maximum log-likelihood at an optimal degree \code{m}
#'  \item \code{pval}  the p-values of change-points for choosing the optimal degrees for the 
#'    marginal densities
#'  \item \code{M} the vector \code{(m1, m2, ... , md)}, where \code{mi} is the largest candidate 
#'    degree when the search stoped for the \code{i}-th marginal density
#'  \item \code{interval} support hyperrectangle \eqn{[a, b]=[a_1, b_1] \times \cdots \times [a_d, b_d]}
#'  \item \code{convergence} An integer code. 0 indicates successful completion(the EM iteration is   
#'    convergent). 1 indicates that the iteration limit \code{maxit} had been reached in the EM iteration;
#' }
#' @examples
#' ## Old Faithful Data
#' \donttest{
#'  a<-c(0, 40); b<-c(7, 110)
#'  ans<- mable.mvar(faithful, M = c(46,19), search =FALSE, method="em",
#'          interval = rbind(a,b), progress=FALSE)
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
#' Wang, T. and Guan, Z.,(2019) Bernstein Polynomial Model for Nonparametric Multivariate Density,    
#'    \emph{Statistics}, Vol. 53, no. 2, 321-338  
#' @export
mable.mvar<-function(x, M0=1L, M, search=TRUE, interval=NULL, mar.deg=TRUE, 
      method=c("cd","em","lmem"), controls = mable.ctrl(), progress=TRUE){
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

  out<-list(interval=interval, xNames=xNames)  
  
  pd<-sum(M) # polynomial degree
  #cat("Maximum polynomial degree=",pd,"\n")
  lk<-rep(0, pd+1)
  p<-rep(0, ceiling((pd/d+1)^d))
  high.dim<-ifelse(method=="lmem", TRUE, FALSE)
  if(search){
  #cat("M0=",M0,"M=",M,"length p=",length(p),"\n")
      pval<-rep(1, pd+1)
      lr<-rep(0, pd+1)
      chpts<-rep(0, pd+1)
      d<-c(d,pd)
      ## Call C mable_cd
      if(method=="cd"){
          res<-.C("mable_cd",
            as.integer(M0), as.integer(M), as.integer(n), as.integer(d), 
            as.double(p), as.double(x), as.integer(controls$maxit), as.double(controls$eps), 
            as.double(controls$sig.level), as.double(pval), as.double(lk), as.double(lr), 
            as.integer(chpts), as.logical(progress))
      }
      else{
      ## Call C mable_mvar
          res<-.C("mable_mvar",
            as.integer(M0), as.integer(M), as.integer(n), as.integer(d), 
            as.double(p), as.double(x), as.integer(controls$maxit), as.double(controls$eps), 
            as.double(controls$sig.level), as.double(pval), as.double(lk), as.double(lr), 
        as.integer(chpts), as.logical(progress), as.integer(convergence), as.logical(high.dim))
        out$convergence=res[[15]] 
      }
      m<-res[[2]]; K<-prod(m+1); 
      out$p<-res[[5]][1:K]; out$m=m; 
      k<-res[[4]];  
      out$khat<-k
      lk<-res[[11]][1:(k+1)]
      out$lk<-lk
      out$lr<-res[[12]][1:(k+1)]
      out$chpts<-res[[13]][1:(k+1)]
      out$mloglik<-lk[out$chpts[k]+1]
  }
  else{
      m<-M; K<-prod(m+1); 
      p<-rep(1, K)/K
      ## Call C mable_m_cd
      if(method=="cd"){
          res<-.C("mable_m_cd",
            as.integer(M), as.integer(n), as.integer(d), as.double(p), as.double(x),
            as.integer(controls$maxit), as.double(controls$eps), as.double(lk))
      }
      else{
      ## Call C mable_m_mvar
          res<-.C("mable_m_mvar",
            as.integer(M), as.integer(n), as.integer(d), as.double(p), as.double(x),
            as.integer(controls$maxit), as.double(controls$eps), as.double(lk), 
            as.logical(progress), as.integer(convergence), as.logical(high.dim))
          out$convergence=res[[10]] 
      }
      out$p<-res[[4]][1:K]; out$m=m; 
      lk<-res[[8]][1]
      out$lk<-lk
      out$mloglik<-lk
  }
  
  message("\n MABLE for ", d,"-dimensional data:")
  if(search){
    message("Model degrees m = (", m[1], append=FALSE)
    for(i in 2:d) message(", ",m[i], append=FALSE)
    if(mar.deg)
      message(") are selected based on marginal data m=", m, "\n")
    else{
      message(") are selected between M0 and M, inclusive, where \n", append=FALSE)
      message("M0 = ", M0, "\nM = ",M,"\n")
    }
    out$M0<-M0
    out$M<-M
  }
  else{
    message("Model degrees are specified: M=", m)
  }
  if(d[1]>1) out$data.type<-"mvar"
  class(out)<-"mable"
  invisible(out)
}



#########################################################
#' Multivariate Mixture Beta Distribution
#' @description Density, distribution function,  and 
#' pseudorandom number generation for the multivariate Bernstein polynomial model, 
#' mixture of multivariate beta distributions, with given mixture proportions 
#' \eqn{p = (p_0, \ldots, p_{K-1})}, given degrees \eqn{m = (m_1, \ldots, m_d)},
#'  and support \code{interval}.
#' @param x a matrix with \code{d} columns or a vector of length \code{d} within 
#'   support hyperrectangle \eqn{[a, b] = [a_1, b_1] \times \cdots \times [a_d, b_d]}
#' @param p a vector of \code{K} values. All components of \code{p} must be 
#'  nonnegative and sum to one for the mixture multivariate beta distribution. See 'Details'. 
#' @param m a vector of degrees, \eqn{(m_1, \ldots, m_d)} 
#' @param n sample size
#' @param interval a vector of two endpoints or a \code{2 x d} matrix, each column containing 
#'    the endpoints of support/truncation interval for each marginal density.
#'    If missing, the i-th column is assigned as \code{c(0,1))}.
#' @details 
#'  \code{dmixmvbeta()} returns a linear combination \eqn{f_m} of \eqn{d}-variate beta densities 
#'  on \eqn{[a, b]}, \eqn{\beta_{mj}(x) = \prod_{i=1}^d\beta_{m_i,j_i}[(x_i-a_i)/(b_i-a_i)]/(b_i-a_i)},   
#'  with coefficients \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}, where
#'  \eqn{[a, b] = [a_1, b_1] \times \cdots \times [a_d, b_d]} is a hyperrectangle, and the  
#'  coefficients are arranged in the column-major order of \eqn{j = (j_1, \ldots, j_d)}, 
#'  \eqn{p_0, \ldots, p_{K-1}},  where \eqn{K = \prod_{i=1}^d (m_i+1)}. 
#'  \code{pmixmvbeta()} returns a linear combination \eqn{F_m} of the distribution
#'  functions of \eqn{d}-variate beta distribution.
#'
#'  If all \eqn{p_i}'s are nonnegative and sum to one, then \code{p}
#'  are the mixture proportions of the mixture multivariate beta distribution.
#' @importFrom stats rbeta
#' @export
dmixmvbeta<-function(x, p, m, interval=NULL){
    d<-length(m)
    if(any(m-floor(m)!=0) || d==0) stop("Invalid argument 'm'.\n")
    km<-c(1,cumprod(m+1))
    K<-length(p)
    if(K==0) stop("Missing mixture proportions 'p' without default.")
    if(K!=km[d+1]) stop("Invalid argument 'p': length must equal to 'prod(m+1)'.\n")
    if(any(p<0)) warning("Argument 'p' has negative component(s).\n")
    #else if(abs(sum(p)-1)>.Machine$double.eps)
    #    warning("Sum of 'p's is not 1.\n")
    if(is.null(interval))  
        interval<-rbind(rep(0,d), rep(1,d))
    else if(is.matrix(interval)){
        if(nrow(interval)!=2 || ncol(interval)!=d)
            stop("Invalid 'interval'.\n")
    }
    else{
        if(length(interval)!=2) 
            stop("Invalid 'interval'.\n")
        else interval<-matrix(rep(interval, d), nrow=2)
    }
    if(is.matrix(x)){
        nx<-nrow(x)
        if(d!=ncol(x)) stop("Wrong dim of 'x'.\n")
        eps<-.Machine$double.eps^.5
        #cat("interval=",interval,"\n")
        if(any(matrix(rep(interval[1,],nx), ncol=d, byrow=TRUE)-x>eps) 
            || any(x-matrix(rep(interval[2,],nx), ncol=d, byrow=TRUE)>eps))
            stop("'x' must be in the hyperrectangle 'interval'.\n")
        x<-t((t(x)-interval[1,])/apply(interval,2,diff))
    }
    else {
        if(d!=length(x)) stop("'x' is a vector. Its length must be the same as 'm'.\n") 
        nx<-1}
    fb<-rep(0, nx)
    density<-TRUE
    res<-.C("mable_mvdf",
      as.integer(d), as.integer(m), as.integer(km), as.integer(nx), as.double(x),  
      as.double(p), as.double(fb), as.logical(density))
    fb<-res[[7]]/prod(interval[2,]-interval[1,])
    return(fb)
}
#' @rdname dmixmvbeta
#' @export
pmixmvbeta<-function(x, p, m, interval=NULL){
    d<-length(m)
    if(any(m-floor(m)!=0) || d==0) stop("Invalid argument 'm'.\n")
    km<-c(1,cumprod(m+1))
    K<-length(p)
    if(K==0) stop("Missing mixture proportions 'p' without default.")
    if(K!=km[d+1]) stop("Invalid argument 'p': length must equal to prod(m+1).\n")
    if(any(p<0)) warning("Argument 'p' has negative component(s).\n")
    #else if(abs(sum(p)-1)>.Machine$double.eps)
    #    warning("Sum of 'p's is not 1.\n")
    if(is.null(interval))  
        interval<-rbind(rep(0L,d), rep(1L,d))
    else if(is.matrix(interval)){
        if(nrow(interval)!=2 || ncol(interval)!=d)
            stop("Invalid 'interval'.\n")
    }
    else{
        if(length(interval)!=2) 
            stop("Invalid 'interval'.\n")
        else interval<-matrix(rep(interval, d), nrow=2)
    }
     if(is.matrix(x)){
        nx<-nrow(x)
        if(d!=ncol(x)) stop("Wrong dim of 'x'.\n")
        eps<-.Machine$double.eps^.5
        if(any(matrix(rep(interval[1,],nx), ncol=d, byrow=TRUE)-x>eps) 
            || any(x-matrix(rep(interval[2,],nx), ncol=d, byrow=TRUE)>eps))
            stop("'x' must be in the hyperrectangle 'interval'.\n")
        x<-t((t(x)-interval[1,])/apply(interval,2,diff))
    }
    else {
        if(d!=length(x)) stop("'x' is a vector. Its length must be the same as 'm'.\n") 
        nx<-1}
#    obj<-list(m=m, p=p, interval=interval)
#    ans<-mvbern.poly(x, obj, density=FALSE)
    fb<-rep(0, nx)
    density<-FALSE
    res<-.C("mable_mvdf",
      as.integer(d), as.integer(m), as.integer(km), as.integer(nx), as.double(x),  
      as.double(p), as.double(fb), as.logical(density))
    fb<-res[[7]] 
    return(fb)
}
#########################################################
# Generating prn from mixture of multivariate beta 
#   with degrees m=(m1,...,md)
#' @rdname dmixmvbeta
#' @export
rmixmvbeta<-function(n, p, m, interval=NULL){
    d<-length(m)
    if(any(m-floor(m)!=0) || d==0) stop("Invalid argument 'm'.\n")
    K<-length(p)
    if(K==0) stop("Missing mixture proportions 'p' without default.")
    if(K!=prod(m+1)) stop("Invalid argument 'p': length must equal to prod(m+1).\n")
    if(any(p<0)) stop("Negative component(s) of argument 'p' is not allowed.\n")
    else if(abs(sum(p)-1)>.Machine$double.eps){
        warning("Sum of 'p's is not 1. Dividing 'p's by the total.\n")
        p<-p/sum(p)
    }
    if(is.null(interval))  
        interval<-rbind(rep(0L,d), rep(1L,d))
    else if(is.matrix(interval)){
        if(nrow(interval)!=2 || ncol(interval)!=d)
            stop("Invalid 'interval'.\n")
    }
    else{
        if(length(interval)!=2) 
            stop("Invalid 'interval'.\n")
        else interval<-matrix(rep(interval, d), nrow=2)
    }
    if(d==1) x<-rmixbeta(n, p, interval)
    else{
        km<-cumprod(m+1)
        # column-major order
        ii<-matrix(0, nrow=K, ncol=d)
        for(k in 1:d) ii[,k]<-rep(0:m[k], each=km[k]/(m[k]+1))
        w<-sample(1:K, n, replace = TRUE, prob = p)
        x<-rep(NULL, d)
        for(i in 1:n){
            tmp<-rbeta(d, shape1 = ii[w[i],]+1, shape2 = m-ii[w[i],]+1)
            x<-rbind(x, interval[1,]+(interval[2,]-interval[1,])*tmp)
        }
    }
    return(x)
}
#########################################################
#' The mixing proportions of marginal distribution from the mixture of 
#' multivariate beta distribution
#' @param p the mixing proportions of the mixture of multivariate beta distribution
#' @param m the model  degrees \code{m=(m1,...,md)} of the mixture of 
#'  multivariate beta distribution
#' @return a list of mixing proportions of all the marginal distributions
#' @export
marginal.p<-function(p, m){
    d<-length(m)
    mp<-list()
    if(d==1) mp[[1]]<-p
    else{
      km<-c(1,cumprod(m+1))
      K <- km[d+1] 
      for(j in 1:d) mp[[j]]<-rep(0,m[j]+1) 
      it<-0 
      while(it<K){
          r <- it 
          for(k in (d-1):1){
              jj <- r%%km[k+1] 
              i <- (r-jj)/km[k+1] 
              mp[[k+1]][i+1] <- mp[[k+1]][i+1]+p[it+1]
              r <- jj 
          }
          mp[[1]][r+1] <- mp[[1]][r+1]+p[it+1]
          it<-it+1
      }
    }
    mp
} 
