#############################################################
## Maximum Approximate Bernstein Likelihood Estimate
##  of Copula Density Function
##
## References:
##  Guan, Z.,(???) Bernstein Polynomial Model for              
##     Nonparametric Estimate of Copula Density     
#############################################################
# setwd("C:\\Users\\zguan\\Documents\\papers\\proportional odds ratio\\parametrization using bernstein polynomials\\multivariate-bernstein-poly-model")
# setwd("C:\\Documents and Settings\\zguan\\My Documents\\papers\\proportional odds ratio\\parametrization using bernstein polynomials\\multivariate-bernstein-poly-model")
# dyn.load("mable-multivar")
## Call C functions by .C
# R --arch x64 CMD SHLIB mable-multivar.c
# R CMD SHLIB mable-multivar.c
#############################################################
#' Maximum Approximate Bernstein Likelihood Estimate
#'  of Copula Density Function
#' @param x an \code{n x d} matrix or \code{data.frame} of multivariate sample of size \code{n}
#'    from d-variate distribution with hyperrectangular specified by \code{interval}. 
#' @param M0 a nonnegative integer or a vector of \code{d} nonnegative integers specify
#'    starting candidate degrees for searching optimal degrees.  
#' @param M a positive integer or a vector of \code{d} positive integers specify  
#'    the maximum candidate or the given model degrees for the joint density.
#' @param unif.mar logical, whether all the marginals distributions are uniform or not. 
#'    If not the pseudo observations will be created using \code{empirical} or \code{mable}  
#'    marginal distributions. 
#' @param pseudo.obs \code{"empirical"}: use empirical distribution to create pseudo,
#'    observations, or \code{"mable"}: use mable of marginal cdfs to create pseudo observations
#' @param interval a vector of two endpoints or a \code{2 x d} matrix, each column containing 
#'    the endpoints of support/truncation interval for each marginal density.
#'    If missing, the i-th column is assigned as \code{extendrange(x[,i])}.
# \code{c(min(x[,i]), max(x[,i]))}.
#'    If \code{unif.mar=TRUE}, then it is \eqn{[0,1]^d}.
#' @param search logical, whether to search optimal degrees between \code{M0} and \code{M} 
#'    or not but use \code{M} as the given model degrees for the joint density.
#' @param mar.deg logical, if TRUE (default), the optimal degrees are selected based on marginal data,   
#'    otherwise, the optimal degrees are chosen by the method of change-point. See details.   
#' @param high.dim logical, data are high dimensional/large sample or not
#'    if TRUE, run a slower version procedure which requires less memory  
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.
#' @param progress if TRUE a text progressbar is displayed
#' @details 
#'   A \eqn{d}-variate copula density \eqn{c(u)} on \eqn{[0, 1]^d} can be approximated 
#'   by a mixture of \eqn{d}-variate beta densities on \eqn{[0, 1]^d}, 
#'   \eqn{\beta_{mj}(x) = \prod_{i=1}^d\beta_{m_i,j_i}(u_i)},
#'   with proportion \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}, 
#'   which satisfy the uniform marginal constraints, the copula (density) has   
#'   uniform marginal cdf (pdf). If \code{search=TRUE} and \code{mar.deg=TRUE}, then the 
#'   optimal degrees are \eqn{(\tilde m_1,\ldots,\tilde m_d)}, where \eqn{\tilde m_i} is 
#'   chosen based on marginal data of \eqn{u_i}, $\eqn{i=1,\ldots,d}. If \code{search=TRUE} 
#'   and \code{mar.deg=FALSE}, then the optimal degrees \eqn{(\hat m_1,\ldots,\hat m_d)}
#'   are chosen using a change-point method based on the joint data. 
#'
#'   For large data and high dimensional density, the search for the model degrees might be 
#'   time-consuming. Thus patience is needed.  
#' @return  A list with components
#' \itemize{
#'  \item \code{m} a vector of the selected optimal degrees by the method of change-point
#'  \item \code{p} a vector of the mixture proportions \eqn{p(j_1, \ldots, j_d)}, arranged in the 
#'   column-major order of \eqn{j = (j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}.
#'  \item \code{mloglik}  the maximum log-likelihood at an optimal degree \code{m}
#'  \item \code{pval}  the p-values of change-points for choosing the optimal degrees for the 
#'    marginal densities
#'  \item \code{M} the vector \code{(m1, m2, ..., md)} at which the search of model degrees stopped.
#'    If \code{mar.deg=TRUE} \code{mi} is the largest candidate degree when the search stoped for 
#'    the \code{i}-th marginal density
#'  \item \code{convergence} An integer code. 0 indicates successful completion(the EM iteration is   
#'    convergent). 1 indicates that the iteration limit \code{maxit} had been reached in the EM iteration;
#'  \item if \code{unif.mar=FALSE}, \code{margin} contains objects of the results of mable fit  
#'       to the marginal data
#' }
#' @examples
#' ## Simulated bivariate data from Gaussian copula
#' \donttest{ 
#'  set.seed(1)
#'  rho<-0.4; n<-1000
#'  x<-rnorm(n)
#'  u<-pnorm(cbind(rnorm(n, mean=rho*x, sd=sqrt(1-rho^2)),x))
#'  res<- mable.copula(u, M = c(3,3), search =FALSE, mar.deg=FALSE,  progress=FALSE)
#'  plot(res, which="density") 
#' }
#' @seealso \code{\link{mable}}, \code{\link{mable.mvar}}
#' @keywords copula distribution nonparametric multivariate 
#' @concept multivariate Bernstein polynomial model for copula
#' @concept copula density estimation
#' @author Zhong Guan <zguan@iu.edu>
#' @references 
#'   Wang, T. and Guan, Z. (2019). Bernstein polynomial model for nonparametric 
#'       multivariate density. Statistics 53(2), 321–338.
#'   Guan, Z.,  Nonparametric Maximum Likelihood Estimation of Copula              
#' @importFrom stats ecdf
#' @importFrom grDevices extendrange 
#' @export
mable.copula<-function(x, M0=1, M, unif.mar=TRUE, pseudo.obs=c("empirical","mable"), 
      interval=NULL, search=TRUE, mar.deg=FALSE, high.dim=FALSE, 
      controls = mable.ctrl(sig.level=0.05), progress=TRUE){
  pseudo.obs<-match.arg(pseudo.obs)
  data.name<-deparse(substitute(x))
  n<-nrow(x)
  d<-ncol(x)
  xNames<-names(x)
  if(is.null(xNames)) 
    for(i in 1:d) xNames[i]<-paste("x",i,sep='')
  #x<-as.matrix(x) 
  if(missing(M) || length(M)==0 || any(M<=0)) stop("'M' is missing or nonpositvie.\n")
  else if(length(M)<d) M<-rep(M,d)
  M<-M[1:d]
  if(min(M)<5 && search) warning("'M' are too small for searching optimal degrees.\n")
  if(length(M0)==1) M0<-rep(M0,d)
  else if(length(M0)!=d){
    if(search){
      if(min(M0)<max(M)) 
        warning("Length of 'M0' does not match 'ncol(u)'. Use least 'M0'.\n")
      else stop("'M0' is/are too big for searching optimal degrees.\n")
      M0<-rep(min(M0),d)}
  }
  #cat("M=(",M[1],",",M[2],")\n")
  if(unif.mar){  
    interval<-matrix(rep(0:1, d), nrow=2) 
    u<-x
  }
  else{
    invalid.interval<-FALSE
    if(is.matrix(interval)){
      if(any(dim(interval)!=c(2,d))) invalid.interval<-TRUE
      else if(any(matrix(rep(interval[1,],n), ncol=d, byrow=TRUE)-x>.Machine$double.eps) 
            || any(x-matrix(rep(interval[2,],n), ncol=d, byrow=TRUE)>.Machine$double.eps))
        invalid.interval<-TRUE
    }
    else{
      if(length(interval)!=2) invalid.interval<-TRUE
      else if(any(apply(x, 2, min)<interval[1])||any(apply(x, 2, max)>interval[2]))
        invalid.interval<-TRUE
      else interval<-matrix(rep(interval, d), nrow=2)
    }
    if(is.null(interval) || invalid.interval){  
      interval<-apply(x, 2, extendrange)
      message("'interval is missing or invalid, set as 'extended range'.\n")
    } 
    u<-NULL
    if(pseudo.obs=="empirical")   
      for(k in 1:d) u<-cbind(u,ecdf(x[,k])(x[,k])*n/(n+1)) 
    if(pseudo.obs=="mable" || mar.deg){
      Ord<-function(i){
        if(any(1:3==i%%10) && !any(11:13==i)) switch(i%%10, "st", "nd", "rd")
        else "th"}
      margin<-list(NULL)  
      #cat("M0=(",M0[1],",",M0[2],")\n")
      for(k in 1:d){
        message("mable fit of the ",k, Ord(k), " marginal data.\n")
        #cat("interval[,",k,"]=(",interval[1,k],",",interval[2,k],")\n")
        margin[[k]]<-mable(x[,k], M=c(M0[k],M[k]), interval=interval[,k], progress =FALSE)
        #margin[[k]]<-mres
        #cat("sum(p)=", sum(mres$p),"\n")
      }
      if(mar.deg){
          for(k in 1:d){
            pval<-margin[[k]]$pval
            M[k]<-margin[[k]]$m
            message("m[",k,"]=", M[k]," with p-value ",  pval[length(pval)],"\n", appendLF=FALSE)
          }
          message("Optimal degrees m = (", M[1], appendLF=FALSE)
          for(i in 2:d) message(", ",M[i],appendLF=FALSE)
          message(") are selected based on marginal data.\n")
          #p<-rep(0, prod(M+1))
          search=FALSE
      }
      if(pseudo.obs=="mable"){
        for(k in 1:d){
            #cat("sum(p)=", sum(margin[[k]]$p),"\n")
            u<-cbind(u,pmixbeta(x[,k], p=margin[[k]]$p, interval=margin[[k]]$interval))
        }
      }
    }  
  }
  m<-rep(0, d)
  pl<-list()
  mlik<-0    
  convergence<-0
  #if(!search) message("Model degrees are selected or prespecified: M=", M)
  pd<-sum(M) # polynomial degree
  #cat("Maximum polynomial degree=",pd,"\n")
  pval<-rep(1, pd+1)
  lk<-rep(-1.0e+10, pd+1)
  lr<-rep(0, pd+1)
  chpts<-rep(0, pd+1)
  p<-rep(0, prod(M+1))
  ## Call C mable_copula
  res<-.C("mable_copula",
    as.integer(M), as.integer(n), as.integer(d), as.integer(search), as.double(p), 
    as.integer(m), as.double(u), as.integer(controls$maxit), as.double(controls$eps), 
    as.double(controls$sig.level), as.double(pval), as.double(lk), as.double(lr), 
    as.integer(chpts), as.logical(progress), as.integer(convergence), as.logical(high.dim))
  m<-res[[6]]
  K<-prod(m+1)
  out<-list(p=res[[5]][1:K], m=m, xNames=xNames, interval=matrix(rep(0:1,d),2), convergence=res[[16]])
  #cat("mable.copula:", out$p,"\n")
  #if(!unif.mar) out$margin<-margin
  #else out$margin<-NULL
  if(search){
      k<-res[[3]]
      lk<-res[[12]][1:(k+1)]
      out$lk<-lk
      out$lr<-res[[13]][1:(k+1)]
      out$chpts<-res[[14]][1:(k+1)]
      out$mloglik<-lk[out$chpts[k]+1]
      out$khat<-k
  }
#  message("\n MABLE of Copula for ", d,"-dimensional data:")
  if(search){
    message("\n Optimal degrees m = (", m[1], appendLF=FALSE)
    for(i in 2:d) message(", ",m[i], appendLF=FALSE)
      message(") are selected by change-point method. \n")
  }
  
  out$M<-M
  out$data.type<-"copula"
  class(out)<-"mable"
  invisible(out)
}

#########################################################
#' Bhattacharyya coefficient and Hellinger correlation
#' @param dcopula a function object defining a 2d copula density function
#' @param ... argument(s) of copula density function
#' @return Bhattacharyya coefficient \code{B} and Hellinger correlation \code{eta}
#' @references Geenens, G. and Lafaye de Micheaux, P. (2022). The Hellinger correlation. 
#'     Journal of the American Statistical Association 117(538), 639–653.             
#' @export
corr.hellinger<-function(dcopula, ...){
  Bhat<-integrate(function(y) {
    sapply(y, function(y) {
      integrate(function(x) sqrt(dcopula(x,y, ...)), 0, 1)$value
    })
  }, 0, 1)$value
  eta<-2*sqrt(Bhat^4+sqrt(4-3*Bhat^4)-2)/Bhat^2
  list(Bhat=Bhat, eta=eta)
}

####################################################################################
#' Estimate of Hellinger Correlation between two random variables and Bootstrap 
#' @param x an \code{n x 2} data matrix of observations of the two random variables
#' @param integral logical, using "integrate()" or not (Riemann sum)
#' @param interval a 2 by 2 matrix, columns are the marginal supports 
#' @param B the number of bootstrap samples and number of Monte Carlo runs for
#'    estimating \code{p.value} of the test for Hellinger correlation = 0 
#'    if \code{test=TRUE}.
#' @param conf.level confidence level
#' @param unif.mar logical, whether all the marginals distributions are uniform or not. 
#'    If not the pseudo observations will be created using \code{empirical} or \code{mable}  
#'    marginal distributions. 
#' @param pseudo.obs \code{"empirical"}: use empirical distribution to form pseudo, 
#'    observations, or \code{"mable"}: use mable of marginal cdfs to form pseudo 
#'    observations
#' @param M0 a nonnegative integer or a vector of \code{d} nonnegative integers specify
#'    starting candidate degrees for searching optimal degrees.  
#' @param M a positive integer or a vector of \code{d} positive integers specify  
#'    the maximum candidate or the given model degrees for the joint density.
#' @param search logical, whether to search optimal degrees between \code{M0} and \code{M} 
#'    or not but use \code{M} as the given model degrees for the joint density.
#' @param mar.deg logical, if TRUE (default), the optimal degrees are selected based on marginal data,   
#'    otherwise, the optimal degrees are chosen by the method of change-point. See details.   
#' @param high.dim logical, data are high dimensional/large sample or not
#'    if TRUE, run a slower version procedure which requires less memory  
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.
#' @param progress if TRUE a text progressbar is displayed
#' @details  This function calls \code{mable.copula()} for estimation of the copula density.    
#' @return 
#'  \itemize{
#'   \item \code{eta} Hellinger correlation
#'   \item \code{CI.eta} Bootstrap confidence interval for 
#'          Hellinger correlation if \code{B}>0.
#'  }
#' @seealso \code{\link{mable}}, \code{\link{mable.mvar}}, \code{\link{mable.copula}}
#' @keywords copula distribution nonparametric multivariate 
#' @concept multivariate Bernstein polynomial model for copula
#' @concept copula density estimation
#' @author Zhong Guan <zguan@iu.edu>
#' @references Guan, Z.,  Nonparametric Maximum Likelihood Estimation of Copula              
#' @importFrom stats integrate quantile sd 
#' @export
mable.hellcorr<-function(x, unif.mar=FALSE, pseudo.obs=c("empirical","mable"), 
     M0=c(1,1), M=c(30,30), search=TRUE, mar.deg=TRUE, high.dim=FALSE, 
     interval=cbind(0:1,0:1), B=200L, conf.level=0.95, integral=TRUE, 
     controls = mable.ctrl(sig.level=0.05), progress=FALSE){
#  bsource<-match.arg(bsource)
  Res<-mable.copula(x, M0=M0, M=M, unif.mar, pseudo.obs, interval, search, 
            mar.deg, high.dim, controls, progress)
  #cat("hellcor:", Res$p,"\n")
  #res<-Res$marg.fit
  # calculate Hellinger Correlation
  eta.hat<-function(res, integral){
    fn<-function(x1, x2) dmixmvbeta(cbind(x1, x2), res$p, res$m, interval=cbind(0:1,0:1))
    if(integral) ans<-corr.hellinger(fn)
    else {
      N<-513
      x1<-seq(0,1,length=N)
      x2 <- x1
      out <- outer(x1, x2, fn)
      Bhat<-mean(sqrt(out)) # Riemann sum
      eta<-2*sqrt(Bhat^4+sqrt(4-3*Bhat^4)-2)/Bhat^2
      ans<-list(Bhat=Bhat, eta=eta)
    }
    ans
  }
  ans<-eta.hat(Res, integral)
  if(B==0) return(ans)
  else {
    eta.bt<-NULL
    n<-nrow(x)
    for(b in 1:B){
      bu<-rmixmvbeta(n, Res$p, Res$m, interval=cbind(0:1,0:1))
      bres<-suppressWarnings(mable.copula(bu, M0=1, M=Res$m, unif.mar=TRUE,
               pseudo.obs, interval=cbind(0:1,0:1), search=FALSE, mar.deg, 
               high.dim, controls, progress=FALSE))
      eta.bt[b]<-eta.hat(bres, integral)$eta
      cat("\r Bootstrap Finished:", format(round(b/B,3)*100, nsmall=1), " % ")
    }

#    else{
#      fxy<-suppressWarnings(mable.mvar(x, M0=M0, M=M, search=search, 
#                 mar.deg=mar.deg, interval=interval, progress=FALSE))
#      for(b in 1:B){
#        bx<-rmixmvbeta(n, fxy$p, fxy$m, interval=interval)
#        bu<-NULL
#        if(pseudo.obs=="empirical")
#          for(k in 1:2){
#            bu<-cbind(bu,ecdf(bx[,k])(bx[,k]))
#          }
#        else
#          for(k in 1:2){
#            bres<-suppressWarnings(mable(bx[,k], M=fxy$m, 
#                        interval=fxy$interval[,k], progress=FALSE))
#            bu<-cbind(bu,pmixbeta(bx[,k], bres$p, interval=fxy$interval[,k]))
#          }
#        Bres<-suppressWarnings(mable.copula(bu, M0=1, M=M, search, mar.deg,  
#                  high.dim, controls, progress=FALSE))
#        eta.bt[b]<-eta.hat(Bres, integral)$eta
#        cat("\r Bootstrap Finished:", format(round(b/B,3)*100, nsmall=1), " % ")
#      }
#    }
    ans$CI.eta<-2*ans$eta-quantile(eta.bt, (1+c(1,-1)*conf.level)/2)
    ans$se<-sd(eta.bt)
    cat("\r Bootstrap Finished:", format(round(b/B,3)*100, nsmall=1), " % \n")
    ans$CI.eta<-pmin(pmax(ans$CI.eta,c(0,0)),c(1,1))
    names(ans$CI.eta)<-names(ans$CI.eta)[c(2,1)]
    return(ans)
  }
}


########################################################
# alias of mable.hellcorr
#' @rdname mable.hellcorr
#' @export
hellcorr<-mable.hellcorr
