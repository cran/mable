
##########################################################
# Maximum Approximate Bernstein Likelihood Estimation
#  for density ratio model f1(x)=f0(x)exp[alpha0+alpha'r(x)]
#          with r(x)=(r1(x),...,r_d(x))
# If support is (a,b) then replace r(x) by r[a+(b-a)x]
##########################################################
#setwd("C:/Users/zguan/Documents/papers/bernstein polynomials/density ratio/C")
#dyn.load("mable-density-ratio-model")
#
# R CMD SHLIB mable-density-ratio-model.c
# R --arch x64 CMD SHLIB mable-dr-model.c
# dyn.load("mable-dr-model")
#
##################################################################
#' MABLE in Desnity Ratio Model
#' @description  Maximum approximate Bernstein/Beta likelihood estimation in a
#'   density ratio model based on two-sample raw data. 
#' @param x,y  original two sample raw data, \code{x}:"Control", \code{y}: "Case".
#' @param M a positive integer or a vector \code{(m0, m1)}. 
#' @param regr regressor vector function \eqn{r(x)=(1,r_1(x),...,r_d(x))} 
#'     which returns n x (d+1) matrix, n=length(x)
#' @param ... additional arguments to be passed to regr
#' @param interval a vector \code{(a,b)} containing the endpoints of 
#'   supporting/truncation interval of x and y.
#' @param alpha initial regression coefficient, missing value is imputed by 
#'     logistic regression
#' @param vb code for vanishing boundary constraints, -1: f0(a)=0 only, 
#'     1: f0(b)=0 only, 2: both, 0: none (default).
#' @param baseline the working baseline, "Control" or "Case", if \code{NULL}  
#'     it is chosen to the one with smaller estimated lower bound for model degree.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#'    and the convergence criterion for EM and Newton iterations. Default is 
#'    \code{\link{mable.ctrl}}. See Details.
#' @param progress logical: should a text progressbar be displayed
#' @param message logical: should warning messages be displayed
#' @details 
#' Suppose that \code{x} ("control") and \code{y} ("case") are independent  
#' samples from f0 and f1 which  samples
#' satisfy f1(x)=f0(x)exp[alpha0+alpha'r(x)] with r(x)=(r1(x),...,r_d(x)). Maximum 
#' approximate Bernstein/Beta likelihood estimates of (alpha0,alpha), f0 and f1 
#' are calculated. If support is (a,b) then replace r(x) by r[a+(b-a)x].
#' For a fixed \code{m}, using the Bernstein polynomial model for baseline \eqn{f_0},
#' MABLEs of \eqn{f_0} and parameters alpha can be estimated by EM algorithm and Newton  
#' iteration. If estimated lower bound \eqn{m_b} for \code{m} based on \code{y}
#' is smaller that that based on \code{x}, then switch \code{x} and \code{y} and
#' \eqn{f_1} is used as baseline. If \code{M=m} or \code{m0=m1=m}, then \code{m} is a 
#' preselected degree. If \code{m0<m1} it specifies the set of consective  
#' candidate model degrees \code{m0:m1} for searching an optimal degree by 
#' the change-point method, where \code{m1-m0>3}.  
#' @return  A list with components
#' \itemize{
#'   \item \code{m} the given or a selected degree by method of change-point
#'   \item \code{p} the estimated vector of mixture proportions \eqn{p = (p_0, \ldots, p_m)}
#'       with the given or selected degree \code{m}
#'   \item \code{alpha} the estimated regression coefficients
#'   \item \code{mloglik}  the maximum log-likelihood at degree \code{m}
#'   \item \code{interval} support/truncation interval \code{(a,b)}
#'   \item \code{baseline} ="control" if \eqn{f_0} is used as baseline, 
#'     or ="case" if \eqn{f_1} is used as baseline.
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1}, if greater than \code{m0}, is the
#'      largest candidate when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references Guan, Z., Maximum Approximate Bernstein Likelihood Estimation of 
#'   Densities in a Two-sample Semiparametric Model
#' @examples
#' \donttest{
#' # Hosmer and Lemeshow (1989): 
#' # ages and the status of coronary disease (CHD) of 100 subjects 
#' x<-c(20, 23, 24, 25, 26, 26, 28, 28, 29, 30, 30, 30, 30, 30, 32,
#' 32, 33, 33, 34, 34, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39,
#' 40, 41, 41, 42, 42, 42, 43, 43, 44, 44, 45, 46, 47, 47, 48, 49,
#' 49, 50, 51, 52, 55, 57, 57, 58, 60, 64)
#' y<-c(25, 30, 34, 36, 37, 39, 40, 42, 43, 44, 44, 45, 46, 47, 48,
#' 48, 49, 50, 52, 53, 53, 54, 55, 55, 56, 56, 56, 57, 57, 57, 57,
#' 58, 58, 59, 59, 60, 61, 62, 62, 63, 64, 65, 69)
#' regr<-function(x) cbind(1,x)
#' chd.mable<-mable.dr(x, y, M=c(1, 15), regr, interval = c(20, 70))
#' chd.mable
#' }
#' @importFrom stats binomial glm 
#' @export

mable.dr<-function(x, y, M, regr, ..., interval = c(0,1), alpha=NULL, vb=0, 
      baseline=NULL, controls = mable.ctrl(),  progress=TRUE, message=FALSE){
    a<-interval[1]
    b<-interval[2]
    eps.em=controls$eps.em
    eps.nt=controls$eps.nt
    maxit.em=controls$maxit.em 
    maxit.nt=controls$maxit.nt
    tini=controls$tini
    level=controls$sig.level
    if(a>=b) stop("'a' must be smaller than 'b'")
    nx<-length(x)
    ny<-length(y)
    if(ny==0 || ny==1) stop("'case' data 'y' are missing or too small.\n")
    x1<-(x-a)/(b-a)
    y1<-(y-a)/(b-a)
    regr <- match.fun(regr)
    ff <- function(x) regr(a+(b-a)*x,...)
    ry<-ff(y1)
    if(is.null(alpha)){
      if(nx==0)stop("Missing 'control' data 'x', 'alpha' must be given.")
      z<-c(rep(0,nx),rep(1,ny))
      alpha=glm(z~regr(c(x,y))[,-1], family= binomial)$coefficients
      alpha[1]=alpha[1]+log(nx/ny)
    }
    d<-length(alpha)-1
    #
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
    if(!is.null(baseline)){
      if(baseline!="Control" && baseline!="Case") stop("Invalid 'baseline'.\n")
      m0<-M[1]}
    else{
      ybar<-mean(y1); sy2<-var(y1); 
      my<-max(0,ceiling(ybar*(1-ybar)/sy2-3)-2)
      if(nx>0){
        xbar<-mean(x1); sx2<-var(x1); 
        mx<-max(0,ceiling(xbar*(1-xbar)/sx2-3)-2)
        m0<-min(mx,my)
      }else{
        m0<-my
        mx<--1}
      if(mx<=my) baseline<-"Control"
      else baseline<-"Case"}
    if(M[1]<m0){
      message("Replace M[1]=",M[1]," by the recommended ", m0,".")
      M[1]=m0; M[2]=max(m0,M[2])}
    k<-M[2]-M[1]
    if(k>0 && k<=3){
      message("Too few candidate model degrees. Add 5 more.")
      M[2]<-M[2]+5}
    #p<-rep(1,M[1]+1)/(M[1]+1)
    level<-controls$sig.level
    ## Call C_mable_dr
    if(baseline=="Control"){
      message("Use 'control' as baseline.")
      wk<-.External("C_mable_dr", ff, rho = environment(),
        as.double(x1), as.double(y1), as.double(ry), as.integer(M), 
        as.double(alpha), as.integer(d), as.integer(nx), as.integer(ny), 
        as.double(eps.em), as.double(eps.nt), as.integer(maxit.em),
        as.integer(maxit.nt),as.double(tini), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
    }else{
      message("Use 'case' as baseline.")
      wk<-.External("C_mable_dr", ff, rho = environment(),
        as.double(y1), as.double(x1), as.double(ff(x1)), as.integer(M), 
        as.double(-alpha), as.integer(d), as.integer(ny), as.integer(nx), 
        as.double(eps.em), as.double(eps.nt), as.integer(maxit.em),
        as.integer(maxit.nt),as.double(tini), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
        #wk$alpha<- -wk$alpha
    }
    res <- c(list(regr=ff), list(x=x, y=y, interval=interval, 
            mloglik=wk$lk[wk$m-M[1]+1], baseline=baseline),
            wk[c("lk", "lr", "p", "wt", "m", "alpha","chpts","pval","M")])
    
    class(res) <- "mable_dr"
    res
}

 
########################################################
#  Bootstrap sd of coefficients in density ratio model #
########################################################
# Internal??

#' Standard errors of coefficients in density ratio model
#' @description  Bootstrap estimates of standard errors for the regression
#'   coefficients which are estimated by maximum approximate Bernstein/Beta 
#'   likelihood estimation method in a density ratio model based on two-sample 
#'   raw data. 
#' @param obj  Class 'mable_dr' object return by \code{mable.dr} or \code{mable.dr.group}  functions
#' @param grouped logical: are data grouped or not.
#' @param B number of bootstrap runs. 
#' @param parallel logical: do parallel or not. 
#' @param ncore number of cores used for parallel computing. Default is half of availables. 
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#'    and the convergence criterion for EM and Newton iterations. Default is 
#'    \code{\link{mable.ctrl}}. See Details.
#' @details Bootstrap method is used based on bootstrap samples generated from
#' the MABLE's of the densities f0 and f1. The bootstrap samples are fitted by
#' the Bernstein polynomial model and the \code{glm()} to obtain bootstrap 
#' versions of coefficient estimates.
#' @return the estimated standard errors 
#' @importFrom stats binomial glm var pnorm
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom foreach foreach %dopar% 
#' @importFrom iterators icount
#' @importFrom tcltk tkProgressBar setTkProgressBar
#' @importFrom doParallel registerDoParallel
# @importFrom pbapply pboptions pblapply 
#' @export


se.coef.dr<-function(obj, grouped=FALSE, B=500L, parallel=FALSE, ncore=NULL,
      controls = mable.ctrl()){
  alpha<-obj$alpha; p<-obj$p
  d<-length(alpha)-1 
  m<-obj$m
  interval<-obj$interval
  a<-interval[1]
  b<-interval[2]    
  regr<-function(x) obj$regr((x-a)/(b-a))
  if(grouped){
    t<-obj$t; n0<-obj$n0; n1<-obj$n1; N<-length(t)-1
    nx<-sum(n0); ny<-sum(n1)
    x<-rep.int((t[-N]+t[-1])/2,times=n0)
    y<-rep.int((t[-N]+t[-1])/2,times=n1) 
  }
  else{
    x<-obj$x; y<-obj$y
    nx<-length(x); ny<-length(y)
    n<-nx+ny}
  alfa.hat<-matrix(0,nrow=d+1,ncol=B)
  delta<-c(rep(0,nx),rep(1,ny)) 
  do.bt<-function(progress=FALSE){
    xx<-a+(b-a)*rmixbeta(n=nx, p=p)
    yy<-rtmixbeta(n=ny, p=p, alpha=alpha, interval, regr=regr)
    if(obj$baseline=="Case"){
      tmp<-xx; xx<-yy; yy<-tmp}
#   cat("Bootstrap: b=",r,"\n")
    if(grouped){
      nn0<-NULL->nn1
      for(i in 1:N){
        nn0[i]<-sum(xx>=t[i] & xx<t[i+1])
        nn1[i]<-sum(yy>=t[i] & yy<t[i+1])}  
      res.b<-suppressMessages(mable.dr.group(t, nn0, nn1, M=m, regr, interval=interval, 
        alpha=alpha, controls=controls, progress=progress, message=FALSE))
          
    }else{   
      res.b<-suppressMessages(mable.dr(x=xx, y=yy, M=m, regr, interval=interval, alpha=alpha, 
        controls=controls, progress=progress, message=FALSE))
    }
    alfa.hat<-res.b$alpha
    if(res.b$baseline=="Case") alfa.hat<--alfa.hat
    alfa.hat
  }
  if(parallel){
    if(is.null(ncore)) ncore<-detectCores()
    cl<-makeCluster(ncore/2)
    registerDoParallel(cl)
    i=0
    #if(txtbar){
    #  pb <- txtProgressBar(min = 1, max = B, style = 3)
    #  THETA <- foreach(i=icount(B), .combine=rbind) %dopar% {
    #    res<-try(do.bt(), TRUE)
    #    setTxtProgressBar(pb, i)
    #    return(res)
    #  }
    #  close(pb)
    #}else{
      THETA <- foreach(i=icount(B), .packages = "tcltk", .combine=rbind) %dopar% {
        res<-try(do.bt(), TRUE)
        if(!exists("pb")) pb <- tkProgressBar("Parallel task", min=1, max=B)
        setTkProgressBar(pb, i)
        #Sys.sleep(0.05)
        return(res)
      }
    #}
    THETA<-as.matrix(THETA)
    stopCluster(cl)
  }else{
  if(requireNamespace("pbapply", quietly = TRUE)){
    op <- pbapply::pboptions(type = "win")
    Theta<- pbapply::pblapply(1:B, function(i) try(do.bt(), TRUE))
    THETA<- unlist(Theta[sapply(Theta, function(x) !inherits(x, "try-error"))])
    pbapply::pboptions(op)
  }else{
    Theta<- lapply(1:B, function(i) try(do.bt(), TRUE))
    THETA<- unlist(Theta[sapply(Theta, function(x) !inherits(x, "try-error"))])
  }
  THETA<-matrix(THETA, nrow=length(THETA)/(d+1), byrow=TRUE)
  }
  if(obj$baseline=="Case") alpha <- -alpha
  Bias<-apply(THETA, 2, mean)-alpha
  Se<-sqrt(apply(THETA, 2, var)+Bias^2)
  Zval<-alpha/Se
  Pval<-2*pnorm(-abs(Zval))
  mable.fit<-cbind(alpha, Se, Zval, Pval)
  colnames(mable.fit)<-c("Estimate","Std. Error","z value", "Pr(>|z|)")
  rownames(mable.fit)<-paste("alpha", 0:d, sep='')
  glm.fit=summary(glm(delta~regr(c(x, y))[,-1], family=binomial))$coef
  glm.fit[1,1]<-glm.fit[1,1]+log(nx/ny)
  list(mable.fit=mable.fit, glm.fit=glm.fit)
}

#' Maximum approximate profile likelihood estimate of the density ratio model 
#' @description Select optimal degree with a given regression coefficients. 
#' @param x,y  original two sample raw data, \code{x}:"Control", \code{y}: "Case".
#' @param M a positive integer or a vector \code{(m0, m1)}. 
#' @param regr regressor vector function \eqn{r(x)=(1,r_1(x),...,r_d(x))} 
#'     which returns n x (d+1) matrix, n=length(x)
#' @param ... additional arguments to be passed to regr
#' @param interval a vector \code{(a,b)} containing the endpoints of 
#'   supporting/truncation interval of x and y.
#' @param alpha a given regression coefficient, missing value is imputed by 
#'     logistic regression
#' @param vb code for vanishing boundary constraints, -1: f0(a)=0 only, 
#'     1: f0(b)=0 only, 2: both, 0: none (default).
#' @param baseline the working baseline, "Control" or "Case", if \code{NULL}  
#'     it is chosen to the one with smaller estimated lower bound for model degree.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#'    and the convergence criterion for EM and Newton iterations. Default is 
#'    \code{\link{mable.ctrl}}. See Details.
#' @param progress logical: should a text progressbar be displayed
#' @param message logical: should warning messages be displayed
#' @details 
#' Suppose that ("control") and \code{y} ("case") are independent samples from   
#' f0 and f1 which satisfy f1(x)=f0(x)exp[alpha0+alpha'r(x)] 
#' with r(x)=(r1(x),...,r_d(x)). Maximum 
#' approximate Bernstein/Beta likelihood estimates of  f0 and f1 are calculated 
#' with a given regression coefficients which are efficient estimates  provided 
#' by other semiparametric methods such as logistic regression.
#' If support is (a,b) then replace r(x) by r[a+(b-a)x].
#' For a fixed \code{m}, using the Bernstein polynomial model for baseline \eqn{f_0},
#' MABLEs of \eqn{f_0} and parameters alpha can be estimated by EM algorithm and Newton  
#' iteration. If estimated lower bound \eqn{m_b} for \code{m} based on \code{y}
#' is smaller that that based on \code{x}, then switch \code{x} and \code{y} and
#' \eqn{f_1} is used as baseline. If \code{M=m} or \code{m0=m1=m}, then \code{m} is a 
#' preselected degree. If \code{m0<m1} it specifies the set of consective  
#' candidate model degrees \code{m0:m1} for searching an optimal degree by 
#' the change-point method, where \code{m1-m0>3}.  
#' @return  A list with components
#' \itemize{
#'   \item \code{m} the given or a selected degree by method of change-point
#'   \item \code{p} the estimated vector of mixture proportions \eqn{p = (p_0, \ldots, p_m)}
#'       with the given or selected degree \code{m}
#'   \item \code{alpha} the given regression coefficients
#'   \item \code{mloglik}  the maximum log-likelihood at degree \code{m}
#'   \item \code{interval} support/truncation interval \code{(a,b)}
#'   \item \code{baseline} ="control" if \eqn{f_0} is used as baseline, 
#'     or ="case" if \eqn{f_1} is used as baseline.
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1}, if greater than \code{m0}, is the
#'      largest candidate when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references Guan, Z., Maximum Approximate Bernstein Likelihood Estimation of 
#'   Densities in a Two-sample Semiparametric Model
#' @importFrom stats binomial glm 
#' @export

maple.dr<-function(x, y, M, regr, ..., interval = c(0,1), alpha=NULL, vb=0, 
      baseline=NULL, controls = mable.ctrl(),  progress=TRUE, message=TRUE){
    a<-interval[1]
    b<-interval[2]
    eps.em=controls$eps.em
    eps.nt=controls$eps.nt
    maxit.em=controls$maxit.em 
    maxit.nt=controls$maxit.nt
    tini=controls$tini
    level=controls$sig.level
    if(a>=b) stop("'a' must be smaller than 'b'")
    x1<-(x-a)/(b-a)
    y1<-(y-a)/(b-a)
    regr <- match.fun(regr)
    ff <- function(x) regr(a+(b-a)*x,...)
    ry<-ff(y1)
    nx<-length(x)
    ny<-length(y)
    if(ny==0 || ny==1) stop("'case' data 'y' are missing or too small.\n")
    if(is.null(alpha)){
      z<-c(rep(0,nx),rep(1,ny))
      alpha=glm(z~ff(c(x1,y1))[,-1], family= binomial)$coefficients
      alpha[1]=alpha[1]+log(nx/ny)
    }
    d<-length(alpha)-1
    ry<-ff(y1)
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
    if(!is.null(baseline)){
      if(baseline!="Control" && baseline!="Case") stop("Invalid 'baseline'.\n")
      m0<-M[1]}
    else{
      xbar<-mean(x1); sx2<-var(x1); 
      mx<-max(0,ceiling(xbar*(1-xbar)/sx2-3)-2)
      ybar<-mean(y1); sy2<-var(y1); 
      my<-max(0,ceiling(ybar*(1-ybar)/sy2-3)-2)
      if(nx>0){
        xbar<-mean(x1); sx2<-var(x1); 
        mx<-max(0,ceiling(xbar*(1-xbar)/sx2-3)-2)
        m0<-min(mx,my)
      }else{
        m0<-my
        mx<--1}
      if(mx<=my) baseline<-"Control"
      else baseline<-"Case"}
    if(M[1]<m0){
      message("Replace M[1]=",M[1]," by the recommended ", m0,".")
      M[1]=m0; M[2]=max(m0,M[2])}
    k<-M[2]-M[1]
    if(k>0 && k<=3){
      message("Too few candidate model degrees. Add 5 more.")
      M[2]<-M[2]+5}
    #p<-rep(1,M[1]+1)/(M[1]+1)
    ## Call maple_dr
    if(baseline=="Control"){
      message("Use 'control' as baseline.")
      wk<-.External("maple_dr", ff, rho = environment(),
        as.double(x1), as.double(y1), as.double(ry), as.integer(M), 
        as.double(alpha), as.integer(d), as.integer(nx),   
        as.integer(ny), as.double(eps.em), as.double(eps.nt), as.integer(maxit.em),
        as.integer(maxit.nt), as.double(tini), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
    }else{
      #baseline<-"Case"
      message("Use 'case' as baseline.")
      wk<-.External("maple_dr", ff, rho = environment(),
        as.double(y1), as.double(x1), as.double(ff(x1)), as.integer(M), 
        as.double(-alpha), as.integer(d), as.integer(ny),   
        as.integer(nx), as.double(eps.em), as.double(eps.nt), as.integer(maxit.em),
        as.integer(maxit.nt), as.double(tini), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
    }
    res <- c(list(regr=ff),list(x=x,y=y, interval=interval, 
          mloglik=wk$lk[wk$m-M[1]+1], baseline=baseline),
          wk[c("lk", "lr", "p", "m", "wt", "alpha", "chpts", "pval", "M")])
    class(res) <- "mable_dr"
    res
}
### MABLE DR MODEL for Grouped Data

#####################################################
#    Density ratio model based on grouped data      #
#####################################################
#' Mable fit of the density ratio model based on grouped data
#' @description Maximum approximate Bernstein/Beta likelihood estimation in a
#'   density ratio model based on two-sample grouped data.
#' @param t cutpoints of class intervals
#' @param n0,n1 frequencies of two sample data grouped by the classes 
#'       specified by \code{t}. \code{n0}:"Control", \code{n1}: "Case".
#' @param M a positive integer or a vector \code{(m0, m1)}. 
#' @param regr regressor vector function \eqn{r(x)=(1,r_1(x),...,r_d(x))} 
#'     which returns n x (d+1) matrix, n=length(x)
#' @param ... additional arguments to be passed to regr
#' @param interval a vector \code{(a,b)} containing the endpoints of 
#'   supporting/truncation interval of x and y.
#' @param alpha a given regression coefficient, missing value is imputed by 
#'     logistic regression
#' @param vb code for vanishing boundary constraints, -1: f0(a)=0 only, 
#'     1: f0(b)=0 only, 2: both, 0: none (default).
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#'    and the convergence criterion for EM and Newton iterations. Default is 
#'    \code{\link{mable.ctrl}}. See Details.
#' @param progress logical: should a text progressbar be displayed
#' @param message logical: should warning messages be displayed
#' @details 
#' Suppose that \code{n0} ("control") and \code{n1} ("case") are frequencies of  
#' independent samples grouped by the classes \code{t} from f0 and f1 which   
#' satisfy f1(x)=f0(x)exp[alpha0+alpha'r(x)] with r(x)=(r1(x),...,r_d(x)). Maximum 
#' approximate Bernstein/Beta likelihood estimates of (alpha0,alpha), f0 and f1 
#' are calculated. If support is (a,b) then replace r(x) by r[a+(b-a)x].
#' For a fixed \code{m}, using the Bernstein polynomial model for baseline \eqn{f_0},
#' MABLEs of \eqn{f_0} and parameters alpha can be estimated by EM algorithm and Newton  
#' iteration. If estimated lower bound \eqn{m_b} for \code{m} based on \code{n1}
#' is smaller that that based on \code{n0}, then switch \code{n0} and \code{n1} and
#' use \eqn{f_1} as baseline. If \code{M=m} or \code{m0=m1=m}, then \code{m} is a 
#' preselected degree. If \code{m0<m1} it specifies the set of consective  
#' candidate model degrees \code{m0:m1} for searching an optimal degree by 
#' the change-point method, where \code{m1-m0>3}.  
#' @importFrom stats binomial glm 
#' @export

mable.dr.group<-function(t, n0, n1, M, regr, ..., interval=c(0,1), alpha=NULL, 
     vb=0, controls=mable.ctrl(), progress=TRUE, message=TRUE){
#    dyn.load("mable-dr-model")
    a<-interval[1]
    b<-interval[2]
    eps.em=controls$eps.em
    eps.nt=controls$eps.nt
    maxit.em=controls$maxit.em 
    maxit.nt=controls$maxit.nt
    tini=controls$tini
    level=controls$sig.level
    nx<-sum(n0)
    ny<-sum(n1)
    N<-length(t)-1
    n<-nx+ny
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
    if(a>=b) stop("'a' must be smaller than 'b'")
    t1<-(t-a)/(b-a)
    if(is.null(n1))stop("'n1' is missing.")
    y1<-rep.int((t1[-N]+t1[-1])/2,times=n1) 
    ybar<-mean(y1); sy2<-var(y1); 
    my<-max(0,ceiling(ybar*(1-ybar)/sy2-3)-2)
    #m0<-max(mx,my)
    #if(M[1]<m0){
    #  message("Replace M[1]=",M[1]," by the recommended ", m0,".")
    #  M[1]=m0; M[2]=max(m0,M[2])}
    if(is.null(n0)){
      if(is.null(alpha)) stop("'alpha' is missing.")
      m0<-my
      mx<--1
    }else{
      x1<-rep.int((t1[-N]+t1[-1])/2,times=n0)
      xbar<-mean(x1); sx2<-var(x1); 
      mx<-max(0,ceiling(xbar*(1-xbar)/sx2-3)-2)
      m0<-min(mx,my)}
    if(M[1]<m0){
      message("Replace M[1]=",M[1]," by the recommended ", m0,".")
      M[1]=m0; M[2]=max(m0,M[2])}
    k<-M[2]-M[1]
    if(k>0 && k<=3){
      message("Too few candidate model degrees. Add 5 more.")
      M[2]<-M[2]+5}
    regr <- match.fun(regr)
    ff <- function(x) regr(a+(b-a)*x,...)
    if(is.null(alpha)){
      # logistic regression
      rho<-ny/nx; z<-c(rep(0,nx),rep(1,ny))
      ry<-ff(y1)
      alpha<-glm(z ~ ff(c(x1,y1))[,-1], family=binomial)$coefficients
      alpha[1]<-alpha[1]-log(rho)
    }
    d<-length(alpha)-1
    ## Call C_mable_dr_group
    if(mx<=my){
      baseline<-"Control"
      message("Use 'control' as baseline.")
      ans<-.External("C_mable_dr_group", ff, rho = environment(),
        as.double(alpha),  as.double(t1), as.integer(n0), 
        as.integer(n1), as.integer(nx), as.integer(ny), as.integer(N), 
        as.integer(M), as.integer(d), as.double(eps.em), as.double(eps.nt),  
        as.integer(maxit.em), as.integer(maxit.nt), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
    }else{
      baseline<-"Case"
      message("Use 'case' as baseline.")
      ans<-.External("C_mable_dr_group", ff, rho = environment(),
        as.double(alpha),  as.double(t1), as.integer(n1), 
        as.integer(n0), as.integer(ny), as.integer(nx), as.integer(N), 
        as.integer(M), as.integer(d), as.double(eps.em), as.double(eps.nt),  
        as.integer(maxit.em), as.integer(maxit.nt), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
    }
    res<-c(list(regr=ff), list(t=t, n0=n0, n1=n1, interval=interval, baseline=baseline,
        mloglik=ans$lk[ans$m-M[1]+1]), 
        ans[c("lk","lr", "alpha", "se", "p", "m", "chpts", "pval", "M")])
    class(res) <- "mable_dr"
    res
}
# Choosing optimal model degree using a fixed alpha estimated by logistic regression
# t: partition of [a,b]
#' Maximum approximate profile likelihood estimate of the density ratio model 
#'  for grouped data with given regression coefficients
#' @description Select optimal degree of Bernstein polynomial model for grouped data
#'   with a given regression coefficients. 
#' @param t cutpoints of class intervals
#' @param n0,n1 frequencies of two sample data grouped by the classes 
#'       specified by \code{t}. \code{n0}:"Control", \code{n1}: "Case".
#' @param M a positive integer or a vector \code{(m0, m1)}. 
#' @param regr regressor vector function \eqn{r(x)=(1,r_1(x),...,r_d(x))} 
#'     which returns n x (d+1) matrix, n=length(x)
#' @param ... additional arguments to be passed to regr
#' @param interval a vector \code{(a,b)} containing the endpoints of 
#'   supporting/truncation interval of x and y.
#' @param alpha a given regression coefficient, missing value is imputed by 
#'     logistic regression
#' @param vb code for vanishing boundary constraints, -1: f0(a)=0 only, 
#'     1: f0(b)=0 only, 2: both, 0: none (default).
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#'    and the convergence criterion for EM and Newton iterations. Default is 
#'    \code{\link{mable.ctrl}}. See Details.
#' @param progress logical: should a text progressbar be displayed
#' @param message logical: should warning messages be displayed
#' @details 
#' Suppose that \code{n0}("control") and \code{n1}("case") are frequencies of  
#' independent samples grouped by the classes \code{t} from f0 and f1 which   
#' satisfy f1(x)=f0(x)exp[alpha0+alpha'r(x)] with r(x)=(r1(x),...,r_d(x)). Maximum 
#' approximate Bernstein/Beta likelihood estimates of f0 and f1 are calculated 
#' with a given regression coefficients which are efficient estimates provided 
#' by other semiparametric methods such as logistic regression.
#'  If support is (a,b) then replace r(x) by r[a+(b-a)x].
#' For a fixed \code{m}, using the Bernstein polynomial model for baseline \eqn{f_0},
#' MABLEs of \eqn{f_0} and parameters alpha can be estimated by EM algorithm and Newton  
#' iteration. If estimated lower bound \eqn{m_b} for \code{m} based on \code{n1}
#' is smaller that that based on \code{n0}, then switch \code{n0} and \code{n1} and
#' use \eqn{f_1} as baseline. If \code{M=m} or \code{m0=m1=m}, then \code{m} is a 
#' preselected degree. If \code{m0<m1} it specifies the set of consective  
#' candidate model degrees \code{m0:m1} for searching an optimal degree by 
#' the change-point method, where \code{m1-m0>3}.  
#' @return  A list with components
#' \itemize{
#'   \item \code{m} the given or a selected degree by method of change-point
#'   \item \code{p} the estimated vector of mixture proportions \eqn{p = (p_0, \ldots, p_m)}
#'       with the given or selected degree \code{m}
#'   \item \code{alpha} the given regression coefficients
#'   \item \code{mloglik}  the maximum log-likelihood at degree \code{m}
#'   \item \code{interval} support/truncation interval \code{(a,b)}
#'   \item \code{baseline} ="control" if \eqn{f_0} is used as baseline, 
#'     or ="case" if \eqn{f_1} is used as baseline.
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1}, if greater than \code{m0}, is the
#'      largest candidate when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references Guan, Z., Application of Bernstein Polynomial Model to Density 
#'  and ROC Estimation in a Semiparametric Density Ratio Model
#' @importFrom stats binomial glm 
#' @export

maple.dr.group<-function(t, n0, n1, M, regr, ..., interval=c(0,1), alpha=NULL, 
     vb=0, controls=mable.ctrl(), progress=TRUE, message=TRUE){
#    dyn.load("mable-dr-model")
    a<-interval[1]
    b<-interval[2]
    eps.em=controls$eps.em
    eps.nt=controls$eps.nt
    maxit.em=controls$maxit.em 
    maxit.nt=controls$maxit.nt
    tini=controls$tini
    level=controls$sig.level
    nx<-sum(n0)
    ny<-sum(n1)
    N<-length(t)-1
    n<-nx+ny
    d<-length(alpha)-1
    if(a>=b) stop("'a' must be smaller than 'b'")
    t1<-(t-a)/(b-a)
    x1<-rep.int((t1[-N]+t1[-1])/2,times=n0)
    y1<-rep.int((t1[-N]+t1[-1])/2,times=n1) 
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
    xbar<-mean(x1); sx2<-var(x1); 
    mx<-max(0,ceiling(xbar*(1-xbar)/sx2-3)-2)
    ybar<-mean(y1); sy2<-var(y1); 
    my<-max(0,ceiling(ybar*(1-ybar)/sy2-3)-2)
    #m0<-max(mx,my)
    #if(M[1]<m0){
    #  message("Replace M[1]=",M[1]," by the recommended ", m0,".")
    #  M[1]=m0; M[2]=max(m0,M[2])}
    m0<-min(mx,my)
    if(M[1]<m0){
      message("Replace M[1]=",M[1]," by the recommended ", m0,".")
      M[1]=m0; M[2]=max(m0,M[2])}
    k<-M[2]-M[1]
    if(k>0 && k<=3){
      message("Too few candidate model degrees. Add 5 more.")
      M[2]<-M[2]+5}
    regr <- match.fun(regr)
    ff <- function(x) regr(a+(b-a)*x,...)
    if(is.null(alpha)){
      # logistic regression
      rho<-ny/nx; z<-c(rep(0,nx),rep(1,ny))
      ry<-ff(y1)
      alpha<-glm(z ~ ff(c(x1,y1))[,-1], family=binomial)$coefficients
      alpha[1]<-alpha[1]-log(rho)
    }
    d<-length(alpha)-1
    ## Call maple_dr_group
    if(mx<=my){
      baseline<-"Control"
      message("Use 'control' as baseline.")
      ans<-.External("maple_dr_group", ff, rho = environment(),
        as.double(alpha), as.double(t), as.integer(n0), 
        as.integer(n1), as.integer(nx),as.integer(ny), as.integer(N), 
        as.integer(M), as.integer(d), as.double(eps.em), as.double(eps.nt), 
        as.integer(maxit.em), as.integer(maxit.nt), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
    }else{
      baseline<-"Case"
      message("Use 'case' as baseline.")
      ans<-.External("maple_dr_group", ff, rho = environment(),
        as.double(alpha), as.double(t), as.integer(n1), 
        as.integer(n0), as.integer(ny),as.integer(nx), as.integer(N), 
        as.integer(M), as.integer(d), as.double(eps.em), as.double(eps.nt), 
        as.integer(maxit.em), as.integer(maxit.nt), as.integer(progress), 
        as.double(level), as.integer(message), as.integer(vb))
    }
    res<-c(list(regr=ff), list(t=t, n0=n0, n1=n1,interval=interval, baseline=baseline, 
      mloglik=ans$lk[ans$m-M[1]+1]),ans[c("lk","lr", "p", "m", "chpts", "pval", "M")])
#    dyn.unload("mable-dr-model")
    class(res) <- "mable_dr"
    res
}
 
########################################################
## Exponentially tilted  mixture of beta distributions #
########################################################
#' Exponentially Tilted Mixture Beta Distribution
#' @description Density, distribution function, quantile function and
#' pseudorandom number generation for the exponentially tilted  mixture of 
#' beta distributions, with shapes \eqn{(i+1, m-i+1)}, \eqn{i = 0, \ldots, m},  
#' given mixture proportions \eqn{p=(p_0,\ldots,p_m)} and support \code{interval}.
#' @param x a vector of quantiles
#' @param u a vector of probabilities
#' @param n sample size
#' @param p a vector of \code{m+1} components of \code{p} must be nonnegative 
#'   and sum to one for mixture beta distribution. See 'Details'.
#' @param alpha regression coefficients  
#' @param interval support/truncation interval \code{[a, b]}.
#' @param regr regressor vector function \eqn{r(x)=(1,r_1(x),...,r_d(x))} 
#'     which returns n x (d+1) matrix, n=length(x)
#' @param ... additional arguments to be passed to regr
#' @return A vector of \eqn{f_m(x; p)} or \eqn{F_m(x; p)} values at \eqn{x}.
#' \code{dmixbeta} returns the density, \code{pmixbeta} returns the cumulative
#'  distribution function, \code{qmixbeta} returns the quantile function, and
#' \code{rmixbeta}  generates pseudo random numbers.
#' @details
#'  The density of the mixture exponentially tilted beta distribution on an 
#'  interval \eqn{[a, b]} can be written \eqn{f_m(x; p)=(b-a)^{-1}\exp(\alpha'r(x))
#'  \sum_{i=0}^m p_i\beta_{mi}[(x-a)/(b-a)]/(b-a)},
#'  where \eqn{p = (p_0, \ldots, p_m)}, \eqn{p_i\ge 0}, \eqn{\sum_{i=0}^m p_i=1} and
#'  \eqn{\beta_{mi}(u) = (m+1){m\choose i}u^i(1-x)^{m-i}}, \eqn{i = 0, 1, \ldots, m},
#'  is the beta density with shapes \eqn{(i+1, m-i+1)}. The cumulative distribution
#' function is \eqn{F_m(x; p) = \sum_{i=0}^m p_i B_{mi}[(x-a)/(b-a);alpha]}, where
#' \eqn{B_{mi}(u ;alpha)}, \eqn{i = 0, 1, \ldots, m}, is the exponentially tilted 
#'  beta cumulative distribution function with shapes \eqn{(i+1, m-i+1)}. 
#' @author Zhong Guan <zguan@iu.edu>
#' @references Guan, Z., Application of Bernstein Polynomial Model to Density 
#'  and ROC Estimation in a Semiparametric Density Ratio Model
#' @keywords distribution  models  nonparametric  smooth
#' @concept Bernstein polynomial model
#' @concept Mixture beta distribution
#' @seealso \code{\link{mable}}
#' @examples
#' # classical Bernstein polynomial approximation
#' a<--4; b<-4; m<-200
#' x<-seq(a,b,len=512)
#' u<-(0:m)/m
#' p<-dnorm(a+(b-a)*u)
#' plot(x, dnorm(x), type="l")
#' lines(x, (b-a)*dmixbeta(x, p, c(a, b))/(m+1), lty=2, col=2)
#' legend(a, dnorm(0), lty=1:2, col=1:2, c(expression(f(x)==phi(x)),
#'                expression(B^{f}*(x))))
#'
#' @importFrom stats runif uniroot rbeta
#' @export

dtmixbeta<-function(x, p, alpha, interval=c(0,1), regr,...){
#    require(mable)
    d<-length(alpha)-1
    m<-length(p)-1
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    if(any(x<a) || any(x>b)) stop("Argument 'x' must be in 'interval'.")
    if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
    if(any(p<0)) message("Argument 'p' has negative component(s).")
    #ff <- function(x) regr((x-a)/(b-a),...)
    f0<-dmixbeta(x, p, interval)
    f1<-as.vector(f0*exp(regr(x)%*%alpha))
    return(f1)
}
# 
# cdf of mixture of exponentially tilted beta distributions
#
#' @rdname dtmixbeta
#' @export

ptmixbeta<-function(x, p, alpha, interval=c(0,1), regr,...){
    d<-length(alpha)-1
    m<-length(p)-1
    nx<-length(x)
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    if(any(x<a) || any(x>b)) stop("Argument 'x' must be in 'interval'.")
    if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
    if(any(p<0)) message("Argument 'p' has negative component(s).")
    x1<-(x-a)/(b-a)
#    regr <- match.fun(regr)
    ff <- function(x) regr(a+(b-a)*x,...)
#    ff <- function(x) regr((x-a)/(b-a),...)
    res<-.External("mixtbeta_cdf", ff, rho = environment(),
      as.double(alpha), as.double(p), as.double(x1), as.integer(d), 
      as.integer(m), as.integer(nx))
    y<-res$y 
    return(y)
}
#########################################################
#' @rdname dtmixbeta
#' @export
qtmixbeta<-function(u, p, alpha, interval=c(0, 1), regr,...){
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    if(any(u<0) || any(u>1)) stop("Argument 'u' must be in [0,1].\n")
    if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
    if(any(p<0)) message("Argument 'p' has negative component(s).")
    m<-length(p)-1
    n<-length(u)
    Q<-a+(b-a)*u
    for(i in 1:n){
      fn<-function(x) ptmixbeta(x, p, alpha, interval, regr,...)-u[i]
      if(u[i]>0 && u[i]<1) Q[i]<-uniroot(fn, interval)$root
    }
    #Q<-a+(b-a)*Q
    
    return(Q)
}

# generating prn from sum(p[i]*beta[m,i](x): i=0,...,m)*exp(a0+a[-1]^r(x))
#' @rdname dtmixbeta
#' @export

rtmixbeta<-function(n, p, alpha, interval=c(0,1), regr, ...){
   m<-length(p)-1
   u<-runif(n)
   a<-interval[1]
   b<-interval[2]
   if(a>=b) stop("a must be smaller than b")
   if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
   if(any(p<0)) message("Argument 'p' has negative component(s).")
   #obj<-list(regr=regr, p=p, alpha=alpha)
   y<-a+(b-a)*u
   for(i in 1:n){
      fun<-function(x) ptmixbeta(x, p, alpha, interval, regr,...)-u[i]
      if(u[i]>0 && u[i]<1) y[i]<-uniroot(fun, interval)$root
   }
   #a+(b-a)*y
   y
}
#rmixbeta(n, p)
#rtmixbeta(n, p, alpha, regr)
