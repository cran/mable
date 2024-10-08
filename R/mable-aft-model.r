##########################################################
#   MABLE of AFT model based on interval-censored data   #
##########################################################
#     AFT model with covariate for interval-censored     #
#  failure time data: S(t|x)=S(t exp(-gamma'(x-x0))|x0)   #
##########################################################
#  tau:  S(tau|x) ~ 0 for all x.  
# Data transformation: dividing all the times by b to obtain y=(y1,y2].
#   Let f(t|x) be the density of the transformed time given X=x.
#  We can approx f(t|x0) and S(t|x0) on [0,1], respectively, by
#    fm(t|x0;p)=sum_{i=0}^m p_i beta_{mi}(t)
#    Sm(t|x0;p)=sum_{i=0}^{m} p_i BS_{mi}(t),   
#       BS_{mi}(t)=int_t^1 beta_{mi}(s)ds.
##################################################################
# R CMD SHLIB mable-aft-model.c
# R --arch x64 CMD SHLIB mable-aft-model.c
# setwd("C:\\Users\\zguan\\Documents\\papers\\bernstein polynomials\\survival function\\C")
#    dyn.load("mable-aft-model")
#' Mable fit of Accelerated Failure Time Model
#' @description Maximum approximate Bernstein/Beta likelihood estimation for
#' accelerated failure time model based on interval censored data.
#' @param formula regression formula. Response must be \code{cbind}. See 'Details'.
#' @param data a data frame containing variables in \code{formula}.
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m0} or \code{m0 = m1 = m},  
#'   then \code{m0} is a preselected degree. If \code{m0 < m1} it specifies the set of 
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.  
#' @param g a \eqn{d}-vector of regression coefficients, default is the zero vector.
#is the regression coefficients
#   with response the logarithm of the midpoint of the observed interval. 
#' @param p an initial coefficients of Bernstein polynomial of degree \code{m0}, 
#'   default is the uniform initial. 
#' @param tau the right endpoint of the support or truncation interval \eqn{[0,\tau)} of the   
#'   baseline density. Default is \code{NULL} (unknown), otherwise if \code{tau} is given 
#'   then it is taken as a known value of \eqn{\tau}.  See 'Details'. 
#' @param x0 a data frame specifying working baseline covariates on the right-hand-side of \code{formula}. See 'Details'.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @details
#' Consider the accelerated failure time model with covariate for interval-censored failure time data: 
#' \eqn{S(t|x) = S(t \exp(\gamma^T(x-x_0))|x_0)}, where \eqn{x} and \eqn{x_0} may
#' contain dummy variables and interaction terms.  The working baseline \code{x0} in arguments
#' contains only the values of terms excluding dummy variables and interaction terms 
#' in the right-hand-side of \code{formula}. Thus \code{g} is the initial guess of 
#' the coefficients \eqn{\gamma} of \eqn{x-x_0} and could be longer than \code{x0}.   
#'   Let \eqn{f(t|x)} and \eqn{F(t|x) = 1-S(t|x)} be the density and cumulative distribution
#' functions of the event time given \eqn{X = x}, respectively.
#' Then \eqn{f(t|x_0)} on a truncation interval \eqn{[0, \tau]} can be approximated by  
#' \eqn{f_m(t|x_0; p) = \tau^{-1}\sum_{i=0}^m p_i\beta_{mi}(t/\tau)},
#' where \eqn{p_i\ge 0}, \eqn{i = 0, \ldots, m}, \eqn{\sum_{i=0}^mp_i=1},  
#' \eqn{\beta_{mi}(u)} is the beta denity with shapes \eqn{i+1} and \eqn{m-i+1}, and
#' \eqn{\tau} is larger than the largest observed time, either uncensored time, or right endpoint of interval/left censored,
#' or left endpoint of right censored time. So we can approximate  \eqn{S(t|x_0)} on \eqn{[0, \tau]} by
#' \eqn{S_m(t|x_0; p) = \sum_{i=0}^{m} p_i \bar B_{mi}(t/\tau)}, where \eqn{\bar B_{mi}(u)} is
#' the beta survival function with shapes \eqn{i+1} and \eqn{m-i+1}.
#'
#' Response variable should be of the form \code{cbind(l, u)}, where \code{(l,u)} is the interval 
#' containing the event time. Data is uncensored if \code{l = u}, right censored 
#' if \code{u = Inf} or \code{u = NA}, and  left censored data if \code{l = 0}.
#' The truncation time \code{tau} and the baseline \code{x0} should be chosen so that 
#' \eqn{S(t|x)=S(t \exp(\gamma^T(x-x_0))|x_0)} on \eqn{[\tau, \infty)} is negligible for
#' all the observed \eqn{x}.
#'
# For general interval censored observations, we keep the 
# right-censored but replace the finite interval with its midpoint and fit the data by 
# \code{aftsrr()} as a right-censored data. 
#' 
#'  The search for optimal degree \code{m} stops if either \code{m1} is reached or the test 
#'  for change-point results in a p-value \code{pval} smaller than \code{sig.level}.
#' @return A list with components
#' \itemize{ 
#'   \item \code{m} the given or selected optimal degree \code{m}
#'   \item \code{p} the estimate of \code{p = (p_0, \dots, p_m)}, the coefficients of Bernstein polynomial of degree \code{m}
#'   \item \code{coefficients} the estimated regression coefficients of the AFT model
#'   \item \code{SE} the standard errors of the estimated regression coefficients 
#'   \item \code{z} the z-scores of the estimated regression coefficients 
#'   \item \code{mloglik} the maximum log-likelihood at an optimal degree \code{m}
#'   \item \code{tau.n} maximum observed time \eqn{\tau_n}
#'   \item \code{tau} right endpoint of trucation interval \eqn{[0, \tau)}
#'   \item \code{x0} the working baseline covariates 
#'   \item \code{egx0} the value of \eqn{e^{\gamma^T x_0}} 
#'   \item \code{convergence} an integer code: 0 indicates a successful completion; 
#'     1 indicates that the search of an optimal degree using change-point method reached  
#'     the maximum candidate degree; 2 indicates that the matimum iterations was reached for
#'     calculating \eqn{\hat p} and \eqn{\hat\gamma} with the selected degree \eqn{m},  
#'     or the divergence of the last EM-like iteration for \eqn{p} or the divergence of
#'     the last (quasi) Newton iteration for \eqn{\gamma}; 3 indicates 1 and 2.   
#'   \item \code{delta} the final \code{delta} if \code{m0 = m1} or the final \code{pval} of the change-point 
#'      for searching the optimal degree \code{m};
#'  }
#'  and, if \code{m0<m1},
#' \itemize{
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last candidate when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references 
#' Guan, Z. (2019) Maximum Approximate Likelihood Estimation in Accelerated Failure Time Model for Interval-Censored Data, 
#' arXiv:1911.07087.
#' @examples \donttest{
#' ## Breast Cosmesis Data
#'   g <- 0.41 #Hanson and  Johnson 2004, JCGS
#'   aft.res<-mable.aft(cbind(left, right)~treat, data=cosmesis, M=c(1, 30), 
#'               g=g, tau=100, x0=data.frame(treat="RCT"))
#'   op<-par(mfrow=c(1,2), lwd=1.5)
#'   plot(x=aft.res, which="likelihood")
#'   plot(x=aft.res, y=data.frame(treat="RT"), which="survival", model='aft', type="l", col=1, 
#'       add=FALSE, main="Survival Function")
#'   plot(x=aft.res, y=data.frame(treat="RCT"), which="survival", model='aft', lty=2, col=1)
#'   legend("bottomleft", bty="n", lty=1:2, col=1, c("Radiation Only", "Radiation and Chemotherapy"))
#'   par(op)
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Accelerated failure time model 
#' @concept interval censoring
#' @seealso \code{\link{maple.aft}}
#' @importFrom stats coef lm reformulate terms
#' @importFrom survival Surv
#' @export
mable.aft<-function(formula, data, M, g=NULL, p=NULL, tau=NULL, x0=NULL,    
                 controls = mable.ctrl(), progress=TRUE){
  data.name<-deparse(substitute(data)) 
  fmla<-Reduce(paste, deparse(formula))
  Dta<-get.mableData(formula, data)
  x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
#  if(is.null(x0)) x0<-rep(0,d)
  delta<-Dta$delta; n<-length(y)
#   vars<-get.facMatrix(formula, data)
#   allvars<-vars$allvars  # all variables 
   allvars<-all.vars(formula)  # all variables 
  if(!is.null(x0)){
    if(length(x0)!=length(allvars)-2 || !is.data.frame(x0)){
        message("Baseline 'x0' must be a dataframe with names matching the RHS of 'formula'.\n")
        message("I will assign an 'x0' for you.\n")
        x0<-NULL}
    data0<-data.frame(cbind(0,1,x0))
    names(data0)<-c(allvars[1:2], names(x0))
    data1<-rbind(data[,names(data0)],data0)
#    data1<-rbind(data[,allvars],data0)
    Dta0<-get.mableData(formula, data1)
    x0<-as.matrix(Dta0$x)[n+1,]
  }
  if(is.null(x0)){
     for(i in 1:d) x0<-rep(0,length(allvars)-2)
   }
  x<-as.matrix(x)
  d<-ncol(x)
#  n<-length(y)
   if(is.null(g)){
     g<-rep(0,d)
    #g<--unname(coef(lm(log((y+min(y2,tau))/2)~ x))[-1])
    #cat("glsq=",g,"\n")
   }
   else if(length(g)!=d) stop("Invalid argument 'g'.")

  b<-max(y2[y2<Inf], y);
  n0<-sum(delta==0)
  n1<-sum(delta==1)
  n<-n0+n1
  ord<-order(delta)
  x<-as.matrix(x[ord,]); y<-y[ord]; y2<-y2[ord]
###
  i2<-(y2<Inf)
  xt<-t(t(x)-x0)
  known.tau<-FALSE
  if(is.null(tau)){ 
    tau<-max((y2*exp(xt%*%g))[i2,],y*exp(xt%*%g))+1/n #???    
  }
  else{
    #if(b>=tau) stop("tau must be greater than the maximum observed time")
    known.tau<-TRUE
    y<-y/tau; y2<-y2/tau
  }
  N<-c(n0,n1)
  dm<-c(d,0)
  conv<-0
  del<-0
  tau<-c(tau,b) # 
  y2[y2==Inf]<-.Machine$double.xmax/2
###
  Eps<-c(controls$eps.nt, controls$eps.em)
  MaxIt<-c(controls$maxit.nt, controls$maxit.em)
  ddell<-diag(0,d)
  if(missing(M) || length(M)==0) stop("'M' is missing.\n")
  else if(length(M)==1) M<-c(M,M)
  else if(length(M)>=2) {
      M<-c(min(M), max(M))
  }
  k<-M[2]-M[1]
  if(k>0 && k<=3) stop("Too few candidate model degrees.")
  ans<-list()
  ans$tau.n<-b
  ans$tau<-tau[1]
  ans$xNames<-Dta$xNames
  if(is.null(p) || length(p)!=M[1]+1) p<-rep(1,M[1]+1)/(M[1]+1)
  if(k==0){
    method<-0
    m<-M[1]    
    dm<-c(d,m)
    ell<-0
    ## Call C mable_aft_m
    res<-.C("mable_aft_m",
      as.double(g), as.double(p), as.integer(dm), as.double(x), as.double(y),  
      as.double(y2), as.double(tau), as.integer(N), as.double(x0), as.double(ell), 
      as.double(ddell), as.double(Eps), as.integer(MaxIt), as.logical(progress), 
      as.integer(conv), as.double(del), as.logical(known.tau), as.integer(method))
    ans$m<-m
    ans$mloglik<-res[[10]][1]
    ans$p<-res[[2]] 
    ans$x0<-res[[9]]
    ans$coefficients<-res[[1]]
    ans$egx0<-exp(sum(res[[1]]*res[[9]]))
    Sig<--n*matrix(res[[11]], nrow=d, ncol=d)
    ans$SE<-sqrt(diag(Sig)/n)
    ans$z<-res[[1]]/ans$SE 
    ans$convergence<-res[[15]]
    ans$delta<-res[[16]] 
    ans$tau<-res[[7]][1]
  }
  else{
    lk<-rep(0, k+1)
    lr<-rep(0, k)    
    pval<-rep(0,k+1)
    chpts<-rep(0,k+1)
    p<-c(p, rep(0, k))
    level<-controls$sig.level # it is also used to return delta
    ## Call C mable_aft
    res<-.C("mable_aft",
      as.integer(M), as.double(g), as.integer(dm), as.double(p), as.double(x),  
      as.double(y), as.double(y2), as.double(tau), as.integer(N), as.double(x0), 
      as.double(lk), as.double(lr), as.double(ddell), as.double(Eps), 
      as.integer(MaxIt), as.logical(progress), as.double(pval), as.integer(chpts), 
      as.double(level), as.integer(conv), as.logical(known.tau))
    M<-res[[1]]
    ans$x0<-res[[10]]
    ans$egx0<-exp(sum(res[[2]]*res[[10]]))
    Sig<--n*matrix(res[[13]], nrow=d, ncol=d)
    ans$coefficients<-res[[2]]
    ans$SE<-sqrt(diag(Sig)/n)
    ans$z<-res[[2]]/ans$SE
    ans$M<-M
    if(k<M[2]-M[1]) stop("k<M[2]-M[1]")
    k<-M[2]-M[1]
    if(k!=res[[3]][1]) stop("k!=res[[3]][1]")
    ans$lk<-res[[11]][1:(k+1)]   
    ans$lr<-res[[12]][1:k]
    ans$m<-res[[3]][2]
    mp1<-ans$m+1
    if(mp1>length(p)) stop("mp1>length(p)")
    if(mp1-M[1]<0 || mp1-M[1]>k+1) stop("error")
    ans$mloglik<-res[[11]][mp1-M[1]]  
    ans$p<-res[[4]][1:mp1]  
    ans$pval<-res[[17]][1:(k+1)] #???
    ans$chpts<-res[[18]][1:(k+1)]+M[1]
    ans$convergence<-res[[20]]
    ans$delta<-res[[19]]  #???  
    ans$tau<-res[[8]][1]
  }
  #ans$egx0<-1/ans$egx0
  ans$model<-"aft"
  ans$callText<-fmla
  ans$data.name<-data.name
#  ans$fmatrices<-vars[[1]]
#  ans$factors<-vars$factors
  ans$allvars<-allvars
#  ans$whichisfactor<-vars$whichisfactor
  class(ans)<-"mable_reg"
  return(ans)
}
######################################################################################
#  Maximum Approximate Profile Likelihood Estimation in AFT model with a given gamma
# M: set of positive integers as candidate degrees of Bernstein poly model
######################################################################################
#' Mable fit of AFT model with given regression coefficients
#' @param formula regression formula. Response must be \code{cbind}.  See 'Details'.
#' @param data a data frame containing variables in \code{formula}.
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m0} or \code{m0 = m1},  
#'   then \code{m0} is a preselected degree. If \code{m0 < m1} it specifies the set of 
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0 > 3}.  
#' @param g the given \eqn{d}-vector of regression coefficients. 
#' @param p an initial coefficients of Bernstein polynomial of degree \code{m0}, 
#'   default is the uniform initial. 
#' @param tau the right endpoint of the support or truncation interval \eqn{[0,\tau)} of the   
#'   baseline density. Default is \code{NULL} (unknown), otherwise if \code{tau} is given 
#'   then it is taken as a known value of \eqn{\tau}.  See 'Details'. 
#' @param x0 a data frame specifying working baseline covariates on the right-hand-side of \code{formula}. See 'Details'.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate profile likelihood estimation of Bernstein
#'  polynomial model in accelerated failure time based on interal 
#'  censored event time data with given regression coefficients which are efficient
#'  estimates provided by other semiparametric methods. 
#' @details
#' Consider the accelerated failure time model with covariate for interval-censored failure time data: 
#' \eqn{S(t|x) = S(t \exp(\gamma^T(x-x_0))|x_0)}, where \eqn{x} and \eqn{x_0} may
#' contain dummy variables and interaction terms.  The working baseline \code{x0} in arguments
#' contains only the values of terms excluding dummy variables and interaction terms 
#' in the right-hand-side of \code{formula}. Thus \code{g} is the initial guess of 
#' the coefficients \eqn{\gamma} of \eqn{x-x_0} and could be longer than \code{x0}. 
#'   Let \eqn{f(t|x)} and \eqn{F(t|x) = 1-S(t|x)} be the density and cumulative distribution
#' functions of the event time given \eqn{X = x}, respectively.
#' Then \eqn{f(t|x_0)} on a support or truncation interval \eqn{[0, \tau]} can be approximated by  
#' \eqn{f_m(t|x_0; p) = \tau^{-1}\sum_{i=0}^m p_i\beta_{mi}(t/\tau)},
#' where \eqn{p_i \ge 0}, \eqn{i = 0, \ldots, m}, \eqn{\sum_{i=0}^mp_i=1},  
#' \eqn{\beta_{mi}(u)} is the beta denity with shapes \eqn{i+1} and \eqn{m-i+1}, and
#' \eqn{\tau} is larger than the largest observed time, either uncensored time, or right endpoint of interval/left censored,
#' or left endpoint of right censored time. We can approximate  \eqn{S(t|x_0)} on \eqn{[0, \tau]} by
#' \eqn{S_m(t|x_0; p) = \sum_{i=0}^{m} p_i \bar B_{mi}(t/\tau)}, where \eqn{\bar B_{mi}(u)} is
#' the beta survival function with shapes \eqn{i+1} and \eqn{m-i+1}.
#'
#' Response variable should be of the form \code{cbind(l, u)}, where \code{(l,u)} is the interval 
#' containing the event time. Data is uncensored if \code{l = u}, right censored 
#' if \code{u = Inf} or \code{u = NA}, and  left censored data if \code{l = 0}.
#' The truncation time \code{tau} and the baseline \code{x0} should be chosen so that 
#' \eqn{S(t|x) = S(t \exp(\gamma^T(x-x_0))|x_0)} on \eqn{[\tau, \infty)} is negligible for
#' all the observed \eqn{x}.
#'
#'  The search for optimal degree \code{m} stops if either \code{m1} is reached or the test 
#'  for change-point results in a p-value \code{pval} smaller than \code{sig.level}.
#' @return A list with components
#' \itemize{ 
#'   \item \code{m} the selected optimal degree \code{m}
#'   \item \code{p} the estimate of \eqn{p=(p_0, \dots, p_m)}, the coefficients of Bernstein polynomial of degree \code{m}
#'   \item \code{coefficients} the given regression coefficients of the AFT model
#'   \item \code{SE} the standard errors of the estimated regression coefficients 
#'   \item \code{z} the z-scores of the estimated regression coefficients 
#'   \item \code{mloglik} the maximum log-likelihood at an optimal degree \code{m}
#'   \item \code{tau.n} maximum observed time \eqn{\tau_n}
#'   \item \code{tau} right endpoint of trucation interval \eqn{[0, \tau)}
#'   \item \code{x0} the working baseline covariates 
#'   \item \code{egx0} the value of \eqn{e^{\gamma^T x_0}} 
#'   \item \code{convergence} an integer code, 1 indicates either the EM-like 
#'     iteration for finding maximum likelihood reached the maximum iteration for at least one \code{m} 
#'     or the search of an optimal degree using change-point method reached the maximum candidate degree,
#'     2 indicates both occured, and 0 indicates a successful completion.  
#'   \item \code{delta} the final \code{delta} if \code{m0 = m1} or the final \code{pval} of the change-point 
#'      for searching the optimal degree \code{m};
#' }
#'  and, if \code{m0<m1},
#' \itemize{
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last candidate when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references 
#' Guan, Z. (2019) Maximum Approximate Likelihood Estimation in Accelerated Failure Time Model for Interval-Censored Data, 
#' arXiv:1911.07087.
#' @examples \donttest{
#' ## Breast Cosmesis Data
#'   g<-0.41 #Hanson and  Johnson 2004, JCGS, 
#'   res1<-maple.aft(cbind(left, right)~treat, data=cosmesis, M=c(1,30),  g=g, 
#'                tau=100, x0=data.frame(treat="RCT"))
#'   op<-par(mfrow=c(1,2), lwd=1.5)
#'   plot(x=res1, which="likelihood")
#'   plot(x=res1, y=data.frame(treat="RT"), which="survival", model='aft', type="l", col=1, 
#'       add=FALSE, main="Survival Function")
#'   plot(x=res1, y=data.frame(treat="RCT"), which="survival", model='aft', lty=2, col=1)
#'   legend("bottomleft", bty="n", lty=1:2, col=1, c("Radiation Only", "Radiation and Chemotherapy"))
#'   par(op)
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Accelerated failure time model 
#' @concept interval censoring
#' @seealso \code{\link{mable.aft}} 
#' @export
maple.aft<-function(formula, data, M, g, tau=NULL, p=NULL, x0=NULL, 
          controls = mable.ctrl(), progress=TRUE){
  data.name<-deparse(substitute(data)) 
  fmla<-Reduce(paste, deparse(formula))
  if(missing(g)) stop("missing argument 'g'.")
  d<-length(g)
  Dta<-get.mableData(formula, data)
  x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
  delta<-Dta$delta; n<-length(y)
#  if(is.null(x0)) x0<-rep(0,d)
#  vars<-get.facMatrix(formula, data)
#  allvars<-vars$allvars  # all variables 
  allvars<-all.vars(formula)  # all variables 
  if(!is.null(x0)){
    if(length(x0)!=length(allvars)-2 || !is.data.frame(x0)){
        message("Baseline 'x0' must be a dataframe with names matching the RHS of 'formula'.\n")
        message("I will assign an 'x0' for you.\n")
        x0<-NULL}
    data0<-data.frame(cbind(0,1,x0))
    names(data0)<-c(allvars[1:2], names(x0))
    data1<-rbind(data[,names(data0)],data0)
#    data1<-rbind(data[,allvars],data0)
    Dta0<-get.mableData(formula, data1)
    x0<-as.matrix(Dta0$x)[n+1,]
  }
  if(is.null(x0)){
    x0<-rep(0,d)   
  }
  b<-max(y2[y2<Inf], y);
  #if(b>tau) stop("tau must be greater than or equal to the maximum observed time")
  #y<-y/tau; y2<-y2/tau
  #y2[y2==Inf]<-.Machine$double.xmax/2 
  #y2[y2==Inf]<-1 # truncation interval is [0, tau)
  x<-as.matrix(x)
  if(d!=ncol(x)) stop("Invalid argument 'g'.")
###
  i2<-(y2<Inf)
  if(is.null(tau)) {
    tau<-max((y2*exp(x%*%g))[i2,],y*exp(x%*%g))
    known.tau<-FALSE
  }
  else{
    if(b>tau) stop("tau must be greater than or equal to the maximum observed time")
    y<-y/tau; y2<-y2/tau 
    known.tau<-TRUE
  }
  tau<-c(tau,b)
  y2[y2==Inf]<-.Machine$double.xmax/2
###
  n<-length(y)
  n0<-sum(delta==0)
  n1<-sum(delta==1)
  N<-c(n0,n1) 
  ord<-order(delta)
  x<-x[ord,]; y<-y[ord]; y2<-y2[ord]
  if(missing(M) || length(M)==0) stop("'M' is missing.\n")
  else if(length(M)==1) M<-c(M,M)
  else if(length(M)>=2){
      M<-c(min(M), max(M))
  }
  k<-M[2]-M[1]
  if(k>0 && k<=3) stop("Too few candidate model degrees.")
  lk<-rep(0, k+1)
  lr<-rep(0, k)
  if(is.null(p)) p<-rep(1, M[1]+1)/(M[1]+1) 
  p<-c(p, rep(0,k))
  ddell<-diag(0,d)
  pval<-rep(0,k+1)
  chpts<-rep(0,k+1)
  level<-controls$sig.level
  dm<-c(d,0)
  conv<-0
  del<-0
  ## Call C mable_aft_gamma
  res<-.C("mable_aft_gamma",
    as.integer(M), as.double(g), as.integer(dm), as.double(x), as.double(y),  
    as.double(y2), as.integer(N), as.double(x0), as.double(lk), as.double(lr), 
    as.double(p), as.double(ddell), as.double(controls$eps), as.integer(controls$maxit), 
    as.logical(progress), as.double(pval), as.integer(chpts), as.double(level), 
    as.integer(conv), as.double(del), as.double(tau), as.logical(known.tau))
  ans<-list()
  M<-res[[1]]
  ans$x0<-res[[8]]
  ans$egx0<-exp(sum(g*x0))
  ans$m<-res[[3]][2]
  ans$tau.n<-b
  ans$tau<-res[[21]][1]
  k<-M[2]-M[1]
  lk<-res[[9]][1:(k+1)]-n0*log(b) 
  ans$mloglik<-lk[ans$m-M[1]+1]
  ans$coefficients<-g
  mp1<-ans$m+1
  ans$p<-res[[11]][1:mp1]
  ans$xNames<-Dta$xNames
  ans$convergence<-res[[19]]
  ans$delta<-res[[20]]
  if(k>0){
    ans$M<-M; ans$lk<-lk; ans$lr<-res[[10]][1:k]
    ans$pval<-res[[16]][1:(k+1)]; ans$chpts<-res[[17]][1:(k+1)]+M[1]}  
  ans$model<-"aft"
  ans$callText<-fmla
  ans$data.name<-data.name
#  ans$fmatrices<-vars[[1]]
#  ans$factors<-vars$factors
  ans$allvars<-allvars
#  ans$whichisfactor<-vars$whichisfactor
  class(ans)<-"mable_reg"
  return(ans)
}


