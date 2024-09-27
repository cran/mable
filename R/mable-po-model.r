############################################################
#                                                          #
#    MABLE for Proportional Odds Regression Model PO       #
#                                                          #
############################################################
#  MABLE for PO Model:                                     #
#        [1-S^eta(t|x)]/S(t|x)                             #
#        ------------------------- = exp[gamma'x]          #
#        [1-S^eta(t|0)]/S(t|0)                             #
#   with covariate for interval-censored failure time data #
#  where S(t|x) is the survival function given covariate x #
############################################################
#  

##################################################################
#' Mable fit of proportional odds rate regression model
#' @param formula regression formula. Response must be \code{cbind}.  See 'Details'.
#' @param data a data frame containing variables in \code{formula}.
#' @param M a positive integer or a vector \code{(m0, m1)}. 
#' If \code{M = m} or \code{m0 = m1 = m},
#'   then \code{m} is a preselected degree. If \code{m0<m1} it specifies the set of
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.
#' @param g an initial guess of \eqn{d}-vector of regression coefficients.  See 'Details'.
#' @param tau right endpoint of support \eqn{[0, \tau]} must be greater than or equal to the
#'        maximum observed time
#' @param x0 a data frame specifying working baseline covariates on the right-hand-side of \code{formula}. See 'Details'.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate Bernstein/Beta likelihood estimation in general 
#'     proportional odds regression model based on interal censored event time data.
#' @details
#' Consider PO model with covariate for interval-censored failure time data:
#' \eqn{[1-S(t|x)]/S(t|x) = \exp[\gamma'(x-x_0)][1-S(t|x_0)]/S(t|x_0)}, 
#'  where \eqn{x_0} satisfies \eqn{\gamma'(x-x_0)\ge 0}, where \eqn{x} and \eqn{x_0} may
#' contain dummy variables and interaction terms.  The working baseline \code{x0} in arguments
#' contains only the values of terms excluding dummy variables and interaction terms 
#' in the right-hand-side of \code{formula}. Thus \code{g} is the initial guess of 
#' the coefficients \eqn{\gamma} of \eqn{x-x_0} and could be longer than \code{x0}.
#'   Let \eqn{f(t|x)} and \eqn{F(t|x) = 1-S(t|x)} be the density and cumulative distribution
#' functions of the event time given \eqn{X = x}, respectively.
#' Then \eqn{f(t|x_0)} on \eqn{[0, \tau]} can be approximated by
#' \eqn{f_m(t|x_0, p) = \tau^{-1}\sum_{i=0}^m p_i\beta_{mi}(t/\tau)},
#' where \eqn{p_i \ge 0}, \eqn{i = 0, \ldots, m}, \eqn{\sum_{i=0}^mp_i = 1},
#' \eqn{\beta_{mi}(u)} is the beta denity with shapes \eqn{i+1} and \eqn{m-i+1}, and
#' \eqn{\tau} is the right endpoint of support interval of the baseline density. 
#'  We can approximate  \eqn{S(t|x_0)} on \eqn{[0,\tau]} by
#' \eqn{S_m(t|x_0; p) = \sum_{i=0}^{m} p_i \bar B_{mi}(t/\tau)}, where
#' \eqn{\bar B_{mi}(u)}, \eqn{i = 0, \ldots, m}, is the beta survival function with shapes
#'  \eqn{i+1} and \eqn{m-i+1}.
#'
#' Response variable should be of the form \code{cbind(l,u)}, where \code{(l,u)} is the interval
#' containing the event time. Data is uncensored if \code{l = u}, right censored
#' if \code{u = Inf} or \code{u = NA}, and  left censored if \code{l = 0}.
#' The associated covariate contains \eqn{d} columns. The baseline \code{x0} should chosen so
#' that \eqn{\gamma'(x-x_0)} is nonnegative for all the observed \eqn{x} and
#' all \eqn{\gamma} in a neighborhood of its true value.
#'
#' A missing initial value of \code{g} is imputed by \code{ic_sp()} of package \code{icenReg} 
#' with \code{model="po"}.
#'  The search for optimal degree \code{m} stops if either \code{m1} is reached or the test
#'  for change-point results in a p-value \code{pval} smaller than \code{sig.level}.
#' This process takes longer than \code{\link{maple.po}} to select an optimal degree.
#' @return A list with components
#' \itemize{
#'   \item \code{m} the selected/preselected optimal degree \code{m}
#'   \item \code{p} the estimate of \eqn{p = (p_0, \dots, p_m)}, the coefficients of
#'             Bernstein polynomial of degree \code{m}
#'   \item \code{coefficients} the estimated regression coefficients of the PO model
#'   \item \code{SE} the standard errors of the estimated regression coefficients
#'   \item \code{z} the z-scores of the estimated regression coefficients
#'   \item \code{mloglik} the maximum log-likelihood at an optimal degree \code{m}
#'   \item \code{tau.n} maximum observed time \eqn{\tau_n}
#'   \item \code{tau} right endpoint of support \eqn{[0, \tau]}
#'   \item \code{x0} the working baseline covariates
#'   \item \code{egx0} the value of \eqn{e^{\gamma'x_0}}
#'   \item \code{convergence} an integer code, 1 indicates either the EM-like iteration 
#'     for finding maximum likelihood reached the maximum iteration for at least one \code{m}
#'     or the search of an optimal degree using change-point method reached the maximum
#'     candidate degree, 2 indicates both occured, and 0 indicates a successful completion.
#'   \item \code{delta} the final \code{delta} if \code{m0 = m1} or the final \code{pval} of
#'      the change-point for searching the optimal degree \code{m};
#'  }
#'  and, if \code{m0<m1},
#' \itemize{
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last candidate degree
#'        when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0,\ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at 
#'        \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references
#' Guan, Z. Maximum Likelihood Estimation in Proportional Odds Regression Model 
#'            Based on Interval-Censored Event-time Data  
#' @examples
#' \donttest{
# # Ovarian Cancer Survival Data
# require(survival)
# require(icenReg)
# futime2<-ovarian$futime
# futime2[ovarian$fustat==0]<-Inf
# ovarian2<-data.frame(age=ovarian$age, futime1=ovarian$futime,
#      futime2=futime2)    
# ova.sp<-ic_sp(cbind(futime1,futime2) ~ age, data = ovarian2, model="po") 
# tau<-2000
# x0<-data.frame(age=min(ovarian$age))
# ova<-mable.po(cbind(futime1, futime2) ~ age, data = ovarian2,              
#      M=c(2,35), g=-ova.sp$coefficients, x0=x0, tau=tau)                             
# op<-par(mfrow=c(2,2))
# plot(ova, which = "likelihood")
# plot(ova, which = "change-point")
# plot(ova, y=data.frame(age=60), which="survival", add=FALSE, type="l",
#       xlab="Days", main="Age = 60")
# plot(ova, y=data.frame(age=65), which="survival", add=FALSE, type="l",
#       xlab="Days", main="Age = 65")
# par(op)
#
#' # Veteran's Administration Lung Cancer Data 
#' require(survival)
#' require(icenReg)
#' require(mable)
#' l<-veteran$time->u
#' u[veteran$status==0]<-Inf
#' veteran1<-data.frame(l=l, u=u, karno=veteran$karno, celltype=veteran$celltype, 
#'            trt=veteran$trt, age=veteran$age, prior=veteran$prior>0) 
#' fit.sp<-ic_sp(cbind(l,u) ~ karno+celltype, data = veteran1,  model="po") 
#' x0<-data.frame(karno=100, celltype="squamous")
#' tau<-2000
#' res<-mable.po(cbind(l,u) ~ karno+celltype, data = veteran1, M=c(1,35),                               
#'      g=-fit.sp$coefficients, x0=x0, tau=tau)                          
#' op<-par(mfrow=c(2,2))
#' plot(res, which = "likelihood")
#' plot(res, which = "change-point")
#' plot(res, y=data.frame(karno=20, celltype="squamous"), which="survival", 
#'       add=FALSE, type="l", xlab="Days", 
#'       main=expression(paste("Survival: ", bold(x)==0)))
#' plot(res, y=data.frame(karno=80, celltype="smallcell"), which="survival", 
#'       add=FALSE, type="l", xlab="Days", 
#'       main=expression(paste("Survival: ", bold(x)==bold(x)[0])))
#' par(op)
#' 
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Generalized proportional odds model
#' @concept interval censoring
#' @seealso \code{\link{maple.ph}}
# @useDynLib mable-ph-model .registration = TRUE
#' @importFrom icenReg ic_sp
#' @export
mable.po<-function(formula, data, M, g=NULL, tau, x0=NULL,
      controls = mable.ctrl(), progress=TRUE){
  data.name<-deparse(substitute(data))
  fmla<-Reduce(paste, deparse(formula))
  Dta<-get.mableData(formula, data)
  x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
  delta<-Dta$delta; n<-length(y)
  if(is.null(g)) g<-ic_sp(formula, data, model = 'po')$coefficients
  d<-length(g)
  allvars<-all.vars(formula) # all variables 
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
    for(i in 1:d) x0[i]<-ifelse(g[i]>=0, min(x[,i]), max(x[,i]))
  }
  b<-max(y2[y2<Inf], y)
  if(b>tau) stop("tau must be greater than or equal to the maximum observed time.")
  if(tau==Inf || tau<0) stop("'tau' must be finite and positive.")
  y2[y2==Inf]<-tau
  y<-y/tau; y2<-y2/tau
  if(!is.matrix(x)) x<-matrix(x, ncol=1)
  if(d!=ncol(x)) stop("Invalid argument 'g'.")
  #n<-length(y)
  n0<-sum(delta==0)
  n1<-sum(delta==1)
  n<-n0+n1
  N<-c(n0,n1)
  dm<-c(d,0)
  ddell<-diag(0,d)
  conv<-0L
  del<-0
  ord<-order(delta)
  x<-as.matrix(x[ord,]); y<-y[ord]; y2<-y2[ord]
  Eps<-c(controls$eps, controls$eps.em, controls$eps.nt, .Machine$double.eps)
  MaxIt<-c(controls$maxit, controls$maxit.em, controls$maxit.nt)
  level<-controls$sig.level
  if(missing(M) || length(M)==0) stop("'M' is missing.\n")
  else if(length(M)==1) M<-c(M,M)
  else if(length(M)>=2) M<-c(min(M), max(M))
  k<-M[2]-M[1]
  if(k>0 && k<=3) stop("Too few candidate model degrees.")
  eta=1
  eta.known=TRUE
  if(k==0){
    m<-M[1]
    dm<-c(d,m)
    ell<-0
    p<-c(rep(1, m+1)/(m+1),0)
    method=0
    ## Call C mable_po_m
    res<-.C("mable_po_m",
      as.double(g), as.double(p), as.integer(dm), as.double(x), as.double(y),
      as.double(y2), as.integer(N), as.double(x0), as.double(ell),
      as.double(ddell), as.double(Eps), as.integer(MaxIt), as.logical(progress),
      as.integer(conv), as.double(del), as.double(eta), as.logical(eta.known), as.integer(method))
    llik<-res[[9]][1]
    x0<-res[[8]]
    gama<-res[[1]]
    egx0<-exp(sum(gama*x0))
    Sig <- -n*matrix(res[[10]], nrow=d, ncol=d)
    se<-sqrt(diag(Sig)/n)
    z<-gama/se 
    ans<-list(m=m, mloglik=llik-n0*log(tau),  p=res[[2]], x0=x0, egx0=egx0, coefficients=gama,
      tau.n=b, tau=tau, SE=se, z=z, xNames=Dta$xNames, 
      convergence=res[[14]], delta=res[[15]], eta=eta)
  }
  else{
    lk<-rep(0, k+1)
    lr<-rep(0, k)
    pval<-rep(0,k+2)
    chpts<-rep(0,k+1)
    p<-rep(0, M[2]+2)   
    #cat("x0=",x0,"g=",g,"\n")
    #cat("x=",x[1:3],"\n")
    ## Call C mable_po
    res<-.C("mable_po",
      as.integer(M), as.double(g), as.integer(dm), as.double(p), #as.double(pi0),
      as.double(x), as.double(y), as.double(y2), as.integer(N), as.double(x0),
      as.double(lk), as.double(lr), as.double(ddell), as.double(Eps),
      as.integer(MaxIt), as.logical(progress), as.double(level), as.double(pval), 
      as.integer(chpts), as.integer(conv), as.double(eta), as.logical(eta.known))
#cat("ok4\n")
    M<-res[[1]]
    k<-M[2]-M[1]
    gama<-res[[2]][1:d]
    m<-res[[3]][2]
    x0<-res[[9]]
    lk<-res[[10]][1:(k+1)]
    lr<-res[[11]][1:k];
    llik<-lk[m-M[1]+1]
#cat("ok5\n")
    Sig <- -n*matrix(res[[12]], nrow=d, ncol=d)
#cat("ok6\n")
    se<-sqrt(diag(Sig)/n) 
    z<-gama/se 
#cat("ok7\n")
    mp1<-m+1
    egx0<-exp(sum(gama*x0))
#cat("ok8\n")
    ans<-list(M=M, lk=lk-n0*log(tau), lr=lr, m=m, tau.n=b, tau=tau, 
      mloglik=llik-n0*log(tau), p=res[[4]][1:mp1], x0=x0, egx0=egx0, 
      coefficients=gama, SE=se, z=z, xNames=Dta$xNames,
      pval=res[[17]][1:(k+1)], delta=res[[17]][k+2], 
      chpts=res[[18]][1:(k+1)]+M[1], convergence=res[[19]], eta=eta)
  }  
  ans$model<-"po"
  ans$callText<-fmla
  ans$data.name<-data.name
  ans$allvars<-allvars
  class(ans)<-"mable_reg"
  return(ans) 
}   


##################################################################
#  Select optimal degree m with given gamma 
##################################################################
#' Mable fit of the PO model with given regression coefficients
#' @param formula regression formula. Response must be \code{cbind}. See 'Details'.
#' @param data a data frame containing variables in \code{formula}.
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m} 
#'  or \code{m0 = m1 = m},  then \code{m} is a preselected degree. 
#'  If \code{m0 < m1} it specifies the set of consective candidate model degrees
#'    \code{m0:m1} for searching an optimal degree, where \code{m1-m0>3}.
#' @param g the given \eqn{d}-vector of regression coefficients
#' @param tau right endpoint of support \eqn{[0, \tau]} must be greater than 
#'      or equal to the maximum observed time
#' @param x0 a data frame specifying working baseline covariates on the right-hand-side of \code{formula}. See 'Details'.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate profile likelihood estimation of Bernstein
#'  polynomial model in proportional odds rate regression  based on interal
#'  censored event time data with given regression coefficients and select 
#'  an optimal degree m if coefficients are efficient
#'  estimates provided by other semiparametric methods. 
#' @details
#' Consider Generalized PO model with covariate for interval-censored failure time data:
#' \eqn{S(t|x) = S(t|x_0)^{\exp(\gamma'(x-x_0))}}, where \eqn{x_0} satisfies  
#'    \eqn{\gamma'(x-x_0)\ge 0}, where \eqn{x} and \eqn{x_0} may
#' contain dummy variables and interaction terms.  The working baseline \code{x0} in arguments
#' contains only the values of terms excluding dummy variables and interaction terms 
#' in the right-hand-side of \code{formula}. Thus \code{g} is the initial guess of 
#' the coefficients \eqn{\gamma} of \eqn{x-x_0} and could be longer than \code{x0}.
#'  Let \eqn{f(t|x)} and 
#'  \eqn{F(t|x) = 1-S(t|x)} be the density and cumulative distribution
#' functions of the event time given \eqn{X = x}, respectively.
#' Then \eqn{f(t|x_0)} on \eqn{[0,\tau_n]} can be approximated by
#' \eqn{f_m(t|x_0; p) = \tau^{-1}\sum_{i=0}^m p_i\beta_{mi}(t/\tau)},
#' where \eqn{p_i \ge 0}, \eqn{i = 0, \ldots, m}, \eqn{\sum_{i=0}^mp_i = 1},
#' \eqn{\beta_{mi}(u)} is the beta denity with shapes \eqn{i+1} and \eqn{m-i+1},
#'  and \eqn{\tau} is the right endpoint of support interval. So we can approximate  
#'    \eqn{S(t|x_0)} on \eqn{[0,\tau]} by
#' \eqn{S_m(t|x_0; p) = \sum_{i=0}^{m} p_i \bar B_{mi}(t/\tau)}, where
#' \eqn{\bar B_{mi}(u)}, \eqn{i = 0, \ldots, m}, is the beta survival function 
#'  with shapes \eqn{i+1} and \eqn{m-i+1}. 
#'
#' Response variable should be of the form \code{cbind(l, u)}, where  
#' \code{(l, u)} is the interval containing the event time. Data are  
#' uncensored if \code{l = u}, right censored if \code{u = Inf} or   
#' \code{u = NA}, and left censored data if \code{l = 0}. The associated  
#' covariate contains \eqn{d} columns. The baseline \code{x0} should chosen 
#' so that \eqn{\gamma'(x-x_0)} is nonnegative for all the observed \eqn{x}.
#'
#' The search for optimal degree \code{m} stops if either \code{m1} 
#' is reached or the test for change-point results in a p-value \code{pval} 
#' smaller than \code{sig.level}.
#' @return a class '\code{mable_reg}' object, a list with components
#' \itemize{
#' \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last 
#'      candidate degree when the search stoped
#' \item \code{m} the selected optimal degree \code{m}
#' \item \code{p} the estimate of \eqn{p = (p_0, \dots, p_m,p_{m+1})}, 
#'       the coefficients of Bernstein polynomial of degree \code{m}
#' \item \code{coefficients} the given regression coefficients of the PH model
#' \item \code{mloglik} the maximum log-likelihood at an optimal degree \code{m}
#' \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#' \item \code{lr} likelihood ratios for change-points evaluated at 
#'    \eqn{m \in \{m_0+1, \ldots, m_1\}}
#' \item \code{tau.n} maximum observed time \eqn{\tau_n}
#' \item \code{tau} right endpoint of support \eqn{[0, \tau)}
#' \item \code{x0} the working baseline covariates
#' \item \code{egx0} the value of \eqn{e^{\gamma'x_0}}
#' \item \code{convergence} an integer code, 0 indicates successful 
#'    completion(the iteration is convergent), 1 indicates that 
#'   the maximum candidate degree had been reached in the calculation;
#' \item \code{delta} the final convergence criterion for EM iteration;
#' \item \code{chpts} the change-points among the candidate degrees;
#' \item \code{pom} the p-value of the selected optimal degree \code{m} 
#'       as a change-point;
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references
#' Guan, Z. et al. (???) Maximum Approximate Bernstein Likelihood Estimation in 
#'   Generalized Proportional Odds Regression Model for Interval-Censored Data 
#' @examples
#' \donttest{
#' ## Simulated Weibull data
#' require(icenReg)
#' set.seed(111)
#' simdata<-simIC_weib(100, model = "po", inspections = 2, 
#'    inspectLength = 2.5, prob_cen=1)
#' sp<-ic_sp(cbind(l, u) ~ x1 + x2, data = simdata, model="po") 
#' gt<--sp$coefficients
#' res0<-maple.po(cbind(l, u) ~ x1 + x2, data = simdata, M=c(1,20), g=gt, tau=6)
#' op<-par(mfrow=c(1,2))
#' plot(res0,  which=c("likelihood","change-point"))
#' par(op)
#' res1<-mable.po(cbind(l, u) ~ x1 + x2, data = simdata, M=c(1,20),
#'    g=gt, tau=6, x0=data.frame(x1=max(simdata$x1),x2=-1))
#' op<-par(mfrow=c(2,2))    
#' plot(res1,  which=c("likelihood","change-point")) 
#' plot(res0, y=data.frame(x1=0,x2=0), which="density", add=FALSE, type="l",
#'     xlab="Time", main="Desnity Function")
#' plot(res1, y=data.frame(x1=0,x2=0), which="density", add=TRUE, lty=2, col=4)
#' lines(xx<-seq(0, 7, len=512), dweibull(xx, 2,2), lty=3, col=2, lwd=1.5)
#' legend("topright", bty="n", lty=1:3, col=c(1,4,2), c(expression(hat(f)[0]),
#'     expression(tilde(f)[0]), expression(f[0])))
#' plot(res0, y=data.frame(x1=0,x2=0), which="survival", add=FALSE, type="l",
#'     xlab="Time", main="Survival Function")
#' plot(res1, y=data.frame(x1=0,x2=0), which="survival", add=TRUE, lty=2, col=4)
#' lines(xx, 1-pweibull(xx, 2, 2), lty=2, col=2)
#' legend("topright", bty="n", lty=1:3, col=c(1,4,2), c(expression(hat(S)[0]),
#'     expression(tilde(S)[0]), expression(S[0])))
#' par(op)
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept proportional odds model
#' @concept interval censoring
#' @seealso \code{\link{mable.po}}
#' @export
maple.po<-function(formula, data, M, g, tau, x0=NULL,  controls = mable.ctrl(), 
progress=TRUE){
  if(missing(g)) stop("missing argument 'g'.")
  d<-length(g)
  data.name<-deparse(substitute(data))
  fmla<-Reduce(paste, deparse(formula))
  Dta<-get.mableData(formula, data)
  x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
  delta<-Dta$delta; n<-length(y)
#  vars<-get.facMatrix(formula, data)
#  allvars<-vars$allvars  # all variables 
  allvars<-all.vars(formula)  # all variables 
  #if(d!=length(allvars)-2) stop("Invalid argument 'g'.")
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
    oxg<-order((x%*%g)[,1])
    x0<-x[oxg[1],]   
  }
  b<-max(y2[y2<Inf], y);
  if(b>tau) stop("tau must be greater than or equal to the maximum observed time")
  if(tau==Inf || tau<0) stop("'tau' must be finite and positive.")
  y2[y2==Inf]<-tau
  y<-y/tau; y2<-y2/tau
  x<-as.matrix(x)
  if(d!=ncol(x)) stop("Invalid argument 'g'.")
#  n<-length(y)
  n0<-sum(delta==0)
  n1<-sum(delta==1)
  N<-c(n0,n1)
  ord<-order(delta)
  x<-as.matrix(x[ord,]); y<-y[ord]; y2<-y2[ord]
  if(missing(M) || length(M)==0) stop("'M' is missing.\n")
  else if(length(M)==1) M<-c(M,M)
  else if(length(M)>=2) M<-c(min(M), max(M))
  k<-M[2]-M[1]
  if(k>0 && k<=3) stop("Too few candidate model degrees.")
  k<-M[2]-M[1];
  lk<-rep(0, k+1)
  lr<-rep(0, k)
  pval<-rep(0,k+1)
  level<-controls$sig.level
  chpts<-rep(0,k+1)
  dm<-c(d,0)
  conv<-0
  del<-c(0,0)
  p<-rep(0, M[2]+2)
  ddell<-diag(0,d)
  eta<-1
  ## Call C mable_po_gamma
  res<-.C("mable_po_gamma",
    as.integer(M), as.double(g), as.integer(dm), as.double(x),
    as.double(y), as.double(y2), as.integer(N), as.double(x0), as.double(lk),
    as.double(lr), as.double(p), as.double(ddell), as.double(controls$eps.em),
    as.integer(controls$maxit.em), as.logical(progress), as.double(level),
    as.double(pval), as.integer(chpts), as.integer(conv), as.double(del), as.double(eta))
  M<-res[[1]]
  k<-M[2]-M[1]
  lk<-res[[9]][1:(k+1)]
  x0<-res[[8]]
  egx0<-exp(sum(g*x0))
  m<-res[[3]][2]
  llik<-lk[m-M[1]+1]
  mp1<-m+1
  ans<-list(m=m, mloglik=llik-n0*log(b), tau.n=b, tau=tau, p=res[[11]][1:mp1], coefficients=g,
  x0=x0, egx0=egx0, convergence=res[[19]],delta=res[[20]][1], xNames=Dta$xNames, eta=eta)
  if(k>0){
    ans$M<-M; ans$lk<-lk; ans$lr<-res[[10]][1:k]; ans$pval<-res[[17]][1:(k+1)];
    ans$chpts<-res[[18]][1:(k+1)]+M[1]; ans$pom<-res[[20]][2];}
  ans$model<-"po"
  ans$callText<-fmla
  ans$data.name<-data.name
  ans$allvars<-allvars
  class(ans)<-"mable_reg"
  return(ans)
}
##############################################
# Generalized PO model with Weibull baseline
##############################################
#' Generalized PO model with Weibull baseline
#' @param formula regression formula. Response must be \code{cbind}. See 'Details'.
#' @param data a dataset
#' @param g initial \eqn{d}-vector of regression coefficients
#' @param scale initial guess of the scale parameter for Weibull baseline
#' @param shape initial guess of the shape parameter for Weibull baseline
#' @param eta the given positive value of \eqn{\eta}. See 'Details'.
#' @param eta.known logical. If \code{TRUE} \code{eta} is the known values of \eqn{\eta},
#'     else \code{eta} is an initial guess of \eqn{\eta}. See 'Details'.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#'
#' @description Maximum likelihood estimation in generalized proportional odds 
#'  rate regression model with Weibull baseline based on interal
#'  censored event time data 
#' @details ???
#' @return a class '\code{mable_reg}' object, a list with components
#' \itemize{
#' \item \code{convergence} an integer code, 0 indicates successful 
#'    completion(the iteration is convergent), 1 indicates that 
#'   the maximum iteration had been reached in the calculation;
#' \item \code{delta} the final convergence criterion for Newton iteration;
#' }
#' @examples
#' \donttest{
#' ## Simulated Weibull data
#' require(icenReg)
#' set.seed(111)
#' simdata<-simIC_weib(100, model = "po", inspections = 2, 
#'    inspectLength = 2.5, prob_cen=1)
#' sp<-ic_sp(cbind(l, u) ~ x1 + x2, data = simdata, model="po") 
#' gt<--sp$coefficients
#' res0<-maple.po(cbind(l, u) ~ x1 + x2, data = simdata, M=c(1,20), g=gt, tau=6)
#' op<-par(mfrow=c(1,2))
#' plot(res0,  which=c("likelihood","change-point"))
#' par(op)
#' res1<-mable.po(cbind(l, u) ~ x1 + x2, data = simdata, M=c(1,20), g=gt, 
#'    tau=6, x0=data.frame(x1=max(simdata$x1),x2=-1))
#' res2<-weib.gpo(cbind(l, u) ~ x1 + x2, data = simdata, g=gt, scale=2, shape=2)  

#' op<-par(mfrow=c(2,2))    
#' plot(res1,  which=c("likelihood","change-point")) 
#' plot(res0, y=data.frame(x1=0,x2=0), which="density", add=FALSE, type="l",
#'     xlab="Time", main="Desnity Function")
#' plot(res1, y=data.frame(x1=0,x2=0), which="density", add=TRUE, lty=2, col=4)
#' lines(xx<-seq(0, 7, len=512), dweibull(xx, 2,2), lty=3, col=2, lwd=1.5)
#' lines(xx, dweibull(xx, res2$shape, res2$scale), lty=5, col=5, lwd=1.5)
#' legend("topright", bty="n", lty=1:3, col=c(1,4,2), c(expression(hat(f)[0]),
#'     expression(tilde(f)[0]), expression(f[0])))
#' plot(res0, y=data.frame(x1=0,x2=0), which="survival", add=FALSE, type="l",
#'     xlab="Time", main="Survival Function")
#' plot(res1, y=data.frame(x1=0,x2=0), which="survival", add=TRUE, lty=2, col=4)
#' lines(xx, 1-pweibull(xx, 2, 2), lty=2, col=2)
#' lines(xx, 1-pweibull(xx, res2$shape, res2$scale), lty=5, col=5, lwd=1.5)
#' legend("topright", bty="n", lty=1:3, col=c(1,4,2), c(expression(hat(S)[0]),
#'     expression(tilde(S)[0]), expression(S[0])))
#' par(op)
#' }
#' @export
weib.gpo<-function(formula, data, g, scale, shape, eta=1, eta.known=TRUE, controls = mable.ctrl(), 
        progress=TRUE){
  data.name<-deparse(substitute(data))
  fmla<-Reduce(paste, deparse(formula))
  vars<-get.facMatrix(formula, data)
  allvars<-vars$allvars  # all variables 
  Dta<-get.mableData(formula, data)
  x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
  delta<-Dta$delta 
  delta[abs(y-y2) < 1e-6]<-0
  x<-as.matrix(x)
  d<-ncol(x)
  if(missing(g)) rep(0,d)
  if(d!=length(g)) stop("Invalid argument 'g'.")
  n<-length(y)
  n0<-sum(delta==0)
  n1<-sum(delta==1)
  N<-c(n0,n1)
  ord<-order(delta)
  x<-as.matrix(x[ord,]); y<-y[ord]; y2<-y2[ord] 
  y2[y2==Inf]<--1 # right censored
  lk=0
  np<-d+2+(!eta.known)
  ddell<-diag(0,np)
  if(eta.known) theta<-c(g, scale, shape, eta)
  else theta<-c(g, eta, scale, shape)###???
  conv<-0
  del<-10
  res<-.C("weib_gpo",
    as.double(theta), as.integer(d), as.double(x), as.integer(n0), as.integer(n1), 
    as.double(y), as.double(y2), as.double(lk), as.double(ddell), 
    as.double(controls$eps), as.integer(controls$maxit), as.logical(progress), 
    as.integer(conv), as.double(del), as.logical(eta.known))
  Sig <- -n*matrix(res[[9]], nrow=np, ncol=np)
  se<-sqrt(diag(Sig)/n) 
  if(eta.known) z<-res[[1]][-(d+3)]/se 
  else z<-res[[1]]/se
  ans<-list(mloglik=res[[8]][1], coefficients=res[[1]][1:d], scale=res[[1]][d+2],
    shape=res[[1]][d+3], z=z, SE=se,  convergence=res[[13]][1], delta=res[[14]][1], 
    xNames=Dta$xNames)  
  if(eta.known) ans$eta=eta
  else ans$eta=res[[1]][d+1] 
  ans$model<-"po"
  ans$callText<-fmla
  ans$data.name<-data.name
#  ans$vars<-vars
#  ans$fmatrices<-vars[[1]]
#  ans$factors<-vars$factors
  ans$allvars<-allvars
#  ans$whichisfactor<-vars$whichisfactor
 class(ans)<-"mable_reg"
  return(ans)
}  

