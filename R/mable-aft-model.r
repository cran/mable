##########################################################
#   MABLE of AFT model based on interval-censored data   #
##########################################################
#     AFT model with covariate for interval-censored     #
#  failure time data: S(t|x)=S(t/exp(gamma'(x-x0))|x0)   #
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
#' @param formula regression formula. Response must be \code{cbind}.  See 'Details'.
#' @param data a dataset
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m} or \code{m0 = m1 = m},  
#'   then \code{m} is a preselected degree. If \code{m0 < m1} it specifies the set of 
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.  
#' @param g initial guess of \eqn{d}-vector of regression coefficients.  See 'Details'. 
#' @param tau a finite truncation time greater than the maximum observed time \eqn{\tau}. See 'Details'.
#' @param x0 a working baseline covariate \eqn{x_0}. See 'Details'. 
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @details
#' Consider the accelerated failure time model with covariate for interval-censored failure time data: 
#' \eqn{S(t|x) = S(t \exp(\gamma'(x-x_0))|x_0)}, where \eqn{x_0} is a baseline covariate.   
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
#' \eqn{S(t|x)=S(t \exp(\gamma'(x-x_0))|x_0)} on \eqn{[\tau, \infty)} is negligible for
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
#'   \item \code{egx0} the value of \eqn{e^{\gamma'x_0}} 
#'   \item \code{convergence} an integer code, 1 indicates either the EM-like iteration for finding 
#'     maximum likelihood reached the maximum iteration for at least one \code{m} or the search of 
#'     an optimal degree using change-point method reached the maximum candidate degree,
#'     2 indicates both occured, and 0 indicates a successful completion.  
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
#' @author Zhong Guan <zguan@iusb.edu>
#' @references 
#' Guan, Z. (2019) Maximum Approximate Likelihood Estimation in Accelerated Failure Time Model for Interval-Censored Data, 
#' arXiv:1911.07087.
#' @examples \donttest{
#' ## Breast Cosmesis Data
#'   require(coxinterval) 
#'   bcos=cosmesis
#'   bcos2<-data.frame(bcos[,1:2], x=1*(bcos$treat=="RCT"))
#'   g <- 0.41 #Hanson and  Johnson 2004, JCGS
#'   aft.res<-mable.aft(cbind(left, right)~x, data=bcos2, M=c(1, 30), g=g, tau=100, x0=1)
#'   op<-par(mfrow=c(1,2), lwd=1.5)
#'   plot(x=aft.res, which="likelihood")
#'   plot(x=aft.res, y=data.frame(x=0), which="survival", model='aft', type="l", col=1, 
#'       add=FALSE, main="Survival Function")
#'   plot(x=aft.res, y=data.frame(x=1), which="survival", model='aft', lty=2, col=1)
#'   legend("bottomleft", bty="n", lty=1:2, col=1, c("Radiation Only", "Radiation and Chemotherapy"))
#'   par(op)
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Accelerated failure time model 
#' @concept interval censoring
#' @seealso \code{\link{maple.aft}}
#' @importFrom stats coef reformulate terms
#' @importFrom survival Surv
#' @export
mable.aft<-function(formula, data, M, g=NULL, tau=1, x0=NULL,    
                   controls = mable.ctrl(), progress=TRUE){
    data.name<-deparse(substitute(data)) 
    fmla<-Reduce(paste, deparse(formula))
    Dta<-get.mableData(formula, data)
    x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
    if(is.null(g)){
        stop("argument 'g' is missing.")
        #status<-1*(y2<Inf)
        #Y<-y 
        #Y[y2<Inf]<-(y[y2<Inf]+y2[y2<Inf])/2
        #rtc.data<-data.frame(Y = Y, status = status, x=x)
        #fmla<-reformulate(attr(terms(formula), "term.labels"), response="Surv(Y, status)")
        #rkest<-aftsrr(fmla, data = rtc.data)
        #g<--as.numeric(coef(rkest))      
    }
    delta<-Dta$delta
    b<-max(y2[y2<Inf], y);
    if(b>=tau) stop("tau must be greater than the maximum observed time")
    y<-y/tau; y2<-y2/tau
#    y2[y2==Inf]<-.Machine$double.xmax/2 
    y2[y2==Inf]<-1 # truncation interval is [0, tau)
    x<-as.matrix(x)
    d<-length(g)
    if(d!=ncol(x)) stop("length of gamma does not match number of covariates.")
    n<-length(y)
    n0<-sum(delta==0)
    n1<-sum(delta==1)
    n<-n0+n1
    N<-c(n0,n1)
    dm<-c(d,0)
    conv<-0
    del<-0
    ord<-order(delta)
    x<-as.matrix(x[ord,]); y<-y[ord]; y2<-y2[ord]
    if(is.null(x0)) x0<-rep(0,d)
    Eps<-c(controls$eps, controls$eps.em)
    MaxIt<-c(controls$maxit, controls$maxit.em)
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
    ans$tau<-tau
    ans$xNames<-Dta$xNames
    if(k==0){
        m<-M[1]    
        dm<-c(d,m)
        ell<-0
        p<-rep(1,m+1)/(m+1)
        ## Call C mable_aft_m
        res<-.C("mable_aft_m",
            as.double(-g), as.double(p), as.integer(dm), as.double(x), as.double(y),  
            as.double(y2), as.integer(N), as.double(x0), as.double(ell), 
            as.double(ddell), as.double(Eps), as.integer(MaxIt), 
            as.logical(progress), as.integer(conv), as.double(del))
        ans$m<-m
        ans$mloglik<-res[[9]][1]
        ans$p<-res[[2]]
        ans$x0<-res[[8]]
        ans$coefficients<-res[[1]]
        ans$egx0<-exp(sum(res[[1]]*res[[8]]))
        Sig<--n*matrix(res[[10]], nrow=d, ncol=d)
        ans$SE<-sqrt(diag(Sig)/n)
        ans$z<-res[[1]]/ans$SE 
        ans$convergence<-res[[14]]
        ans$delta<-res[[15]] 
    }
    else{
        lk<-rep(0, k+1)
        lr<-rep(0, k)    
        pval<-rep(0,k+1)
        chpts<-rep(0,k+1)
        level<-controls$sig.level
        p<-rep(0, M[2]+1)
        ## Call C mable_aft
        res<-.C("mable_aft",
            as.integer(M), as.double(-g), as.integer(dm), as.double(p), as.double(x),  
            as.double(y), as.double(y2), as.integer(N), as.double(x0), as.double(lk), 
            as.double(lr), as.double(ddell), as.double(Eps), as.integer(MaxIt), 
            as.logical(progress), as.double(pval), as.integer(chpts), as.double(level), 
            as.integer(conv))
        M<-res[[1]]
        ans$M<-M
        k<-M[2]-M[1]
        ans$coefficients<-res[[2]]
        ans$x0<-res[[9]]
        ans$egx0<-exp(sum(res[[2]]*res[[9]]))
        Sig<--n*matrix(res[[12]], nrow=d, ncol=d)
        ans$SE<-sqrt(diag(Sig)/n)
        ans$lk<-res[[10]][1:(k+1)]-n0*log(b) 
        ans$lr<-res[[11]][1:k]
        ans$m<-res[[3]][2]
        ans$mloglik<-res[[10]][ans$m-M[1]+1]-n0*log(b)
        mp1<-ans$m+1
        ans$p<-res[[4]][1:mp1]  
        ans$z<-res[[2]]/ans$SE
        ans$pval<-res[[16]][1:(k+1)]
        ans$chpts<-res[[17]][1:(k+1)]+M[1]
        ans$convergence<-res[[19]]
        ans$delta<-res[[16]][k+1]   
    }
    ans$coefficients<--ans$coefficients
    ans$egx0<-1/ans$egx0
    ans$model<-"aft"
    ans$callText<-fmla
    ans$data.name<-data.name
    class(ans)<-"mable_reg"
    return(ans)
}
####################################################################################
#  Maximum Approximate Profile Likelihood Estimation in AFT model with a given gamma
# M: set of positive integers as candidate degrees of Bernstein poly model
##################################################################
#' Mable fit of AFT model with given regression coefficients for AFT model
#' @param formula regression formula. Response must be \code{cbind}.  See 'Details'.
#' @param data a dataset
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m} or \code{m0 = m1 = m},  
#'   then \code{m} is a preselected degree. If \code{m0 < m1} it specifies the set of 
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0 > 3}.  
#' @param g the given \eqn{d}-vector of regression coefficients 
#' @param tau a truncation time greater than or equal to the maximum observed time \eqn{\tau}. See 'Details'. 
#' @param x0 a working baseline covariate \eqn{x_0}. See 'Details'. 
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate profile likelihood estimation of Bernstein
#'  polynomial model in accelerated failure time based on interal 
#'  censored event time data with given regression coefficients which are efficient
#'  estimates provided by other semiparametric methods. 
#' @details
#' Consider the accelerated failure time model with covariate for interval-censored failure time data: 
#' \eqn{S(t|x) = S(t \exp(\gamma'(x-x_0))|x_0)}, where \eqn{x_0} is a baseline covariate.   
#'   Let \eqn{f(t|x)} and \eqn{F(t|x) = 1-S(t|x)} be the density and cumulative distribution
#' functions of the event time given \eqn{X = x}, respectively.
#' Then \eqn{f(t|x_0)} on a truncation interval \eqn{[0, \tau]} can be approximated by  
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
#' \eqn{S(t|x) = S(t \exp(\gamma'(x-x_0))|x_0)} on \eqn{[\tau, \infty)} is negligible for
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
#'   \item \code{egx0} the value of \eqn{e^{\gamma'x_0}} 
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
#' @author Zhong Guan <zguan@iusb.edu>
#' @references 
#' Guan, Z. (2019) Maximum Approximate Likelihood Estimation in Accelerated Failure Time Model for Interval-Censored Data, 
#' arXiv:1911.07087.
#' @examples \donttest{
#' ## Breast Cosmesis Data
#'   require(coxinterval) 
#'   bcos=cosmesis
#'   bcos2<-data.frame(bcos[,1:2], x=1*(bcos$treat=="RCT"))
#'   g<-0.41 #Hanson and  Johnson 2004, JCGS, 
#'   res1<-maple.aft(cbind(left, right)~x, data=bcos2, M=c(1,30),  g=g, tau=100, x0=1)
#'   op<-par(mfrow=c(1,2), lwd=1.5)
#'   plot(x=res1, which="likelihood")
#'   plot(x=res1, y=data.frame(x=0), which="survival", model='aft', type="l", col=1, 
#'       add=FALSE, main="Survival Function")
#'   plot(x=res1, y=data.frame(x=1), which="survival", model='aft', lty=2, col=1)
#'   legend("bottomleft", bty="n", lty=1:2, col=1, c("Radiation Only", "Radiation and Chemotherapy"))
#'   par(op)
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Accelerated failure time model 
#' @concept interval censoring
#' @seealso \code{\link{mable.aft}} 
#' @export
maple.aft<-function(formula, data, M, g, tau=1, x0=NULL, 
            controls = mable.ctrl(), progress=TRUE){
    data.name<-deparse(substitute(data)) 
    fmla<-Reduce(paste, deparse(formula))
    Dta<-get.mableData(formula, data)
    x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
    delta<-Dta$delta
    b<-max(y2[y2<Inf], y);
    if(b>tau) stop("tau must be greater than or equal to the maximum observed time")
    y<-y/tau; y2<-y2/tau
#    y2[y2==Inf]<-.Machine$double.xmax/2 
    y2[y2==Inf]<-1 # truncation interval is [0, tau)
    x<-as.matrix(x)
    d<-length(g)
    if(d!=ncol(x)) stop("length of gamma and number of covariates do not match.")
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
    p<-rep(0, M[1]+k+1) 
    if(is.null(x0)) x0<-rep(0,d)
    ddell<-diag(0,d)
    pval<-rep(0,k+1)
    chpts<-rep(0,k+1)
    level<-controls$sig.level
    dm<-c(d,0)
    conv<-0
    del<-0
    ## Call C mable_aft_gamma
    res<-.C("mable_aft_gamma",
        as.integer(M), as.double(-g), as.integer(dm), as.double(x), as.double(y),  
        as.double(y2), as.integer(N), as.double(x0), as.double(lk), as.double(lr), 
        as.double(p), as.double(ddell), as.double(controls$eps), as.integer(controls$maxit), 
        as.logical(progress), as.double(pval), as.integer(chpts), as.double(level), 
        as.integer(conv), as.double(del))
    ans<-list()
    M<-res[[1]]
    ans$x0<-res[[8]]
    ans$egx0<-exp(sum(g*x0))
    ans$m<-res[[3]][2]
    ans$tau.n<-b
    ans$tau<-tau
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
    class(ans)<-"mable_reg"
    return(ans)
}

