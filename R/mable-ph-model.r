############################################################
#                                                          #
#    MABLE for Cox Proportional Hazard Regression Model    #
#                                                          #
############################################################
#  Cox's PH model with covariate for interval-censored     #
#  failure time data: S(t|x)=S(t|x0)^exp(gamma'(x-x0))     #
############################################################
# R CMD SHLIB mable-ph-model.c
# R --arch x64 CMD SHLIB mable-ph-model.c
# setwd("C:\\Users\\zguan\\Documents\\papers\\bernstein polynomials\\survival function\\C")
# dyn.load("mable-ph-model")
#
##################################################################
#' Mable fit of Cox's proportional hazards regression model 
#' @param formula regression formula. Response must be \code{cbind}.  See 'Details'.
#' @param data a dataset
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m} or \code{m0 = m1 = m},  
#'   then \code{m} is a preselected degree. If \code{m0<m1} it specifies the set of 
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.  
#' @param g initial guess of \eqn{d}-vector of regression coefficients.  See 'Details'. 
#' @param pi0 Initial guess of \eqn{\pi(x_0) = F(\tau_n|x_0)}. Without right censored data, \code{pi0 = 1}. See 'Details'.
#' @param tau right endpoint of support \eqn{[0, \tau)} must be greater than or equal to the maximum observed time
#' @param x0 a working baseline covariate. See 'Details'. 
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate Bernstein/Beta likelihood estimation in Cox's proportional hazards regression model 
#' based on interal censored event time data.
#' @details
#' Consider Cox's PH model with covariate for interval-censored failure time data: 
#' \eqn{S(t|x) = S(t|x_0)^{\exp(\gamma'(x-x_0))}}, where \eqn{x_0} satisfies \eqn{\gamma'(x-x_0)\ge 0}.   
#'   Let \eqn{f(t|x)} and \eqn{F(t|x) = 1-S(t|x)} be the density and cumulative distribution
#' functions of the event time given \eqn{X = x}, respectively.
#' Then \eqn{f(t|x_0)} on \eqn{[0, \tau_n]} can be approximated by  
#' \eqn{f_m(t|x_0, p) = \tau_n^{-1}\sum_{i=0}^m p_i\beta_{mi}(t/\tau_n)},
#' where \eqn{p_i \ge 0}, \eqn{i = 0, \ldots, m}, \eqn{\sum_{i=0}^mp_i = 1-p_{m+1}},  
#' \eqn{\beta_{mi}(u)} is the beta denity with shapes \eqn{i+1} and \eqn{m-i+1}, and 
#' \eqn{\tau_n} is the largest observed time, either uncensored time, or right endpoint of interval/left censored,
#' or left endpoint of right censored time. So we can approximate  \eqn{S(t|x_0)} on \eqn{[0, \tau_n]} by
#' \eqn{S_m(t|x_0; p) = \sum_{i=0}^{m+1} p_i \bar B_{mi}(t/\tau_n)}, where 
#' \eqn{\bar B_{mi}(u)}, \eqn{i = 0, \ldots, m}, is the beta survival function with shapes 
#'  \eqn{i+1} and \eqn{m-i+1}, \eqn{\bar B_{m,m+1}(t) = 1}, \eqn{p_{m+1} = 1-\pi(x_0)}, and
#' \eqn{\pi(x_0) = F(\tau_n|x_0)}. For data without right-censored time, \eqn{p_{m+1} = 1-\pi(x_0) = 0}.
#'
#' Response variable should be of the form \code{cbind(l, u)}, where \code{(l, u)} is the interval 
#' containing the event time. Data is uncensored if \code{l = u}, right censored 
#' if \code{u = Inf} or \code{u = NA}, and  left censored data if \code{l = 0}.
#' The associated covariate contains \eqn{d} columns. The baseline \code{x0} should chosen so that 
#' \eqn{\gamma'(x-x_0)} is nonnegative for all the observed \eqn{x} and 
#' all \eqn{\gamma} in a neighborhood of its true value.
#'
#' A missing initial value of \code{g} is imputed by \code{ic_sp()} of package \code{icenReg}. 
#'
#'  The search for optimal degree \code{m} stops if either \code{m1} is reached or the test 
#'  for change-point results in a p-value \code{pval} smaller than \code{sig.level}.
#' This process takes longer than \code{\link{maple.ph}} to select an optimal degree.  
#' @return A list with components
#' \itemize{ 
#'   \item \code{m} the selected/preselected optimal degree \code{m}
#'   \item \code{p} the estimate of \eqn{p = (p_0, \dots, p_m, p_{m+1})}, the coefficients of Bernstein polynomial of degree \code{m}
#'   \item \code{coefficients} the estimated regression coefficients of the PH model
#'   \item \code{SE} the standard errors of the estimated regression coefficients 
#'   \item \code{z} the z-scores of the estimated regression coefficients 
#'   \item \code{mloglik} the maximum log-likelihood at an optimal degree \code{m}
#'   \item \code{tau.n} maximum observed time \eqn{\tau_n}
#'   \item \code{tau} right endpoint of support \eqn{[0, \tau)}
#'   \item \code{x0} the working baseline covariates 
#'   \item \code{egx0} the value of \eqn{e^{\gamma'x_0}} 
#'   \item \code{convergence} an integer code, 1 indicates either the EM-like 
#'     iteration for finding maximum likelihood reached the maximum iteration for at least one \code{m} 
#'     or the search of an optimal degree using change-point method reached the maximum candidate degree,
#'     2 indicates both occured, and 0 indicates a successful completion.  
#'   \item \code{delta} the final \code{delta} if \code{m0 = m1} or the final \code{pval} of the change-point 
#'      for searching the optimal degree \code{m};
#'  }
#'  and, if \code{m0<m1},
#' \itemize{
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last candidate degree when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0,\ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iusb.edu>
#' @references 
#' Guan, Z. Maximum Approximate Bernstein Likelihood Estimation in Proportional Hazard Model for Interval-Censored Data, 
#' Statistics in Medicine. 2020; 1–21. https://doi.org/10.1002/sim.8801.
#' @examples
#' \donttest{
#'    # Ovarian Cancer Survival Data
#'    require(survival)
#'    futime2<-ovarian$futime
#'    futime2[ovarian$fustat==0]<-Inf
#'    ovarian2<-data.frame(age=ovarian$age, futime1=ovarian$futime, 
#'         futime2=futime2)
#'    ova<-mable.ph(cbind(futime1, futime2) ~ age, data = ovarian2, 
#'         M=c(2,35), g=.16, x0=35)
#'    op<-par(mfrow=c(2,2))
#'    plot(ova, which = "likelihood")
#'    plot(ova, which = "change-point")
#'    plot(ova, y=data.frame(age=60), which="survival", add=FALSE, type="l", 
#'          xlab="Days", main="Age = 60")
#'    plot(ova, y=data.frame(age=65), which="survival", add=FALSE, type="l", 
#'          xlab="Days", main="Age = 65")
#'    par(op)
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Cox proportional hazards model 
#' @concept interval censoring
#' @seealso \code{\link{maple.ph}} 
# @useDynLib mable-ph-model .registration = TRUE
#' @importFrom icenReg ic_sp
#' @export
mable.ph<-function(formula, data, M, g=NULL, pi0=NULL, tau=Inf, x0=NULL, 
        controls = mable.ctrl(), progress=TRUE){
    data.name<-deparse(substitute(data)) 
    fmla<-Reduce(paste, deparse(formula))
    Dta<-get.mableData(formula, data)
    x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
    delta<-Dta$delta
    if(is.null(g)){
        g<-ic_sp(formula, data, model = 'ph')$coefficients
    }
    if(is.null(pi0)) pi0<-mean(y2<Inf)
    b<-max(y2[y2<Inf], y) 
    if(b>tau) stop("tau must be greater than or equal to the maximum observed time")
    y<-y/b; y2<-y2/b
    if(tau==Inf) y2[y2==Inf]<-.Machine$double.xmax/2
    else y2[y2==Inf]<-tau
    if(!is.matrix(x)) x<-matrix(x, ncol=1)
    d<-length(g)
    if(d!=ncol(x)) stop("length of gama does not match number of covariates.")
    n<-length(y)
    n0<-sum(delta==0)
    n1<-sum(delta==1)
    n<-n0+n1
    N<-c(n0,n1) 
    dm<-c(d,0)
    ddell<-diag(0,d)
    conv<-0
    del<-0
    ord<-order(delta)
    x<-as.matrix(x[ord,]); y<-y[ord]; y2<-y2[ord]
    if(is.null(x0)) 
        for(i in 1:d) x0[i]<-ifelse(g[i]>=0, min(x[,i]), max(x[,i]))
    Eps<-c(controls$eps, controls$eps.em, controls$eps.nt, .Machine$double.eps)
    MaxIt<-c(controls$maxit, controls$maxit.em, controls$maxit.nt)
    level<-controls$sig.level
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
    k<-M[2]-M[1]
    if(k>0 && k<=3) stop("Too few candidate model degrees.")
    if(k==0){  
        m<-M[1]
        dm<-c(d,m)
        ell<-0  
        p<-rep(1, m+2)/(m+2)
        ## Call C mable_ph_m
        res<-.C("mable_ph_m",
            as.double(g), as.double(p), as.integer(dm), as.double(x), as.double(y),  
            as.double(y2), as.integer(N), as.double(x0), as.double(ell), 
            as.double(ddell), as.double(Eps), as.integer(MaxIt),  
            as.logical(progress), as.integer(conv), as.double(del)) 
        llik<-res[[9]][1]
        x0<-res[[8]]
        gama<-res[[1]]
        egx0<-exp(sum(gama*x0))
        Sig <- -n*matrix(res[[10]], nrow=d, ncol=d)
        se<-sqrt(diag(Sig)/n)
        ans<-list(m=m, mloglik=llik-n0*log(b),  p=res[[2]], x0=x0, egx0=egx0, coefficients=gama, 
            tau.n=b, tau=tau, SE=se, z=gama/se, xNames=Dta$xNames, convergence=res[[14]],
            delta=res[[15]])
    }
    else{
        lk<-rep(0, k+1)
        lr<-rep(0, k)
        pval<-rep(0,k+1)
        chpts<-rep(0,k+1)
        p<-rep(0, M[2]+2)
        ## Call C mable_ph
        res<-.C("mable_ph",
            as.integer(M), as.double(g), as.integer(dm), as.double(p), as.double(pi0),
            as.double(x), as.double(y), as.double(y2), as.integer(N), as.double(x0), 
            as.double(lk), as.double(lr), as.double(ddell), as.double(Eps), 
            as.integer(MaxIt), as.logical(progress), as.double(level),
            as.double(pval), as.integer(chpts), as.integer(conv))
        M<-res[[1]]
        k<-M[2]-M[1]
        gama<-res[[2]]
        #x0<-res[[10]]
        Sig <- -n*matrix(res[[13]], nrow=d, ncol=d)
        se<-sqrt(diag(Sig)/n)
        lr<-res[[12]][1:k]; 
        lk<-res[[11]][1:(k+1)]
        m<-res[[3]][2]  
        llik<-lk[m-M[1]+1]
        mp2<-m+2
        egx0<-exp(sum(gama[1:d]*x0))
        ans<-list(M=M, lk=lk-n0*log(b), lr=lr, m=m, tau.n=b, tau=tau, mloglik=llik-n0*log(b), 
            p=res[[4]][1:mp2], x0=x0, egx0=egx0,  SE=se, z=gama/se, 
            coefficients=gama[1:d], convergence=res[[20]], xNames=Dta$xNames,
            pval=res[[18]][1:(k+1)], chpts=res[[19]][1:(k+1)]+M[1], delta=res[[18]][k+1])
    } 
    ans$model<-"ph"
    ans$callText<-fmla
    ans$data.name<-data.name
    class(ans)<-"mable_reg"
    return(ans)
}
##################################################################
#  Select optimal degree m with a given gamma
##################################################################
#' Mable fit of the PH model with given regression coefficients 
#' @param formula regression formula. Response must be \code{cbind}.  See 'Details'.
#' @param data a dataset
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m} or \code{m0 = m1 = m},  
#'   then \code{m} is a preselected degree. If \code{m0 < m1} it specifies the set of 
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.  
#' @param g the given \eqn{d}-vector of regression coefficients 
#' @param pi0 Initial guess of \eqn{\pi(x_0) = F(\tau_n|x_0)}. Without right censored data, \code{pi0 = 1}. See 'Details'.
#' @param tau right endpoint of support \eqn{[0, \tau)} must be greater than or equal to the maximum observed time
#' @param x0 a working baseline covariate. See 'Details'. 
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate profile likelihood estimation of Bernstein
#'  polynomial model in Cox's proportional hazards regression  based on interal 
#'  censored event time data with given regression coefficients which are efficient
#'  estimates provided by other semiparametric methods.
#' @details
#' Consider Cox's PH model with covariate for interval-censored failure time data: 
#' \eqn{S(t|x) = S(t|x_0)^{\exp(\gamma'(x-x_0))}}, where \eqn{x_0} satisfies \eqn{\gamma'(x-x_0)\ge 0}.   
#'   Let \eqn{f(t|x)} and \eqn{F(t|x) = 1-S(t|x)} be the density and cumulative distribution
#' functions of the event time given \eqn{X = x}, respectively.
#' Then \eqn{f(t|x_0)} on \eqn{[0,\tau_n]} can be approximated by  
#' \eqn{f_m(t|x_0; p) = \tau_n^{-1}\sum_{i=0}^m p_i\beta_{mi}(t/\tau_n)},
#' where \eqn{p_i \ge 0}, \eqn{i = 0, \ldots, m}, \eqn{\sum_{i=0}^mp_i = 1-p_{m+1}},
#' \eqn{\beta_{mi}(u)} is the beta denity with shapes \eqn{i+1} and \eqn{m-i+1}, and
#' \eqn{\tau_n} is the largest observed time, either uncensored time, or right endpoint of interval/left censored,
#' or left endpoint of right censored time. So we can approximate  \eqn{S(t|x_0)} on \eqn{[0, \tau_n]} by
#' \eqn{S_m(t|x_0; p) = \sum_{i=0}^{m+1} p_i \bar B_{mi}(t/\tau_n)}, where 
#' \eqn{\bar B_{mi}(u)}, \eqn{i = 0, \ldots, m}, is the beta survival function with shapes 
#'  \eqn{i+1} and \eqn{m-i+1}, \eqn{\bar B_{m,m+1}(t) =  1}, \eqn{p_{m+1} = 1-\pi(x_0)}, and
#' \eqn{\pi(x_0) = F(\tau_n|x_0)}. For data without right-censored time, \eqn{p_{m+1} = 1-\pi(x_0) = 0.}  
#'
#' Response variable should be of the form \code{cbind(l, u)}, where \code{(l, u)} is the interval 
#' containing the event time. Data is uncensored if \code{l = u}, right censored 
#' if \code{u = Inf} or \code{u = NA}, and  left censored data if \code{l = 0}.
#' The associated covariate contains \eqn{d} columns. The baseline \code{x0} should chosen so that 
#' \eqn{\gamma'(x-x_0)} is nonnegative for all the observed \eqn{x}.
#'
#'  The search for optimal degree \code{m} stops if either \code{m1} is reached or the test 
#'  for change-point results in a p-value \code{pval} smaller than \code{sig.level}.
#' @return a class '\code{mable_reg}' object, a list with components
#' \itemize{ 
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last candidate degree when the search stoped
#'   \item \code{m} the selected optimal degree \code{m}
#'   \item \code{p} the estimate of \eqn{p = (p_0, \dots, p_m,p_{m+1})}, the coefficients of Bernstein polynomial of degree \code{m}
#'   \item \code{coefficients} the given regression coefficients of the PH model
#'   \item \code{mloglik} the maximum log-likelihood at an optimal degree \code{m}
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{tau.n} maximum observed time \eqn{\tau_n}
#'   \item \code{tau} right endpoint of support \eqn{[0, \tau)}
#'   \item \code{x0} the working baseline covariates 
#'   \item \code{egx0} the value of \eqn{e^{\gamma'x_0}} 
#'   \item \code{convergence} an integer code. 0 indicates successful completion(the iteration is 
#'    convergent). 1 indicates that the maximum candidate degree had been reached in the calculation;
#'   \item \code{delta} the final convergence criterion for EM iteration;
#'   \item \code{chpts} the change-points among the candidate degrees;
#'   \item \code{pom} the p-value of the selected optimal degree \code{m} as a change-point;
#'  }  
#' @author Zhong Guan <zguan@iusb.edu>   
#' @references 
#' Guan, Z. (2019) Maximum Approximate Bernstein Likelihood Estimation in Proportional Hazard Model for Interval-Censored Data, 
#' arXiv:1906.08882 .
#' @examples
#' \donttest{
#'  ## Simulated Weibull data
#'    require(icenReg) 
#'    set.seed(123)
#'    simdata<-simIC_weib(70, inspections = 5, inspectLength = 1)
#'    sp<-ic_sp(cbind(l, u) ~ x1 + x2, data = simdata)
#'    res0<-maple.ph(cbind(l, u) ~ x1 + x2, data = simdata, M=c(2,20), 
#'         g=sp$coefficients, tau=7)
#'    op<-par(mfrow=c(1,2))
#'    plot(res0,  which=c("likelihood","change-point"))
#'    par(op)
#'    res1<-mable.ph(cbind(l, u) ~ x1 + x2, data = simdata, M=res0$m, 
#'       g=c(.5,-.5), tau=7)
#'    op<-par(mfrow=c(1,2))
#'    plot(res1, y=data.frame(x=0, x2=0), which="density", add=FALSE, type="l", 
#'        xlab="Time", main="Desnity Function")
#'    lines(xx<-seq(0, 7, len=512), dweibull(xx, 2,2), lty=2, col=2)
#'    legend("topright", bty="n", lty=1:2, col=1:2, c("Estimated","True"))
#'    plot(res1, y=data.frame(x=0, x2=0), which="survival", add=FALSE, type="l", 
#'        xlab="Time", main="Survival Function")
#'    lines(xx, 1-pweibull(xx, 2, 2), lty=2, col=2)
#'    legend("topright", bty="n", lty=1:2, col=1:2, c("Estimated","True"))
#'    par(op)
#' }
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Cox proportional hazards model 
#' @concept interval censoring
#' @seealso \code{\link{mable.ph}}
#' @export
maple.ph<-function(formula, data, M, g, pi0=NULL, tau=Inf, x0=NULL, 
        controls = mable.ctrl(), progress=TRUE){
    data.name<-deparse(substitute(data)) 
    fmla<-Reduce(paste, deparse(formula))
    Dta<-get.mableData(formula, data)
    x<-Dta$x;  y<-Dta$y; y2<-Dta$y2
    delta<-Dta$delta
    if(is.null(pi0)) pi0<-mean(y2<Inf)
    b<-max(y2[y2<Inf], y);
    if(b>tau) stop("tau must be greater than or equal to the maximum observed time")
    y<-y/b; y2<-y2/b
    if(tau==Inf) y2[y2==Inf]<-.Machine$double.xmax/2
    else y2[y2==Inf]<-tau
    x<-as.matrix(x)
    d<-length(g)
    if(d!=ncol(x)) stop("length of g and number of covariates do not match.")
    n<-length(y)
    n0<-sum(delta==0)
    n1<-sum(delta==1)
    N<-c(n0,n1) 
    ord<-order(delta)
    x<-as.matrix(x[ord,]); y<-y[ord]; y2<-y2[ord]
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
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
    if(is.null(x0)) 
        for(i in 1:d) x0[i]<-ifelse(g[i]>=0, min(x[,i]), max(x[,i]))
    ddell<-diag(0,d)
    ## Call C mable_ph_gamma
    res<-.C("mable_ph_gamma",
        as.integer(M), as.double(g), as.integer(dm), as.double(pi0), as.double(x), 
        as.double(y), as.double(y2), as.integer(N), as.double(x0), as.double(lk), 
        as.double(lr), as.double(p), as.double(ddell), as.double(controls$eps.em),
        as.integer(controls$maxit.em), as.logical(progress), as.double(level),
        as.double(pval), as.integer(chpts), as.integer(conv), as.double(del))
    M<-res[[1]]
    k<-M[2]-M[1]
    lk<-res[[10]][1:(k+1)]
    x0<-res[[9]]
    egx0<-exp(sum(g*x0))
    m<-res[[3]][2]  
    llik<-lk[m-M[1]+1]
    mp2<-m+2
    ans<-list(m=m, mloglik=llik-n0*log(b), tau.n=b, tau=tau, p=res[[12]][1:mp2], coefficients=g, 
    x0=x0, egx0=egx0, convergence=res[[20]],delta=res[[21]][1], xNames=Dta$xNames)
    if(k>0){
        ans$M<-M; ans$lk<-lk-n0*log(b); ans$lr<-res[[11]][1:k]; ans$pval<-res[[18]][1:(k+1)];
        ans$chpts<-res[[19]][1:(k+1)]+M[1]; ans$pom<-res[[21]][2];}
    ans$model<-"ph"
    ans$callText<-fmla
    ans$data.name<-data.name
    class(ans)<-"mable_reg"
    return(ans)
}
###############################################
#  Mable fit based on interval censored data without covariate
# M: set of positive integers as candidate degrees of Bernstein poly model
#   p: coefficients of Bernstein poly approx of the baseline density f(.)
#       In real data analysis we use histogram or kernel density based on "survfit"
#  ...: extra arguments of fn
##################################################################
#' Mable fit based on one-sample interval censored data
#' @param data a dataset either \code{data.frame} or an \code{n x 2} matrix.
#' @param M an positive integer or a vector \code{(m0, m1)}. If \code{M = m} or \code{m0 = m1 = m},  
#'   then \code{m} is a preselected degree. If \code{m0 < m1} it specifies the set of 
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.  
#' @param pi0 Initial guess of \eqn{\pi = F(\tau_n)}. Without right censored data, \code{pi0 = 1}. See 'Details'.
#' @param tau right endpoint of support \eqn{[0, \tau)} must be greater than or equal to the maximum observed time
#' @param IC information criterion(s) in addition to Bayesian information criterion (BIC). Current choices are  
#'  "aic" (Akaike information criterion) and/or 
#'  "qhic" (Hannan–Quinn information criterion). 
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate Bernstein/Beta likelihood estimation of density and 
#'  cumulative/survival distributions functions  based on interal censored event time data.
#' @details
#'  Let \eqn{f(t)} and \eqn{F(t) = 1 - S(t)} be the density and cumulative distribution
#'  functions of the event time, respectively. Then \eqn{f(t)} on \eqn{[0, \tau_n]} can be
#'  approximated by \eqn{f_m(t; p) = \tau_n^{-1}\sum_{i=0}^m p_i\beta_{mi}(t/\tau_n)},
#'  where \eqn{p_i \ge 0}, \eqn{i = 0, \ldots, m}, \eqn{\sum_{i=0}^mp_i = 1-p_{m+1}},
#'  \eqn{\beta_{mi}(u)} is the beta denity with shapes \eqn{i+1} and \eqn{m-i+1}, and
#'  \eqn{\tau_n} is the largest observed time, either uncensored time, or right endpoint of 
#'  interval/left censored, or left endpoint of right censored time. We can approximate  
#'  \eqn{S(t)} on \eqn{[0, \tau]} by \eqn{S_m(t; p) = \sum_{i=0}^{m+1} p_i \bar B_{mi}(t/\tau)},  
#'  where  \eqn{\bar B_{mi}(u)}, \eqn{i = 0, \ldots, m}, is the beta survival function with shapes 
#'  \eqn{i+1} and \eqn{m-i+1}, \eqn{\bar B_{m,m+1}(t) = 1}, \eqn{p_{m+1} = 1 - \pi}, and
#'  \eqn{\pi = F(\tau_n)}. For data without right-censored time, \eqn{p_{m+1} = 1-\pi=0}.  
#'  The search for optimal degree \code{m} is stoped if either \code{m1} is reached or the test 
#'  for change-point results in a p-value \code{pval} smaller than \code{sig.level}. 
#' 
#' Each row of \code{data}, \code{(l, u)}, is the interval containing the event time. 
#' Data is uncensored if \code{l = u}, right censored if \code{u = Inf} or \code{u = NA},  
#' and left censored data if \code{l = 0}.
#' @return a class '\code{mable}' object with components
#' \itemize{ 
#'   \item \code{p} the estimated \code{p} with degree \code{m}
#'     selected by the change-point method 
#'   \item \code{mloglik} the maximum log-likelihood at an optimal degree \code{m}
#'   \item \code{interval} support/truncation interval \code{(0, b)}
#'   \item \code{M} the vector \code{(m0,m1)}, where \code{m1} is the last candidate when the search stoped
#'   \item \code{m} the selected optimal degree by the method  of change-point 
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m  \in  \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{tau.n} maximum observed time \eqn{\tau_n}
#'   \item \code{tau} right endpoint of support \eqn{[0, \tau)}
#'   \item \code{ic} a list containing the selected information criterion(s)
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#'   \item \code{convergence} an integer code. 0 indicates successful completion(the iteration is   
#'    convergent). 1 indicates that the maximum candidate degree had been reached in the calculation;
#'   \item \code{delta} the final \code{pval} of the change-point for selecting the optimal degree \code{m};
#' }
#' @seealso \code{\link{mable.group}}
#' @author Zhong Guan <zguan@iusb.edu>
#' @references 
#' Guan, Z. (2019) Maximum Approximate Bernstein Likelihood Estimation in Proportional Hazard Model for Interval-Censored Data, 
#' arXiv:1906.08882 .
#' @examples
#' \donttest{
#'  library(coxinterval) 
#'  bcos=cosmesis
#'  bc.res0<-mable.ic(bcos[bcos$treat=="RT",1:2], M=c(1,50), IC="none")
#'  bc.res1<-mable.ic(bcos[bcos$treat=="RCT",1:2], M=c(1,50), IC="none")
#'  op<-par(mfrow=c(2,2),lwd=2)
#'  plot(bc.res0, which="change-point", lgd.x="right")
#'  plot(bc.res1, which="change-point", lgd.x="right")
#'  plot(bc.res0, which="survival", add=FALSE, xlab="Months", ylim=c(0,1), main="Radiation Only")
#'  legend("topright", bty="n", lty=1:2, col=1:2, c(expression(hat(S)[CP]),
#'                expression(hat(S)[BIC])))
#'  plot(bc.res1, which="survival", add=FALSE, xlab="Months", main="Radiation and Chemotherapy")
#'  legend("topright", bty="n", lty=1:2, col=1:2, c(expression(hat(S)[CP]),
#'                expression(hat(S)[BIC])))
#'  par(op)
#' }
#' @export
mable.ic<-function(data, M, pi0=NULL, tau=Inf, IC=c("none", "aic", "hqic", "all"),
              controls = mable.ctrl(), progress=TRUE){
  IC <- match.arg(IC, several.ok=TRUE)    
  xNames<-deparse(substitute(data))    
  if(ncol(data)>2) stop("data contains too many columns.")
  y<-as.numeric(data[,1]);y2<-as.numeric(data[,2])
  y2[is.na(y2)]<-Inf
  delta<-1*(y<y2)# rvar[,3]
  b<-max(y2[y2<Inf], y);
  if(b>tau) stop("tau must be greater than or equal to the maximum observed time")
  if(is.null(pi0)) pi0<-mean(y2<Inf)
  else if(pi0<=0) stop("pi0 must be in (0, 1]")
  y<-y/b; y2<-y2/b
  if(tau==Inf) y2[y2==Inf]<-.Machine$double.xmax/2
  else y2[y2==Inf]<-tau
  n<-length(y)
  n0<-sum(delta==0)
  n1<-sum(delta==1)
  n<-n0+n1
  N<-c(n0,n1) 
  ord<-order(delta)
  y<-y[ord]; y2<-y2[ord]
  if(missing(M) || length(M)==0) stop("'M' is missing.\n")
  else if(length(M)==1) M<-c(M,M)
  else if(length(M)>=2) {
    M<-c(min(M), max(M))
  }
  k<-M[2]-M[1]
  lk<-rep(0, k+1)
  lr<-rep(0, k)
  bic<-rep(0,k+1)
  pval<-rep(0,k+1)
  level<-controls$sig.level
  chpts<-rep(0,k+1)
  eps<-c(controls$eps, .Machine$double.eps)
  convergent<-0
  #if(!any(y2>1)) pi0<-1
  p<-c(rep(pi0, M[2]+1)/(M[2]+1), 1-pi0)
  optim<-0
  del<-0
  ## Call C mable_ic
  res<-.C("mable_ic", as.integer(M), as.double(pi0), as.double(y), as.double(y2), 
    as.integer(N), as.double(lk), as.double(lr), as.double(p), as.double(eps), 
    as.integer(controls$maxit), as.logical(progress), as.double(pval),  
    as.double(bic), as.integer(chpts), as.integer(optim), as.double(level), 
    as.integer(convergent), as.double(del))
  M<-res[[1]]
  k<-M[2]-M[1]
  lk<-res[[6]][1:(k+1)]-n0*log(b)
  m<-res[[15]]; 
  mllik<-lk[m-M[1]+1]
  ans<-list(m=m, mloglik=mllik-n0*log(b), tau.n=b, tau=tau, interval=c(0, b), 
      convergence=res[[17]], delta=res[[12]][k+1])
  if(k==0) ans$p<-res[[8]][1:(m[1]+2)]
  if(k>0){
    bic<-res[[13]][1:(k+1)]
    ic<-list()
    ic$BIC<-bic
    if(!any(IC=="none")) {
      d<-2*(lk-bic)/log(n)
      if(any(IC=="aic")|| any(IC=="all")){
        ic$AIC<-((log(n)-2)*lk+2*bic)/log(n)#lk-2*(lk-bic)/log(n)
        #aic<-lk-d-(d^2+d)/(n-d-1)
      }
      if(any(IC=="qhic")|| any(IC=="all")){
        ic$QHC<-lk-2*(lk-bic)*log(log(n))/log(n)}
    }
    ans$lr<-res[[7]][1:k] 
    ans$p<-res[[8]][1:(m+2)]         
    #ans$p<-list(p.cp=res[[8]][1:(m[1]+2)], p.bic=res[[8]][(m[1]+3):(m[1]+m[2]+4)])          
    ans$M<-M; ans$lk<-lk; ans$pval<-res[[12]][1:(k+1)]; ans$ic<-ic;
    ans$chpts<-res[[14]][1:(k+1)]+M[1] 
  }  
  ans$xNames<-xNames 
  ans$data.type<-"icen"
  class(ans)<-"mable"
  return(ans)
}

##############################################
#' Plot mathod for class 'mable_reg'
##############################################
#' @param x a class 'mable_reg' object return by functions such as \code{mable.ph} which contains 
#'  \code{M}, \code{coefficients}, \code{p}, \code{m}, \code{x0}, \code{tau.n}, \code{tau} 
#'  \code{lk}, \code{lr}.
#' @param y a new data.frame of covariate value(s) as row(s), whose columns are
#'          arranged in the same order as in the \code{formula} called by the function
#'          that returned the object \code{x}.
#' @param newdata a new data.frame (ignored if \code{y} is included), imputed
#'          by the working baseline \code{x0} if both missing. 
#' @param ntime number of evaluations of density, survival or cumulative distribution
#'         curve to be plotted.
#' @param xlab x-axis label 
#' @param which indicates which graphs to plot, options are 
#'  "survival", "likelihood", "change-point", "density", or "all". If not "all", 
#'  \code{which} can contain more than one options.
#' @param add logical add to an existing plot or not
#' @param ... additional arguments to be passed to the base plot function
#' @author Zhong Guan <zguan@iusb.edu>
#' @importFrom stats dexp pexp 
#' @export  
plot.mable_reg<-function(x, y, newdata =NULL, ntime=512, xlab="Time",
      which=c("survival", "likelihood", "change-point", "density", "all"), add=FALSE,...){
  which <- match.arg(which, several.ok=TRUE)
  #model<-substr(x$model, 7,9) 
  if(any(which=="likelihood")||any(which=="all")) {
    if(is.null(x$lk)) stop("Cannot plot 'likelihood'.")
    if(!add) plot(x$M[1]:x$M[2], x$lk, type="p", xlab="m",ylab="Loglikelihood",
        main="Loglikelihood")
    else points(x$M[1]:x$M[2], x$lk, pch=1,...)
    segments(x$m[1], min(x$lk), x$m[1], x$mloglik[1], lty=2)
    axis(1, x$m[1], as.character(x$m[1]), col=2)
  }
  if(any(which=="change-point")||any(which=="all")) {
    if(is.null(x$lr)) stop("Cannot plot 'likelihood ratios'.")
    if(!add) plot((x$M[1]+1):x$M[2], x$lr, type="p", xlab="m",ylab="Loglikelihood Ratio",
        main="Change-Point")
    else points((x$M[1]+1):x$M[2], x$lr, pch=1, ...)
    segments(x$m[1], 0, x$m[1], max(x$lr), lty=2)
    axis(1, x$m[1], as.character(x$m[1]), col=2)
  }
  if(any(which=="density")||any(which=="survival")||any(which=="all")) {
    if(missing(y) && is.null(newdata)) {
      y<-data.frame(t(x$x0))
      message("missing y and newdata, assigned as x=x0")}
    else if(missing(y)) y<-newdata
    nr=dim(y)[1]
    tau.n<-x$tau.n
    tau<-x$tau
    #cat("y=",y[1,],"\n")
    gama<-x$coefficients
    p<-x$p
    m<-x$m[1]
    #cat("p=",p,"\n")
    #cat("m=",m,"\n")
    if(tau==Inf) tau<-tau.n
    time<-seq(0, tau, length=ntime)
    for(i in 1:nr){
      if(x$model == "ph"){	          
        xlb<-time[time<=tau.n]
        xgb<-time[time>tau.n]
        rate<-(m+1)*p[m+1]/p[m+2]/tau.n
        Sb<-c(1-suppressMessages(pmixbeta(xlb, p=p[-(m+2)], c(0, tau.n))), p[m+2]*(1-pexp(xgb-tau.n,rate)))
        fb<-c(suppressMessages(dmixbeta(xlb, p=p[-(m+2)], c(0, tau.n))), p[m+2]*dexp(xgb-tau.n, rate))
        egxt<-as.vector(exp(sum(y[i,]*gama))/x$egx0)
        fbx<-egxt*fb*Sb^(egxt-1)              
        Sbx<-Sb^egxt}
      if(x$model == "aft"){	
        Sbx<-1-pmixbeta(time, p=p, c(0, tau))
        fbx<-dmixbeta(time, p=p, c(0, tau))
        egxt<-as.vector(exp(sum(y[i,]*gama))/x$egx0)
        time<-time/egxt
        tau<-tau.n}
      if(x$model == "po"){	          
        xlb<-time[time<=tau.n]
        xgb<-time[time>tau.n]
        rate<-(m+1)*p[m+1]/p[m+2]/tau.n
        Sb<-c(1-suppressMessages(pmixbeta(xlb, p=p[-(m+2)], c(0, tau.n))), p[m+2]*(1-pexp(xgb-tau.n,rate)))
        fb<-c(suppressMessages(dmixbeta(xlb, p=p[-(m+2)], c(0, tau.n))), p[m+2]*dexp(xgb-tau.n, rate))
        egxt<-as.vector(exp(sum(y[i,]*gama))/x$egx0)
        eta<-x$eta
        fbx<-egxt*fb/(egxt+(1-egxt)*Sb^eta)^(1+1/eta)              
        Sbx<-Sb/(egxt+(1-egxt)*Sb^eta)^(1/eta)}
      if(any(which=="density")||any(which=="all"))
        if(!add && i==1) plot(time, fbx, xlab=xlab, ylab="Density", xlim=c(0,tau),...)
        else lines(time, fbx, xlab=xlab, ylab="Density",  xlim=c(0,tau), ...)
      if(any(which=="survival")||any(which=="all")) # default is "survival" for plotting "all" 
        if(!add && i==1) plot(time, Sbx, xlab=xlab, ylab="Probability", xlim=c(0,tau),  ylim=c(0,1), ...)
        else lines(time, Sbx, xlab=xlab, ylab="Probability",xlim=c(0,tau),  ylim=c(0,1),...)  
    }
  }
}
######################
#' readingCall: an internal function of 'icenReg' package  
#' @keywords internal
#' @noRd
readingCall <- function(mf){
  m <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf$formula <- quote(formula)
  mf$data <- quote(data)
  mf$na.action <- quote(na.pass)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  ans <- list(mt = mt, mf = mf)
  return(ans)
}
######################
#' Getting data from formula and data  
#' @keywords internal
#' @importFrom stats model.matrix model.response 
#' @noRd
get.mableData<-function(formula, data){
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    callInfo <- readingCall(mf)
    mf <- callInfo$mf
    mt <- callInfo$mt
    rvar <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)
    if (is.matrix(x)) xNames <- colnames(x)
    else {
        xNames <- as.character(formula[[3]])
        x<-matrix(x, ncol=1)}
    if ("(Intercept)" %in% colnames(x)) {
        ind <- which(colnames(x) == "(Intercept)")
        x <- x[, -ind]
        xNames <- xNames[-ind]
    } # This should not happen.
    y<-rvar[,1]
    y2<-rvar[,2]
    y2[is.na(y2)]<-Inf
    delta<-1*(y<y2)# rvar[,3]
    if (sum(is.na(x)) > 0) 
        stop("NA's not allowed in covariates")
    callText <- mf$formula
    out<-list(x=x, y=y, y2=y2, delta=delta, callText=callText, xNames=xNames)
    return(out)
}
######################################################
# Wrap all MABLE for regression fit in one function
######################################################
#' Mable fit of semiparametric regression model based on interval censored data 
#' @description Wrapping all code{mable} fit of regression models in one function.
#' Using maximum approximate Bernstein/Beta likelihood
#' estimation to fit semiparametric regression models: Cox ph model,
#' proportional odds(po) model, accelerated failure time model, and so on.
#' @param formula regression formula. Response must be of the form \code{cbind(l, u)}.  See 'Details'.
#' @param data a dataset
#' @param model the model to fit. Current options are "\code{ph}"
#'  (Cox PH) or "\code{aft}" (accelerated failure time model)
#' @param M a vector \code{(m0, m1)} specifies the set of consective integers as candidate degrees
#' @param g  an initial guess of the regression coefficients 
#' @param pi0 Initial guess of \eqn{\pi(x_0) = F(\tau_n|x_0)}. Without right censored data, \code{pi0 = 1}. See 'Details'.
#' @param tau right endpoint of support \eqn{[0, \tau)} must be greater than or equal to the maximum observed time
#' @param x0 a working baseline covariate. See 'Details'. 
#' @param eta the given positive value of \eqn{\eta}. Used when \code{model="po"}.
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit 
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @author Zhong Guan <zguan@iusb.edu>
#' @return A 'mable_reg' class object
#' @details For "\code{ph}"  model a missing initial guess of the regression coefficients 
#'    \code{g} is obtained by \code{ic_sp()} of package \code{icenReg}. For "\code{aft}" model a
#'    missing \code{g} is imputed by the rank estimate \code{aftsrr()} of package \code{aftgee} 
#'    for right-censored data. For general interval censored observations, we keep the 
#'    right-censored but replace the finite interval with its midpoint and fit the data by 
#'    \code{aftsrr()} as a right-censored data.  
#' 
#' @keywords distribution models nonparametric regression smooth survival
#' @concept Cox proportional hazards model 
#' @concept accelerated failure time model 
#' @concept interval censoring
#' @seealso \code{\link{mable.aft}}, \code{\link{mable.ph}} 
#' @export
mable.reg<-function(formula, data, model=c("ph","aft"), M, g=NULL, pi0=NULL, tau=Inf, 
        x0=NULL, eta=1, controls = mable.ctrl(), progress = TRUE){
    model=match.arg(model)
    out<-switch(model,
        ph = mable.ph(formula, data, M, g, pi0, tau, x0, controls, progress), 
        #po = mable.po(formula, data, M, g, pi0, tau, x0, eta, controls, progress), 
        aft = mable.aft(formula, data, M, g, tau, x0, controls, progress))
    class(out)<-"mable_reg" # 
    return(out)
} 
