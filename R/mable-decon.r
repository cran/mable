#############################################################
##                 MABLE Deconvolution                     ##
##   Maximum Approximate Bernstein Likelihood Estimation   ##
##              for Density Deconvolution                  ##
#############################################################
## Y = X + e, X~F unknown, e~G known/unknown, X and e are independent
## Yi = Xi + ei, i=1,...,n. Y1,...,Yn are observations.
## We want to estimate f=F' based on Y1,...,Yn.
#  setwd("C:\\Users\\zguan\\Documents\\papers\\bernstein polynomials\\deconvolution\\density deconvolution\\C")
#  setwd("C:\\Documents and Settings\\zguan\\My Documents\\papers\\deconvolution\\density deconvolution\\C")
#  dyn.load("em-algorithm")
## Call C functions by .External
# R CMD SHLIB mable-decon.c
# R --arch x64 CMD SHLIB mable-decon.c
######################################################################################################
#' Mable deconvolution with a known error density
#' @param y vector of observed data values
#' @param gn error density function if known, default is NULL if unknown
#' @param ... additional arguments to be passed to gn
#' @param M a vector \code{(m0, m1)} specifies the set of consective candidate model degrees, \code{M = m0:m1}.
#'    If \code{gn} is unknown then \code{M} a 2 x 2 matrix whose rows \code{(m0,m1)} and \code{(k0,k1)}
#'   specify lower and upper bounds for degrees \code{m} and \code{k}, respectively.
#' @param interval a finite vector \code{(a,b)}, the endpoints of supporting/truncation interval
#'   if \code{gn} is known. Otherwise, it is a 2 x 2 matrix whose rows \code{(a,b)} and \code{(a1,b1)}
#'   specify supporting/truncation intervals of \code{X} and \eqn{\epsilon}, respectively. See Details.
#' @param IC information criterion(s) in addition to Bayesian information criterion (BIC). Current choices are
#'  "aic" (Akaike information criterion) and/or
#'  "qhic" (Hannan–Quinn information criterion).
#' @param vanished logical whether the unknown error density vanishes at both end-points of \code{[a1,b1]}
# @param new logical whether to use new c procedure or not
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate Bernstein/Beta likelihood estimation in additive density deconvolution model
#' with a known error density.
#' @details
#' Consider the additive measurement error model \eqn{Y = X + \epsilon}, where
#' \eqn{X} has an unknown distribution \eqn{F} on a known support \code{[a,b]}, \eqn{\epsilon} has a known or unknown distribution \eqn{G},
#' and \eqn{X} and \eqn{\epsilon} are independent. We want to estimate density \eqn{f = F'}
#' based on independent observations, \eqn{y_i = x_i + \epsilon_i}, \eqn{i = 1, \ldots, n}, of \eqn{Y}.
#' We approximate \eqn{f} by a Bernstein polynomial model on \code{[a,b]}. If \eqn{g=G'} is unknown on
#' a known support \code{[a1,b1]}, then we approximate \eqn{g} by a Bernstein polynomial model on
#' \code{[a1,b1]}, \eqn{a1<0<b1}. We assume \eqn{E(\epsilon)=0}. AIC and BIC methods are used to
#' select model degrees \code{(m,k)}.
#' @return A \code{mable} class object with components, if \eqn{g} is known,
#' \itemize{
#'      \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last candidate degree when the search stoped
#'      \item \code{m} the selected optimal degree \code{m}
#'      \item \code{p} the estimate of \code{p = (p_0, \dots, p_m)}, the coefficients of Bernstein polynomial of degree \code{m}
#'      \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'      \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'      \item \code{convergence} An integer code. 0 indicates an optimal degree
#'        is successfully selected in \code{M}. 1 indicates that the search stoped at \code{m1}.
#'       \item \code{ic} a list containing the selected information criterion(s)
#'      \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'      \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' if \eqn{g} is unknown,
#' \itemize{
#'      \item \code{M} the 2 x 2 matrix with rows \code{(m0, m1)} and \code{(k0,k1)}
#'      \item \code{nu_aic} the selected optimal degrees \code{(m,k)} using AIC method
#'      \item \code{p_aic} the estimate of \code{p = (p_0, \dots, p_m)}, the coefficients
#'       of Bernstein polynomial model for \eqn{f} of degree \code{m} as in \code{nu_aic}
#'      \item \code{q_aic} the estimate of \code{q = (q_0, \dots, q_k)}, the coefficients
#'       of Bernstein polynomial model for \eqn{g} of degree \code{k} as in \code{nu_aic}
#'      \item \code{nu_bic} the selected optimal degrees \code{(m,k)} using BIC method
#'      \item \code{p_bic} the estimate of \code{p = (p_0, \dots, p_m)}, the coefficients
#'       of Bernstein polynomial model for \eqn{f} of degree \code{m} as in \code{nu_bic}
#'      \item \code{q_bic} the estimate of \code{q = (q_0, \dots, q_k)}, the coefficients
#'       of Bernstein polynomial model for \eqn{g} of degree \code{k} as in \code{nu_bic}
#'      \item \code{lk} matrix of log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'           and \eqn{k \in \{k_0, \ldots, k_1\}}
#'      \item \code{aic} a matrix containing the Akaike information criterion(s) at
#'        \eqn{m \in \{m_0, \ldots, m_1\}} and \eqn{k \in \{k_0, \ldots, k_1\}}
#'      \item \code{bic} a matrix containing the Bayesian information criterion(s) at
#'        \eqn{m \in \{m_0, \ldots, m_1\}} and \eqn{k \in \{k_0, \ldots, k_1\}}
#' }
#' @author Zhong Guan <zguan@iu.edu>
#' @references
#' Guan, Z., (2019) Fast Nonparametric Maximum Likelihood Density Deconvolution Using Bernstein Polynomials, Statistica Sinica,
#' doi:10.5705/ss.202018.0173
#' @examples
#' \donttest{
#'  # A simulated normal dataset
#'  set.seed(123)
#'  mu<-1; sig<-2; a<-mu-sig*5; b<-mu+sig*5;
#'  gn<-function(x) dnorm(x, 0, 1)
#'  n<-50;
#'  x<-rnorm(n, mu, sig); e<-rnorm(n); y<-x+e;
#'  res<-mable.decon(y, gn, interval = c(a, b), M = c(5, 50))
#'  op<-par(mfrow = c(2, 2),lwd = 2)
#'  plot(res, which="likelihood")
#'  plot(res, which="change-point", lgd.x="topright")
#'  plot(xx<-seq(a, b, length=100), yy<-dnorm(xx, mu, sig), type="l", xlab="x",
#'      ylab="Density", ylim=c(0, max(yy)*1.1))
#'  plot(res, which="density", types=c(2,3), colors=c(2,3))
#'  # kernel density based on pure data
#'  lines(density(x), lty=4, col=4)
#'  legend("topright", bty="n", lty=1:4, col=1:4,
#'  c(expression(f), expression(hat(f)[cp]), expression(hat(f)[bic]), expression(tilde(f)[K])))
#'  plot(xx, yy<-pnorm(xx, mu, sig), type="l", xlab="x", ylab="Distribution Function")
#'  plot(res, which="cumulative",  types=c(2,3), colors=c(2,3))
#'  legend("bottomright", bty="n", lty=1:3, col=1:3,
#'      c(expression(F), expression(hat(F)[cp]), expression(hat(F)[bic])))
#'  par(op)
#' }
#' @concept Additive measurement error
#' @concept Bernstein polynomial model
#' @concept Density deconvolution
#' @export
mable.decon<-function(y, gn=NULL, ...,  M, interval=c(0, 1),
    IC=c("none", "aic", "hqic", "all"), vanished = TRUE,
    controls=mable.ctrl(maxit.em=1e5, eps.em=1e-5, maxit.nt=1e2, eps.nt=1e-10),
    progress=TRUE){
    cat("Starting Deconvolution ... \n")
    level<-controls$sig.level
    yName<-deparse(substitute(y))
    if(!is.null(gn)){
        if(M[2]-M[1]<=2) stop("Too few candidate model degrees.")
        IC <- match.arg(IC, several.ok=TRUE)
        a<-interval[1]
        b<-interval[2]
        if(a>=b) stop("'a' must be smaller than 'b'")
        if(any(y<a) | any(y>b)) stop("'interval' must contains all 'y'.")
        y<-(y-a)/(b-a)
        gn <- match.fun(gn)
        ff <- function(x) (b-a)*gn((b-a)*x,...)
        ## Call mable_decon
        wk<-.External("mable_decon", ff, rho = environment(), as.double(y),
            as.integer(M), as.double(controls$eps), as.integer(controls$maxit),
            as.logical(progress), as.double(level), as.double(.Machine$double.eps))
        res <- wk[c("lk", "lr", "p", "m", "pval", "bic", "chpts", "M")]
        res$interval<-interval
        m<-res$m
        lk<-res$lk
        bic<-res$bic
        res$mloglik<-lk[res$m-M[1]+1]
        #res$p<-list(p.cp=res$p[1:(m[1]+1)], p.bic=res$p[(m[1]+2):(m[1]+m[2]+2)])
        res$convergence<-1*(M[2]==res$M[2])
        ic<-list()
        ic$BIC<-bic
        n<-length(y)
        if(!any(IC=="none")){
            d<-2*(lk-bic)/log(n)
            if(any(IC=="aic")|| any(IC=="all")){
                ic$AIC<-((log(n)-2)*lk+2*bic)/log(n)#lk-2*(lk-bic)/log(n)
                #aic<-lk-d-(d^2+d)/(n-d-1)
            }
            if(any(IC=="qhic")|| any(IC=="all")){
                ic$QHC<-lk-2*(lk-bic)*log(log(n))/log(n)}
        }
        res$ic<-ic
        res$bic<-NULL
    }
    else{
        if(!is.matrix(interval)||!all(dim(interval)==c(2,2)))
            stop("'interval' must be a 2x2 matrix.")
        a<-interval[1,1]
        b<-interval[1,2]
        a1<-interval[2,1]
        b1<-interval[2,2]
        if(a>=b) stop("'a' must be smaller than 'b'")
        if(a1>=b1 || a1>=b1 || b1<=0) stop("'a1' and 'b1' must satisfy a1<0<b1")
        if(any(y<a+a1) || any(y>b+b1)) stop("All 'y' must be in [a+a1, b+b1].")
        res0<-mable(y, M=c(2,max(100,m)), c(a+a1, b+b1), progress =FALSE)
        y<-(y-a)/(b-a)
        interval<-c(a1,b1)/(b-a)
        if(!is.matrix(M)||!all(dim(M)==c(2,2)))
            stop("'M' must be a 2x2 matrix.")
        m0<-M[1,1]
        m<-M[1,2]
        k0<-M[2,1]
        if(vanished && k0<=2) {
            cat("Note: For vanished error density, k0 should be grater than 2. Change k0 to 3.")
            k0<-3}
        k<-M[2,2]
        ## Call optim_decon
        wk<-.External("optim_decon", as.double(y), as.double(interval), as.logical(vanished),
            as.integer(M), as.double(res0$p), as.integer(res0$m), as.double(controls$eps.em),
            as.integer(controls$maxit.em), as.logical(progress),
            as.double(controls$eps.nt), as.integer(controls$maxit.nt))
        res <- wk[c("lk", "np", "D", "aic", "bic", "nu_d", "nu_aic", "nu_bic",
                    "p_d", "q_d", "p_aic", "q_aic", "p_bic", "q_bic")]
        nr<-m-m0+1
        nc<-k-k0+1
        res$pi<-res0$p
        res$v<-res0$m
        res$interval<-res0$interval
        res$lk<-matrix(res$lk, nrow=nr, ncol=nc)
        res$np<-matrix(res$np, nrow=nr, ncol=nc)
        res$D<-matrix(res$D, nrow=nr, ncol=nc)
        res$aic<-matrix(res$aic, nrow=nr, ncol=nc)
        res$bic<-matrix(res$bic, nrow=nr, ncol=nc)
        res$lk_aic<-res$lk[res$nu_aic[1]-m0+1, res$nu_aic[2]-k0+1]
        res$lk_bic<-res$lk[res$nu_bic[1]-m0+1, res$nu_bic[2]-k0+1]
    }
    res$xNames<-yName
    res$data.type<-"noisy"
    class(res) <- "mable"
    res
}
 
