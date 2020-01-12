#############################################################
##                 MABLE Deconvolution                     ##
##   Maximum Approximate Bernstein Likelihood Estimation   ##
##              for Density Deconvolution                  ##
#############################################################
## Y = X + e, X~F unknown, e~G known, X and e are independent
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
#' @param y observaed data
#' @param gn error density function
#' @param ... additional arguments to be passed to gn
#' @param M a vector \code{(m0, m1)} specifies the set of consective candidate model degrees, \code{M=m0:m1}.
#' @param interval a finite vector containing the endpoints of supporting/truncation interval
#' @param IC information criterion(s) in addition to Bayesian information criterion (BIC). Current choices are
#'  "aic" (Akaike information criterion) and/or
#'  "qhic" (Hannan–Quinn information criterion).
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and other control options. Default is \code{\link{mable.ctrl}}.
#' @param progress if \code{TRUE} a text progressbar is displayed
#' @description Maximum approximate Bernstein/Beta likelihood estimation in additive density deconvolution model
#' with a known error density.
#' @details
#' Consider the additive measurement error model \eqn{Y = X + \epsilon}, where
#' \eqn{X} has an unknown distribution \eqn{F}, \eqn{\epsilon} has a known distribution \eqn{G},
#' and \eqn{X} and \eqn{\epsilon} are independent. We want to estimate density \eqn{f=F'}
#' based on independent observations, \eqn{y_i = x_i + \epsilon_i}, \eqn{i=1,\ldots,n}, of \eqn{Y}.
#' @return A \code{mable} class object with components
#' \itemize{
#'      \item \code{M} the vector \code{(m0,m1)}, where \code{m1} is the last candidate degree when the search stoped
#'      \item \code{m} the selected optimal degree \code{m}
#'      \item \code{p} the estimate of \code{p=(p_0,\dots,p_m)}, the coefficients of Bernstein polynomial of degree \code{m}
#'      \item \code{lk} log-likelihoods evaluated at \eqn{m\in\{m_0,\ldots, m_1\}}
#'      \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m\in\{m_0+1,\ldots, m_1\}}
#'      \item \code{convergence} An integer code. 0 indicates an optimal degree
#'        is successfully selected in \code{M}). 1 indicates that the search stoped at \code{m1}.
#'       \item \code{ic} a list containing the selected information criterion(s)
#'      \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'      \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iusb.edu>
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
#'  res<-mable.decon(y, gn, interval=c(a,b), M=c(5, 50))
#'  op<-par(mfrow=c(2,2),lwd=2)
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
mable.decon<-function(y, gn, ...,  M, interval=c(0, 1),
    IC=c("none", "aic", "hqic", "all"),
    controls=mable.ctrl(maxit=50000, eps=1e-7), progress=TRUE){
    #dyn.load("mable-decon")
    if(M[2]-M[1]<=2) stop("Too few candidate model degrees.")
    IC <- match.arg(IC, several.ok=TRUE)
    gn <- match.fun(gn)
    yName<-deparse(substitute(y))
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("'a' must be smaller than 'b'")
    if(any(y<a) | any(y>b)) stop("'interval' must contains all 'y'.")
    y<-(y-a)/(b-a)
#    gn<-function(x) (b-a)*gn(a+(b-a)*x,...)
    ff <- function(x) (b-a)*gn((b-a)*x,...)
    level=controls$sig.level
    ## Call mable_decon
    wk<-.External("mable_decon", ff, rho = environment(), as.double(y),
        as.integer(M), as.double(controls$eps), as.integer(controls$maxit),
        as.logical(progress), as.double(level), as.double(.Machine$double.eps))
    res <- wk[c("lk", "lr", "p", "m", "pval", "bic", "chpts", "M")]
    res$support<-interval
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
    res$xNames<-yName
    res$ic<-ic
    res$bic<-NULL
    res$data.type<-"noisy"
    class(res) <- "mable"
    res
}
