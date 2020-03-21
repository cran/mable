#setwd("C:\\Documents and Settings\\zguan\\My Documents\\papers\\bernstein polynomials\\C")
#setwd("C:\\Users\\zguan\\Documents\\papers\\bernstein polynomials\\C")
#dyn.load("mable")
#dyn.unload("mable")
## Call C functions by .C
# R CMD SHLIB mable.c
# R --arch x64 CMD SHLIB mable.c
################################################
#          One-sample raw data
################################################
#' Mable fit of one-sample raw data with an optimal or given degree.
#' @param x a (non-empty) numeric vector of data values.
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m} or \code{m0 = m1 = m},
#'   then \code{m} is a preselected degree. If \code{m0<m1} it specifies the set of
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.
#' @param interval a vector containing the endpoints of supporting/truncation interval
#' @param IC information criterion(s) in addition to Bayesian information criterion (BIC). Current choices are
#'  "aic" (Akaike information criterion) and/or
#'  "qhic" (Hannan–Quinn information criterion).
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.
#' @param progress if TRUE a text progressbar is displayed
#' @description  Maximum approximate Bernstein/Beta likelihood estimation based on
#'  one-sample raw data with an optimal selected by the change-point method among \code{m0:m1}
#'  or a preselected model degree \code{m}.
#' @details
#'  Any continuous density function \eqn{f} on a known closed supporting interval \eqn{[a,b]} can be
#'  estimated by Bernstein polynomial \eqn{f_m(x; p) = \sum_{i=0}^m p_i\beta_{mi}[(x-a)/(b-a)]/(b-a)},
#'  where \eqn{p = (p_0, \ldots, p_m)}, \eqn{p_i \ge 0}, \eqn{\sum_{i=0}^m p_i = 1} and
#'  \eqn{\beta_{mi}(u) = (m+1){m\choose i}u^i(1-x)^{m-i}}, \eqn{i = 0, 1, \ldots, m},
#'  is the beta density with shapes \eqn{(i+1, m-i+1)}.
#'  For each \code{m}, the MABLE of the coefficients \code{p}, the mixture proportions, are
#'  obtained using EM algorithm. The EM iteration for each candidate \code{m} stops if either
#'  the total absolute change of the log likelihood and the coefficients of Bernstein polynomial
#'  is smaller than \code{eps} or the maximum number of iterations \code{maxit} is reached.
#'
#'  If \code{m0<m1}, an optimal model degree is selected as the change-point of the increments of
#'  log-likelihood, log likelihood ratios, for \eqn{m \in \{m_0, m_0+1, \ldots, m_1\}}. Alternatively,
#'  one can choose an optimal degree based on the BIC (Schwarz, 1978) which are evaluated at
#'  \eqn{m \in \{m_0, m_0+1, \ldots, m_1\}}. The search for optimal degree \code{m} is stoped if either
#'  \code{m1} is reached with a warning or the test for change-point results in a p-value \code{pval}
#'  smaller than \code{sig.level}.  The BIC for a given degree \code{m} is calculated as in
#'  Schwarz (1978) where the dimension of the model is \eqn{d = \#\{i: \hat p_i\ge\epsilon,
#'   i = 0, \ldots, m\} - 1} and a default \eqn{\epsilon} is chosen as \code{.Machine$double.eps}.
#' @return  A list with components
#' \itemize{
#'   \item \code{m} the selected/given optimal degree by methods of change-point
#'   \item \code{p} the estimated vector of mixture proportions \eqn{p = (p_0, \ldots, p_m)}
#'       with the selected/given optimal degree \code{m}
#'   \item \code{mloglik}  the maximum log-likelihood at degree \code{m}
#'   \item \code{interval} support/truncation interval \code{(a,b)}
#'   \item \code{convergence} An integer code. 0 indicates successful completion (all the EM iterations are convergent and an optimal degree
#'     is successfully selected in \code{M}). Possible error codes are
#'    \itemize{
#'       \item 1, indicates that the iteration limit \code{maxit} had been reached in at least one EM iteration;
#'       \item 2, the search did not finish before \code{m1}.
#'     }
#'   \item \code{delta} the convergence criterion \code{delta} value
#'  }
#'  and, if \code{m0<m1},
#' \itemize{
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1}, if greater than \code{m0}, is the
#'      largest candidate when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{ic} a list containing the selected information criterion(s)
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @note Since the Bernstein polynomial model of degree \eqn{m} is nested in the model of
#' degree \eqn{m+1}, the maximum likelihood is increasing in \eqn{m}. The change-point method
#' is used to choose an optimal degree \eqn{m}.
#' @author Zhong Guan <zguan@iusb.edu>
#' @references Guan, Z. (2016) Efficient and robust density estimation using Bernstein type polynomials. \emph{Journal of Nonparametric Statistics}, 28(2):250-271.
#' @examples
#' \donttest{
#' # Vaal Rive Flow Data
# load("Vaal.Flow.rdata")
#'  data(Vaal.Flow)
#'  x<-Vaal.Flow$Flow
#'  res<-mable(x, M = c(2,100), interval = c(0, 3000), controls =
#'         mable.ctrl(sig.level = 1e-8, maxit = 2000, eps = 1.0e-9))
#'  op<-par(mfrow = c(1,2),lwd = 2)
#'  layout(rbind(c(1, 2), c(3, 3)))
#'  plot(res, which = "likelihood", cex = .5)
#'  plot(res, which = c("change-point"), lgd.x = "topright")
#'  hist(x, prob = TRUE, xlim = c(0,3000), ylim = c(0,.0022), breaks = 100*(0:30),
#'   main = "Histogram and Densities of the Annual Flow of Vaal River",
#'   border = "dark grey",lwd = 1,xlab = "x", ylab = "f(x)", col  = "light grey")
#'  lines(density(x, bw = "nrd0", adjust = 1), lty = 4, col = 4)
#'  lines(y<-seq(0, 3000, length = 100), dlnorm(y, mean(log(x)),
#'                    sqrt(var(log(x)))), lty = 2, col = 2)
#'  plot(res, which = "density", add = TRUE)
#'  legend("top", lty = c(1, 2, 4), col = c(1, 2, 4), bty = "n",
#'  c(expression(paste("MABLE: ",hat(f)[B])),
#'         expression(paste("Log-Normal: ",hat(f)[P])),
#'                expression(paste("KDE: ",hat(f)[K]))))
#'  par(op)
#' }
#' \donttest{
#' # Old Faithful Data
#'  library(mixtools)
#'  x<-faithful$eruptions
#'  a<-0; b<-7
#'  v<-seq(a, b,len = 512)
#'  mu<-c(2,4.5); sig<-c(1,1)
#'  pmix<-normalmixEM(x,.5, mu, sig)
#'  lam<-pmix$lambda; mu<-pmix$mu; sig<-pmix$sigma
#'  y1<-lam[1]*dnorm(v,mu[1], sig[1])+lam[2]*dnorm(v, mu[2], sig[2])
#'  res<-mable(x, M = c(2,300), interval = c(a,b), controls  =
#'         mable.ctrl(sig.level = 1e-8, maxit = 2000, eps = 1.0e-7))
#'  op<-par(mfrow = c(1,2),lwd = 2)
#'  layout(rbind(c(1, 2), c(3, 3)))
#'  plot(res, which = "likelihood")
#'  plot(res, which = "change-point")
#'  hist(x, breaks = seq(0,7.5,len = 20), xlim = c(0,7), ylim = c(0,.7),
#'      prob  = TRUE,xlab = "t", ylab = "f(t)", col  = "light grey",
#'      main = "Histogram and Density of
#'                Duration of Eruptions of Old Faithful")
#'  lines(density(x, bw = "nrd0", adjust = 1), lty = 4, col = 4, lwd = 2)
#'  plot(res, which = "density", add = TRUE)
#'  lines(v, y1, lty = 2, col = 2, lwd = 2)
#'  legend("topright", lty = c(1,2,4), col = c(1,2,4), lwd = 2, bty = "n",
#'       c(expression(paste("MABLE: ",hat(f)[B](x))),
#'          expression(paste("Mixture: ",hat(f)[P](t))),
#'          expression(paste("KDE: ",hat(f)[K](t)))))
#'  par(op)
#' }
#'
#' @keywords distribution  models  nonparametric  smooth univar
#' @concept Bernstein polynomial model
# @useDynLib mable, .registration = TRUE
#' @export
mable<-function(x, M, interval = c(0,1), IC = c("none", "aic", "hqic", "all"),
        controls = mable.ctrl(), progress = TRUE){
    IC <- match.arg(IC, several.ok = TRUE)
    n<-length(x)
    xName<-deparse(substitute(x))
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    x<-(x-a)/(b-a)
    convergent<-0
    del<-controls$eps
    eps<-c(controls$eps, .Machine$double.eps)
    llik<-0;
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
    k<-M[2]-M[1]
    if(k>0 && k<=3) stop("Too few candidate model degrees.")
    if(k==0){
        m<-M[1]
        p<-rep(1,m+1)/(m+1)
        ## Call C mable_em
        res<-.C("mable_em",
          as.integer(m), as.integer(n), as.double(p), as.double(x), as.integer(controls$maxit),
          as.double(controls$eps),  as.double(llik), as.logical(convergent), as.double(del))
        ans<-list(p=res[[3]], m=m, mloglik=res[[7]], interval=c(a,b), convergent=!res[[8]], del=res[[9]])
    }
    else{
        #p<-rep(1, 2*M[2]+2)
        p<-rep(1, M[2]+1)
        lk<-rep(0, k+1)
        lr<-rep(0,k)
        bic<-rep(0,k+1)
        pval<-rep(0,k+1)
        level<-controls$sig.level
        chpts<-rep(0,k+1)
        optim<-0#c(0,0)
        ## Call C mable_optim
        res<-.C("mable_optim",
          as.integer(M), as.integer(n), as.double(p), as.double(x), as.integer(controls$maxit),
          as.double(eps), as.double(lk), as.double(lr), as.integer(optim),
          as.double(pval), as.double(bic), as.integer(chpts), as.double(controls$tini),
          as.logical(progress), as.integer(convergent), as.double(del), as.double(level))
        M<-res[[1]]; k<-M[2]-M[1]
        m<-res[[9]]
        mloglik<-res[[7]][m-M[1]+1]
        lk<-res[[7]][1:(k+1)]
        bic<-res[[11]][1:(k+1)]
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
        ans<-list(p=res[[3]][1:(m[1]+1)],
            m=m, mloglik=mloglik, lk=lk, lr=res[[8]][1:k], M=M,
            interval=c(a,b), pval=res[[10]][1:(k+1)], ic=ic,
            chpts=res[[12]][1:(k+1)]+M[1], convergence=res[[15]], delta=res[[16]])
    }
    ans$xNames<-xName
    ans$data.type<-"raw"
    class(ans)<-"mable"
    return(ans)
}
################################################
#        One-sample grouped data
################################################
#' Mable fit of one-sample grouped data by an optimal or a preselected model degree
#' @param x vector of frequencies
#' @param breaks class interval end points
#' @param M a positive integer or a vector \code{(m0, m1)}. If \code{M = m} or \code{m0 = m1 = m},
#'   then \code{m} is a preselected degree. If \code{m0<m1} it specifies the set of
#'   consective candidate model degrees \code{m0:m1} for searching an optimal degree,
#'   where \code{m1-m0>3}.
#' @param interval a vector containing the endpoints of support/truncation interval
#' @param IC information criterion(s) in addition to Bayesian information criterion (BIC). Current choices are
#'  "aic" (Akaike information criterion) and/or
#'  "qhic" (Hannan–Quinn information criterion).
#' @param controls Object of class \code{mable.ctrl()} specifying iteration limit
#' and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.
#' @param progress if TRUE a text progressbar is displayed
#' @description  Maximum approximate Bernstein/Beta likelihood estimation based on
#'  one-sample grouped data with an optimal selected by the change-point method among \code{m0:m1}
#'  or a preselected model degree \code{m}.
#' @details
#'  Any continuous density function \eqn{f} on a known closed supporting interval \eqn{[a, b]} can be
#'  estimated by Bernstein polynomial \eqn{f_m(x; p) = \sum_{i=0}^m p_i\beta_{mi}[(x-a)/(b-a)]/(b-a)},
#'  where \eqn{p = (p_0, \ldots, p_m)}, \eqn{p_i\ge 0}, \eqn{\sum_{i=0}^m p_i=1} and
#'  \eqn{\beta_{mi}(u) = (m+1){m\choose i}u^i(1-x)^{m-i}}, \eqn{i = 0, 1, \ldots, m},
#'  is the beta density with shapes \eqn{(i+1, m-i+1)}.
#'  For each \code{m}, the MABLE of the coefficients \code{p}, the mixture proportions, are
#'  obtained using EM algorithm. The EM iteration for each candidate \code{m} stops if either
#'  the total absolute change of the log likelihood and the coefficients of Bernstein polynomial
#'  is smaller than \code{eps} or the maximum number of iterations \code{maxit} is reached.
#'
#'  If \code{m0<m1}, an optimal model degree is selected as the change-point of the increments of
#'  log-likelihood, log likelihood ratios, for \eqn{m \in \{m_0, m_0+1, \ldots, m_1\}}. Alternatively,
#'  one can choose an optimal degree based on the BIC (Schwarz, 1978) which are evaluated at
#'  \eqn{m \in \{m_0, m_0+1, \ldots, m_1\}}. The search for optimal degree \code{m} is stoped if either
#'  \code{m1} is reached with a warning or the test for change-point results in a p-value \code{pval}
#'  smaller than \code{sig.level}.  The BIC for a given degree \code{m} is calculated as in
#'  Schwarz (1978) where the dimension of the model is \eqn{d=\#\{i: \hat p_i \ge \epsilon,
#'  i = 0, \ldots, m\} - 1} and a default \eqn{\epsilon} is chosen as \code{.Machine$double.eps}.
#' @return A list with components
#' \itemize{
#'   \item \code{m} the given/selected optimal degree by the method of change-point
#'   \item \code{p} the estimated \code{p} with degree \code{m}
#'   \item \code{mloglik}  the maximum log-likelihood at degree \code{m}
#'   \item \code{interval} supporting interval \code{(a, b)}
#'   \item \code{convergence} An integer code. 0 indicates successful completion
#'    (all the EM iterations are convergent and an optimal degree
#'     is successfully selected in \code{M}). Possible error codes are
#'    \itemize{
#'       \item 1, indicates that the iteration limit \code{maxit} had been
#'           reached in at least one EM iteration;
#'       \item 2, the search did not finish before \code{m1}.
#'     }
#'   \item \code{delta} the convergence criterion \code{delta} value
#'  }
#'  and, if \code{m0<m1},
#' \itemize{
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1}, if greater than \code{m0}, is the
#'      largest candidate when the search stoped
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{ic} a list containing the selected information criterion(s)
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @author Zhong Guan <zguan@iusb.edu>
#' @references Guan, Z. (2017) Bernstein polynomial model for grouped continuous data.
#'         \emph{Journal of Nonparametric Statistics}, 29(4):831-848.
#' @examples
#' \donttest{
#' ## Chicken Embryo Data
# load("chicken.embryo.rdata")
#'  data(chicken.embryo)
#'  a<-0; b<-21
#'  day<-chicken.embryo$day
#'  nT<-chicken.embryo$nT
#'  Day<-rep(day,nT)
#'  res<-mable.group(x=nT, breaks=a:b, M=c(2,100), interval=c(a, b), IC="aic",
#'     controls=mable.ctrl(sig.level=1e-6,  maxit=2000, eps=1.0e-7))
#'  op<-par(mfrow=c(1,2), lwd=2)
#'  layout(rbind(c(1, 2), c(3, 3)))
#'  plot(res, which="likelihood")
#'  plot(res, which="change-point")
#'  fk<-density(x=rep((0:20)+.5, nT), bw="sj", n=101, from=a, to=b)
#'  hist(Day, breaks=seq(a,b,  length=12), freq=FALSE, col="grey",
#'           border="white", main="Histogram and Density Estimates")
#'  plot(res, which="density",types=1:2, colors=1:2)
#'  lines(fk, lty=2, col=2)
#'  legend("topright", lty=c(1:2), c("MABLE", "Kernel"), bty="n", col=c(1:2))
#'  par(op)
#' }
# @useDynLib mable mable_optim_group .registration = TRUE
#' @keywords distribution  models  nonparametric  smooth univar
#' @concept Bernstein polynomial model
#' @seealso \code{\link{mable.ic}}
#' @export
mable.group<-function(x, breaks, M, interval=c(0, 1), IC=c("none", "aic", "hqic", "all"),
            controls = mable.ctrl(), progress=TRUE){
    IC <- match.arg(IC, several.ok=TRUE)
    N<-length(x)
    xName<-deparse(substitute(x))
    if(missing(M) || length(M)==0) stop("'M' is missing.\n")
    else if(length(M)==1) M<-c(M,M)
    else if(length(M)>=2) {
        M<-c(min(M), max(M))
    }
    k<-M[2]-M[1]
    if(k>0 && k<=3) stop("Too few candidate model degrees.")
    a<-interval[1]
    b<-interval[2]
    t<-(breaks-a)/(b-a)
    if(any(t<0) | any(t>1)) stop("'breaks' must be in 'interval'")
    convergence<-0
    del<-controls$eps
    eps<-c(del, .Machine$double.eps)
    if(k==0){
        llik<-0
        m<-M[1]
        p<-rep(1,m+1)/(m+1)
        ## Call C mable_em_group
        res<-.C("mable_em_group",
          as.integer(m), as.integer(x), as.integer(N), as.double(p), as.double(t),
          as.integer(controls$maxit), as.double(controls$eps),  as.double(llik),
          as.logical(convergence), as.double(del))
        ans<-list(p=res[[4]], m=m, mloglik=res[[8]], interval=c(a,b), convergence=!res[[9]], del=res[[10]])
    }
    else{
        lk<-rep(0, k+1)
        lr<-rep(0, k)
        bic<-rep(0,k+1)
        pval<-rep(0,k+1)
        level<-controls$sig.level
        chpts<-rep(0,k+1)
        #p<-rep(1, 2*M[2]+2)
        #optim<-c(0,0)
        p<-rep(1, M[2]+1)
        optim<-0
        ## Call C mable_optim_group
        res<-.C("mable_optim_group",
          as.integer(M), as.integer(N), as.double(p), as.double(t), as.integer(x),
          as.integer(controls$maxit), as.double(eps), as.double(lk),
          as.double(lr), as.integer(optim), as.logical(progress),
          as.integer(convergence), as.double(del), as.double(controls$tini),
          as.double(bic), as.double(pval), as.integer(chpts), as.double(level))
        M<-res[[1]]
        k<-M[2]-M[1]
        optim<-res[[10]]
        m<-optim+1
        lk<-res[[8]][1:(k+1)]
        mloglik=lk[optim-M[1]+1]
        bic<-res[[15]][1:(k+1)]
        n<-sum(x)
        ic<-list()
        ic$BIC<-bic
        if(!any(IC=="none")){
            d<-2*(lk-bic)/log(n)
            if(any(IC=="aic")|| any(IC=="all")){
                ic$AIC<-((log(n)-2)*lk+2*bic)/log(n)#lk-2*(lk-bic)/log(n)
                #aic<-lk-d-(d^2+d)/(n-d-1)
            }
            if(any(IC=="qhic")|| any(IC=="all")){
                ic$QHC<-lk-2*(lk-bic)*log(log(n))/log(n)}
        }
        ans<-list(p=res[[3]][1:m[1]],
            mloglik=mloglik, m=optim, lk=lk, lr=res[[9]][1:k], M=M, interval=c(a,b),
            convergence=res[[12]], del=res[[13]], ic=ic, pval=res[[16]][1:(k+1)],
            chpts=res[[17]][1:(k+1)]+M[1])
    }
    ans$xNames<-xName
    ans$data.type<-"grp"
    class(ans)<-"mable"
    return(ans)
}
##############################################
#' Mixture Beta Distribution
#' @description Density, distribution function, quantile function and
#' pseudorandom number generation for the Bernstein polynomial model,
#' mixture of beta distributions, with shapes \eqn{(i+1, m-i+1)}, \eqn{i = 0, \ldots, m},
#'  given mixture proportions \eqn{p = (p_0, \ldots, p_m)} and support \code{interval}.
#' @param x a vector of quantiles
#' @param u a vector of probabilities
#' @param n sample size
#' @param p a vector of \code{m+1} values. The \code{m+1} components of \code{p}
#'   must be nonnegative and sum to one for mixture beta distribution. See 'Details'.
#' @param interval support/truncation interval \code{[a, b]}.
#' @return A vector of \eqn{f_m(x; p)} or \eqn{F_m(x; p)} values at \eqn{x}.
#' \code{dmixbeta} returns the density, \code{pmixbeta} returns the cumulative
#'  distribution function, \code{qmixbeta} returns the quantile function, and
#' \code{rmixbeta}  generates pseudo random numbers.
#' @details
#'  The density of the mixture beta distribution on an interval \eqn{[a, b]} can be written as a
#'  Bernstein polynomial \eqn{f_m(x; p) = \sum_{i=0}^m p_i\beta_{mi}[(x-a)/(b-a)]/(b-a)},
#'  where \eqn{p = (p_0, \ldots, p_m)}, \eqn{p_i\ge 0}, \eqn{\sum_{i=0}^m p_i=1} and
#'  \eqn{\beta_{mi}(u) = (m+1){m\choose i}u^i(1-x)^{m-i}}, \eqn{i = 0, 1, \ldots, m},
#'  is the beta density with shapes \eqn{(i+1, m-i+1)}. The cumulative distribution
#' function is \eqn{F_m(x; p) = \sum_{i=0}^m p_i B_{mi}[(x-a)/(b-a)]}, where
#' \eqn{B_{mi}(u)}, \eqn{i = 0, 1, \ldots, m}, is the beta cumulative distribution function
#'  with shapes \eqn{(i+1, m-i+1)}. If \eqn{\pi = \sum_{i=0}^m p_i<1}, then \eqn{f_m/\pi}
#'   is a truncated desity on \eqn{[a, b]} with cumulative distribution function
#'  \eqn{F_m/\pi}. The argument \code{p} may be any numeric vector of \code{m+1}
#'  values when \code{pmixbeta()} and and \code{qmixbeta()} return the integral
#'  function \eqn{F_m(x; p)} and its inverse, respectively, and \code{dmixbeta()}
#'  returns a Bernstein polynomial \eqn{f_m(x; p)}.
#' @author Zhong Guan <zguan@iusb.edu>
#' @references
#' Bernstein, S.N. (1912), Demonstration du theoreme de Weierstrass fondee sur le calcul des probabilities,
#' Communications of the Kharkov Mathematical Society, 13, 1–2.
#'
#' Guan, Z. (2016) Efficient and robust density estimation using Bernstein type polynomials. \emph{Journal of Nonparametric Statistics}, 28(2):250-271.
#'
#' Guan, Z. (2017) Bernstein polynomial model for grouped continuous data. \emph{Journal of Nonparametric Statistics}, 29(4):831-848.
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
#' @importFrom stats dunif punif qunif runif uniroot rbeta
#' @export
dmixbeta<-function(x, p, interval=c(0, 1)){
    n<-length(x)
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    if(any(x<a) || any(x>b)) stop("Argument 'x' must be in 'interval'.")
    if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
    if(any(p<0)) warning("Argument 'p' has negative component(s).\n")
    m<-length(p)-1
    if(m==0) y<-dunif(x, a, b)
    else{
        cdf<-FALSE
        u<-(x-a)/(b-a)
        res<-.C("mable_approx",
          as.double(u), as.double(p), as.integer(m), as.integer(n), as.logical(cdf))
        y<-res[[1]]/(b-a)
    }
    return(y)
}
#' @rdname dmixbeta
#' @export
pmixbeta<-function(x, p, interval=c(0, 1)){
    n<-length(x)
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    if(any(x<a) || any(x>b)) stop("Argument 'x' must be in 'interval'.")
    if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
    if(any(p<0)) warning("Argument 'p' has negative component(s).\n")
    m<-length(p)-1
    if(m==0) y<-punif(x, a, b)
    else{
        cdf<-TRUE
        u<-(x-a)/(b-a)
        res<-.C("mable_approx",
          as.double(u), as.double(p), as.integer(m), as.integer(n), as.logical(cdf))
        y<-res[[1]]
    }
    return(y)
}
#' @rdname dmixbeta
#' @export
#########################################################
# another method is Newton-Raphson as 'mable-roc.r'
#########################################################
qmixbeta<-function(u, p, interval=c(0, 1)){
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    if(any(u<0) || any(u>1)) stop("Argument 'u' must be in [0,1].\n")
    if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
    if(any(p<0)) warning("Argument 'p' has negative component(s).\n")
    m<-length(p)-1
    if(m==0) Q<-qunif(u, a, b)
    else{
        n<-length(u)
        Q<-u
        for(i in 1:n){
            fn<-function(x) pmixbeta(x, p)-u[i]
            Q[i]<-uniroot(fn, c(0,1))$root
        }
        Q<-a+(b-a)*Q
    }
    return(Q)
}
#' @rdname dmixbeta
#' @export
#########################################################################
# generating prn from sum(p[i]*beta[m,i][(x-a)/(b-a)]/(b-a): i=0,...,m)
#########################################################################
rmixbeta<-function(n, p, interval=c(0, 1)){
    a<-interval[1]
    b<-interval[2]
    if(a>=b) stop("a must be smaller than b")
    if(length(p)==0) stop("Missing mixture proportions 'p' without default.")
    if(any(p<0)) stop("Negative component(s) of argument 'p'is not allowed.\n")
    if(abs(sum(p)-1)>.Machine$double.eps){
        warning("Sum of 'p's is not 1. Dividing 'p's by the total.\n")
            p<-p/sum(p)
    }
    m<-length(p)-1
    x<-rep(0,n)
    if(m==0) x<-runif(n, a, b)
    else{
        w<-sample(0:m, n, replace = TRUE, prob = p)
        x<-rbeta(n, shape1 = w+1, shape2 = m-w+1)
        x<-a+(b-a)*x
        #res<-.C("rbeta_mi", as.integer(n), as.integer(m), as.integer(w), as.double(x))
        #x<-a+(b-a)*res[[4]]
    }
    return(x)
}
##############################################
#' Plot mathod for class 'mable'
##############################################
#' @param x  Class "mable" object return by \code{mablem}, \code{mable}, \code{mablem.group} or \code{mable.group} functions
#'      which contains \code{p}, \code{mloglik}, and \code{M = m0:m1}, \code{lk}, \code{lr},
#' @param which indicates which graphs to plot, options are
#'  "density", "cumulative", "likelihood", "change-point", "all". If not "all",
#'  \code{which} can contain more than one options.
#' @param add logical add to an existing plot or not
#' @param lgd.x,lgd.y coordinates of position where the legend is displayed
#' @param nx  number of evaluations of density, or cumulative distribution curve to be plotted.
#' @param ...  additional arguments to be passed to the base plot function
#' @importFrom graphics axis lines plot segments points title legend par persp
#' @export
plot.mable<-function(x, which=c("density", "cumulative", "survival", "likelihood",
                "change-point", "all"), add = FALSE,
                lgd.x=NULL, lgd.y=NULL, nx=512,...){
    phat<-x$p
    m<-x$m
    if(is.null(x$dim)) dim<-1
    else dim<-x$dim
    if(dim>2) stop("There is no method to 'plot' the object.\n")
    support<-x$interval
    if(is.null(lgd.x)) lgd.x="topright"
    which <- match.arg(which, several.ok=TRUE)
    if(dim==1){
        a<-support[1]
        b<-support[2]
        xx<-seq(a, b,len=nx)
        if(any(which=="all")){
            add=FALSE
            op<-par(mfrow=c(2,2))}
    if(any(which=="likelihood")|| any(which=="all")) {
        if(is.null(x$lk)) stop("Cannot plot 'likelihood'.")
        M<-x$M[1]:x$M[2]
        if(!add) plot(M, x$lk, type="p", xlab="m",ylab="Loglikelihood",
            main="Loglikelihood")
        else points(M, x$lk, pch=1)
        segments(x$m, min(x$lk), x$m, x$mloglik, lty=2, col=1:2)
        axis(1, x$m,  as.character(x$m))
    }
    if(any(which=="change-point")||any(which=="all")) {
        if(is.null(x$lr)) stop("Cannot plot 'likelihood ratios'.")
        M<-x$M[1]:x$M[2]
        ymin<-min(x$lr)
        ymax<-max(x$lr)
        if(!add){
            plot(M, 0*M, type="n", xlab="m", ylab="", ylim=c(ymin,ymax))
            points(M[-1], x$lr, col=1, pch=1, ylab="Likelihood Ratio")
            segments(m, ymin, m, max(x$lr), lty=2, col=1)
            axis(1, m,  as.character(m), col=1)
            title("LR of Change-Point")}
        else{
            points(M[-1], x$lr, ...)
            segments(m, ymin, m, max(x$lr), ...)
            axis(1, m,  as.character(m), ...)}
    }
    if(any(which=="density")||any(which=="all")){
        yy<-dmixbeta(xx, phat[1:(m+1)], c(a, b))
        if(!add)
            plot(xx, yy, type="l",  ylab="Density", xlim=c(a, b), ...)
        else
            lines(xx, yy, ...)
    }
    if(any(which=="cumulative")||any(which=="survival")||any(which=="all")){
        if(any(which=="survival")||any(which=="all"))
            yy<-1-pmixbeta(xx, phat[1:(m+1)], c(a, b))
        else
            yy<-pmixbeta(xx, phat[1:(m+1)], c(a, b))
        if(!add) {
            if(any(which=="survival")||any(which=="all")) plot(xx, yy,
                type="l", ylab="Survival Probability", xlim=c(a,b), col=1, ...)
            else plot(xx, yy, type="l", ylab="Cumulative Probability", xlim=c(a,b), col=1, ...)}
        else lines(xx, yy, ...)}
        if(any(which=="all")) par(op)
    }
    else{
        a<-support[,1]
        b<-support[,2]
        if(!any(which=="cumulative")&& !any(which=="density")&& !any(which=="all")){
            cat("Only 'density' and 'cumulative' distribution can be plotted.\n")
            which<-"all"
        }
        if(any(which=="all")) op<-par(mfrow=c(2,1))
        if(any(which=="density")||any(which=="all")){
            fn<-function(x1, x2) dmixmvbeta(cbind(x1, x2), phat, m, interval=support)
            x1 <- seq(a[1], b[1], length= 40)
            x2 <- seq(a[2], b[2], length= 40)
            z1 <- outer(x1, x2, fn)
            persp(x1, x2, z1, theta = 30, phi = 20, expand = 0.5, col = "lightblue",
              ltheta = 90, shade = 0.1, ticktype = "detailed", main = expression(paste("MABLE ",hat(f))),
              xlab = x$xNames[1], ylab = x$xNames[2], zlab = "Joint Density")
        }
        if(any(which=="cumulative")||any(which=="all")){
            fn<-function(x1, x2) pmixmvbeta(cbind(x1, x2), phat, m, interval=support)
            x1 <- seq(a[1], b[1], length= 40)
            x2 <- seq(a[2], b[2], length= 40)
            z2 <- outer(x1, x2, fn)
            persp(x1, x2, z2, theta = 30, phi = 20, expand = 0.5, col = "lightblue",
              ltheta = 90, shade = 0.1, ticktype = "detailed", main = expression(paste("MABLE ",hat(F))),
              xlab = x$xNames[1], ylab = x$xNames[2], zlab = "Joint CDF")
        }
        if(any(which=="all")) par(op)
    }
}
##############################################
#' Summary mathods for classes 'mable' and 'mable_reg'
##############################################
#' @param object  Class "mable" or 'mable_reg' object return by \code{mable} or \code{mable.xxxx}  functions
#' @param ...  for future methods
#' @description Produces a summary of a mable fit.
#' @return Invisibly returns its argument, \code{object}.
#' @examples
#' \donttest{
#'   # Vaal Rive Flow Data
#'   data(Vaal.Flow)
#'   res<-mable(Vaal.Flow$Flow, M = c(2,100), interval = c(0, 3000),
#'      controls = mable.ctrl(sig.level = 1e-8, maxit = 2000, eps = 1.0e-9))
#'   summary(res)
#' }
#' \donttest{
#' ## Breast Cosmesis Data
#'   require(interval)
#'   data(bcos)
#'   bcos2<-data.frame(bcos[,1:2], x=1*(bcos$treatment=="RadChem"))
#'   aft.res<-mable.aft(cbind(left, right)~x, data=bcos2, M=c(1, 30), tau=100, x0=1)
#'   summary(aft.res)
#' }
#' @importFrom stats printCoefmat
#' @export
summary.mable<-function(object, ...){
    obj<-object
    cl<-class(obj)
    ans<-list()
    if(obj$data.type=="mvar") ans$m<-obj$m
    else ans$m<-obj$m
    ans$mloglik<-obj$mloglik
    ans$interval<-obj$interval
    if(is.null(obj$dim)) ans$dim<-1
    else  ans$dim<-obj$dim
    ans$p<-obj$p
    ans$pval<-obj$pval[length(obj$pval)]
    switch(obj$data.type,
        raw=cat("Call: mable() for raw data"),
        grp=cat("Call: mable.group() for grouped data"),
        mvar=cat("Call: mable.mvar() for multivariate data"),
        icen=cat("Call: mable.ic() for interval censored data"),
        noisy=cat("Call: mable.decon() for noisy data"))
    cat("\nObj Class:", cl,"\n")
    cat("Input Data:", obj$xNames,"\n")
    cat("Dimension of data:", ans$dim,"\n")
    if(obj$data.type=="mvar"){
        cat("Optimal degrees:\n")
        prnt<-matrix(ans$m, nrow=1)
        rownames(prnt)<-"m"
        colnames(prnt)<-obj$xNames
        printCoefmat(prnt)
    }
    else cat("Optimal degree m:", ans$m,"\n")
    cat("P-value of Change-point:", ans$pval,"\n")
    cat("Maximum loglikelihood:", ans$mloglik,"\n")
    cat("MABLE of p: can be retrieved using name 'p' \n")
    cat("Note: the optimal model degrees are selected by the method of change-point.\n")
    invisible(ans)
}
#' @rdname summary.mable
#' @export
summary.mable_reg<-function(object, ...){
    obj<-object
    cl<-class(obj)
    ans<-list()
    ans$m<-obj$m
    ans$mloglik<-obj$mloglik
    ans$interval<-obj$interval
    ans$baseline<-obj$x0
    ans$names<-obj$xNames
    ans$coefficients<-obj$coefficients
    ans$dim<-1
    ans$dimcov<-length(obj$x0)
    ans$se<-obj$SE
    ans$z<-obj$z
    ans$p<-obj$p
    ans$pval<-obj$pval[length(obj$pval)]
    cat(paste("Call: ",obj$model,"(",obj$callText,")", sep=''))
    cat("\nData:", obj$data.name,"\n")
    cat("Obj Class:", cl,"\n")
    cat("Dimension of response:", ans$dim,"\n")
    cat("Dimension of covariate:", ans$dimcov,"\n")
    cat("Optimal degree m:", ans$m,"\n")
    cat("P-value:", ans$pval,"\n")
    cat("Maximum loglikelihood:", ans$mloglik,"\n")
    cat("MABLE of p: can be retrieved using name 'p' \n")
    cat("\n")
    prnt<-cbind(ans$coefficients, ans$se, ans$z, ans$baseline)
    colnames(prnt)<-c("Estimate", "Std.Err", "Z value", "Baseline x0")
    rownames(prnt)<-ans$names
    printCoefmat(prnt)
    cat("\nNote: the optimal model degree is selected by the method of change-point.\n")
    invisible(ans)
}
####################################################
#' Choosing optimal model degree by gamma change-point method
#' @description Choose an optimal degree using gamma change-point model with two
#'   changing shape and scale parameters.
#' @param obj a class "mable" or 'mable_reg' object containig a vector \code{M = (m0, m1)}, \code{lk},
#'    loglikelihoods evaluated evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#' @return a list with components
#' \itemize{
#'   \item \code{m} the selected optimal degree \code{m}
#'   \item \code{M} the vector \code{(m0, m1)}, where \code{m1} is the last candidate when the search stoped
#'   \item \code{mloglik}  the maximum log-likelihood at degree \code{m}
#'   \item \code{interval} support/truncation interval \code{(a, b)}
#'   \item \code{lk} log-likelihoods evaluated at \eqn{m \in \{m_0, \ldots, m_1\}}
#'   \item \code{lr} likelihood ratios for change-points evaluated at \eqn{m \in \{m_0+1, \ldots, m_1\}}
#'   \item \code{pval} the p-values of the change-point tests for choosing optimal model degree
#'   \item \code{chpts} the change-points chosen with the given candidate model degrees
#' }
#' @examples \donttest{
#'  # simulated data
#'  p<-c(1:5,5:1)
#'  p<-p/sum(p)
#'  x<-rmixbeta(100, p)
#'  res1<-mable(x, M=c(2, 50), IC="none")
#'  m1<-res1$m[1]
#'  res2<-optim.gcp(res1)
#'  m2<-res2$m
#'  op<-par(mfrow=c(1,2))
#'  plot(res1, which="likelihood", add=FALSE)
#'  plot(res2, which="likelihood")
#'  #segments(m2, min(res1$lk), m2, res2$mloglik, col=4)
#'  plot(res1, which="change-point", add=FALSE)
#'  plot(res2, which="change-point")
#'  par(op)
#' }
#' @export
optim.gcp<-function(obj){
    M<-obj$M
    lk<-obj$lk
    k<-M[2]-M[1]
    lr<-rep(0, k)
    pval<-rep(0, k+1)
    chpts<-rep(0, k+1)
    m<-M[1]
    res<-.C("optim_gcp", as.integer(M), as.double(lk), as.double(lr), as.integer(m),
        as.double(pval), as.integer(chpts))
    ans<-list(m=res[[4]], M=M, mloglik=lk[res[[4]]-M[1]+1], lk=lk, lr=res[[3]],
        interval=obj$interval, pval=res[[5]], chpts=res[[6]]+M[1])
    class(ans)<-"mable"
    return(ans)
}
