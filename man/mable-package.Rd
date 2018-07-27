\name{mable-package}
\alias{mable-package}
\alias{mable}
\docType{package}
\title{MABLE: Maximum Approximate Bernstein Likelihood Estimation}
\description{
%More about what it does (maybe more than one line)
%~~ A concise (1-5 lines) description of the package ~~
Fit raw or grouped data from a nonparametric density \eqn{f} on \eqn{[0,1]} by an approximate Bernstein polynomial model
\eqn{f_m(x;p)=\sum_{i=0}^m p_i\beta_{mi}(x)}, where \eqn{\beta_{mi}} is the density of beta distribution with shapes \eqn{(i+1, m-i+1)}, \eqn{i=0,\ldots,m}, and find maximum approximate Bernstein likelihood estimator of \eqn{p=(p_0,\ldots,p_m)}. If the support of \eqn{f} is \eqn{[a,b]} then the data \eqn{x_i} are transformed to \eqn{(x_i-a)/(b-a)}.
}
\details{
\tabular{ll}{
Package: \tab mable\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2018-07-11\cr
License: GPL (>= 3) \cr
}
Functions: \cr
\code{\link{bern.approx}:} compute the Bernstein polynomial approximate of density or distribution function.\cr
\code{\link{mable}:} fit raw data by the approximate Bernstein polynomial model with a given degree \eqn{m}.\cr
\code{\link{mable.optim}:} choose an optimal model degree \eqn{m} and find the maximum approximate Bernstein likelihood estimator of \eqn{p} based on raw data.\cr
\code{\link{mable.optim.group}:} choose an optimal model degree \eqn{m} and find the maximum approximate Bernstein likelihood estimator of \eqn{p} based on grouped data.\cr
\code{\link{plotmable}:} plot the maximum approximate Bernstein likelihood estimate of density function \eqn{f} or distribution function \eqn{F}.\cr
}
\author{Zhong Guan <zguan@iu.edu>
}
\references{
Guan, Z. (2016) Efficient and robust density estimation using Bernstein type
              polynomials. \emph{Journal of Nonparametric Statistics}, 28(2):250-271.
              https://doi.org/10.1080/10485252.2016.1163349

Guan, Z. (2017) Bernstein polynomial model for grouped continuous data.
              \emph{Journal of Nonparametric Statistics}, 29(4):831-848.
              https://doi.org/10.1080/10485252.2017.1374384
}
\keyword{ package }
%\seealso{}
\examples{
\donttest{
# Beta distribution
n<-50; m0<-2; m1<-20
set.seed(1234567)
x<-rbeta(n,7,5)
res<-mable.optim(x, m0, m1)
op<-par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(m0:m1, res$Llik, type="b", xlab="m", ylab="Loglikelihook l(m)")
segments(m<-res$opti.m, 0, m, res$llik, lty=2, col=2)
plot((m0+1):m1, lr<-res$lr, type="b", xlab="m",
                         ylab="LR for Change-Point: R(m)")
segments(m, 0, m, max(lr), lty=2, col=2)
plot(y<-seq(0,1,len=512), dbeta(y,7,5), type="l", xlab="x",
            ylab="f(x)", ylim=c(0, 3.1))
plotmable(res, lty=2, col=2)
lines(density(x, n=512, from=0, to=1), lty=4, col=4)
legend(0, 2.7, bty="n", cex=.8, lty=c(1,2,4),col=c(1,2,4),
           c("TRUE: f",expression(paste("MABLE:",hat(f))),
                 expression(paste("KDE:",tilde(f)))))
plot(y, pbeta(y,7,5), type="l", xlab="x", ylab="F(x)",
                main="CDF of beta(7,5)")
plotmable(res, density=FALSE, lty=2,col=2)
lines(ecdf(x), do.points =FALSE, verticals = TRUE, lty=4, col=4)
legend(0, 1, bty="n", cex=.8, lty=c(1,2,4),col=c(1,2,4),
     c("TRUE: F", expression(paste("MABLE:",hat(F))),
                            expression(paste("ECDF:",F[n]))))
par(op)

# Normal distribution
n<-30; m0<-2; m1<-30
set.seed(13579)
x<-rnorm(n); a<--4; b<-4
res<-mable.optim(x, m0, m1, a, b)
op<-par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(m0:m1, res$Llik, type="b", xlab="m", ylab="Loglikelihook l(m)")
segments(m<-res$opti.m, 0, m, res$llik, lty=2, col=2)
plot((m0+1):m1, lr<-res$lr, type="b", xlab="m",
                              ylab="LR for Change-Point:  R(m)")
segments(m, 0, m, max(lr), lty=2, col=2)
plot(y<-seq(a,b,len=512), dnorm(y), type="l", xlab="x", ylab="f(x)")
plotmable(res, lty=2, col=2)
lines(density(x, n=512, from=a, to=b), lty=4, col=4)
legend(-4, .4, bty="n", cex=.8, lty=c(1,2,4),col=c(1,2,4),
    c("TRUE: f", expression(paste("MABLE:",hat(f))),
                       expression(paste("KDE:",tilde(f)))))
plot(y, pnorm(y), type="l", xlab="x", ylab="F(x)")
plotmable(res, density=FALSE, lty=2, col=2)
lines(ecdf(x), do.points =FALSE, verticals = TRUE, lty=4, col=4)
legend(-4, 1, bty="n", cex=.8, lty=c(1,2,4),col=c(1,2,4),
    c("TRUE: F", expression(paste("MABLE:",hat(F))),
                             expression(paste("ECDF:",F[n]))))
par(op)
}
}
