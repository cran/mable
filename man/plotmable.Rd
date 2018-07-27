\name{plotmable}
\alias{plotmable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Maximum Approximate Bernstein Likelihood Estimates of \eqn{f} or \eqn{F}}
\description{Takes an output of \code{mable.em()} or \code{mable.optim()} and returns various graphical output for the nonparametric estimate.}
\usage{
plotmable(mable.fit, density=TRUE, nx=512, add = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mable.fit}{mable object containing \eqn{\hat p} and (a,b) of \code{mable.em()} or \code{mable.optim()}.}
  \item{density}{logical; whether a density or a distribution function to be calculated. Default is TRUE.}
  \item{nx}{number of points to plot.}
  \item{add}{logical; if TRUE only add to an existing plot. Default is TRUE.}
  \item{...}{other arguments of \code{plot()}.}
}
\references{
Guan, Z. (2016) Efficient and robust density estimation using Bernstein type
              polynomials. \emph{Journal of Nonparametric Statistics}, 28(2): 250-271.

Guan, Z. (2017) Bernstein polynomial model for grouped continuous data.
              \emph{Journal of Nonparametric Statistics}, 29(4): 831-848.
%% ~put references to the literature/web site here ~
}
\author{Zhong Guan <zguan@iu.edu>
%%  ~~who you are~~
}
 %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mable.em}}, \code{\link{mable.optim}}, \code{\link{mable.em.group}},
\code{\link{mable.optim.group}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
# standard normal
a<--4; b<-4;
x<-rnorm(30)
bfit<-mable.em(x, m=20, a, b)
# density estimate
plot(xx<-seq(a,b,len=512), dnorm(xx), type="l", xlab="x", ylab="f(x)")
plotmable(bfit, lty=2,col=2)
# cdf estimate
plot(xx<-seq(a,b,len=512), pnorm(xx), type="l", xlab="x", ylab="F(x)")
plotmable(bfit, density=FALSE, lty=2,col=2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~models }
\keyword{ ~nonparametric }
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
