\name{mable.em}
\alias{mable.em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit raw data with the Bernstein polynomial model with a given degree m
%%  ~~function to do ... ~~
}
\description{It returns the MABLE of \eqn{p} and the log-likelihood.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{mable.em(x, m, a = 0, b = 1, maxit = 50000, eps = 1e-09)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{raw sample data
%%     ~~Describe \code{x} here~~
}
  \item{m}{a given model degree
%%     ~~Describe \code{m} here~~
}
  \item{a}{left endpoint of supporting interval
%%     ~~Describe \code{a} here~~
}
  \item{b}{right endpoint of supporting interval
%%     ~~Describe \code{b} here~~
}
  \item{maxit}{maximum number of iterations
%%     ~~Describe \code{maxit} here~~
}
  \item{eps}{convergence criterion for iteration
%%     ~~Describe \code{eps} here~~
}
}
\details{It returns the MABLE of \eqn{p} and the log-likelihood. The MABLE of the coefficients \eqn{p}, the mixture proportions, are obtained using EM algorithm.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{phat }{MABLE of \eqn{p}, the coefficients of Bernstein polynomial}
  \item{llik }{log-likelihood}
  \item{support }{supporting interval (a,b)}
%% ...
}
\references{
Guan, Z. (2016) Efficient and robust density estimation using Bernstein type
              polynomials. \emph{Journal of Nonparametric Statistics}, 28(2):250-271.
%% ~put references to the literature/web site here ~
}
\author{Zhong Guan <zguan@iu.edu>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mable.optim}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Old Faithful Data
  library(mixtools)
  x<-faithful$eruptions
  a<-0;  b<-7
  u<-seq(0,1,len=512); v<-(b-a)*u+a
  mu<-c(2,4.5); sig<-c(1,1)
  pmixem<-normalmixEM(x,.5,mu, sig)
  y1<-pmixem$lambda[1]*dnorm(v,pmixem$mu[1], pmixem$sigma[1])
      +pmixem$lambda[2]*dnorm(v,pmixem$mu[2],pmixem$sigma[2])
  mhat<-94 # preselected optimal degree m
  res<-mable.em(x, m=mhat, a, b, maxit=2000, eps=1.0e-4)
  hist(x, breaks=seq(0,7.5,len=20), xlim=c(0,7), ylim=c(0,.7), 
      probability =TRUE, xlab="t", ylab="f(t)", col ="light grey", 
      main="Histogram and Density of Duration of Old Faithful")
  lines(density(x, bw = "nrd0", adjust = 1), lty=4, col = 4, lwd=2)
  lines(v, bern.approx(u, p=res$phat)/(b-a), lty=1, col = 1, lwd=2)
  lines(v, y1, lty=5, col=3, lwd=2)
  legend(6,.4, lty=c(5,1,4), col=c(3,1,4), lwd=2, bty="n", 
      c(expression(hat(f)[P](t)), expression(hat(f)[B](t)), 
      expression(hat(f)[K](t))))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~models }
\keyword{ ~nonparametric }
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
