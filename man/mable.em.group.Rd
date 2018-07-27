\name{mable.em.group}
\alias{mable.em.group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit grouped data with the Bernstein polynomial model with a given degree m
%%  ~~function to do ... ~~
}
\description{It returns the MABLE of \eqn{p} and the log-likelihood.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{mable.em.group(x, breaks, m, a=0, b=1, maxit = 50000, 
                       eps = 1.0e-9)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of frequencies
%%     ~~Describe \code{x} here~~
}
  \item{breaks}{end points of class intervals that  partition [0,1] after transformation
%%     ~~Describe \code{m} here~~
}
  \item{m}{a given model degree \eqn{m}
%%     ~~Describe \code{a} here~~
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
\details{The MABLE of the coefficients \eqn{p}, the mixture proportions, are obtained using EM algorithm.}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{phat }{MABLE of \eqn{p}, the coefficients of Bernstein polynomial}
  \item{llik }{log-likelihood}
  \item{support }{supporting interval (a,b)}
%% ...
}
\references{
Guan, Z. (2017) Bernstein polynomial model for grouped continuous data.
              \emph{Journal of Nonparametric Statistics}, 29(4):831-848.
%% ~put references to the literature/web site here ~
}
\author{Zhong Guan <zguan@iu.edu>
%%  ~~who you are~~
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ }
\examples{
# Chicken Embryo Data
data(chicken.embryo)
a<-0; b<-21
t<-((a:b)-a)/(b-a)
day<-chicken.embryo$day
nT<-chicken.embryo$nT
Day<-rep(day,nT)
res<-mable.em.group(x=nT, t, m=13)
phat<-res$phat; mhat<-13
xx<-seq(0,1,len=100)
fb<-bern.approx(x=xx, p=phat)/(b-a)
fk<-density(x=rep((0:20)+.5, nT), bw="sj", n=101, from=a, to=b)
hist(Day, breaks=seq(a,b,  length=12), freq=FALSE,
         main="Histogram and Density Estimates")
lines(a+(b-a)*xx, fb, lty=2, col=2)
lines(fk, lty=3, col=3)
legend(14, .2, lty=1:3, c("Histogram", "MABLE", "Kernel"),  
                           bty="n", col=1:3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~models }
\keyword{ ~nonparametric }
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
