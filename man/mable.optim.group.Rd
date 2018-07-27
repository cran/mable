\name{mable.optim.group}
\alias{mable.optim.group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit grouped data with the Bernstein polynomial model with an optimal degree m
%%  ~~function to do ... ~~
}
\description{It also returns the MABLE of \eqn{p} and the log-likelihood.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mable.optim.group(x, breaks, m0, m1, a=0, b=1,
                            maxit = 50000, eps = 1e-09)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{vector of frequencies
%%     ~~Describe \code{x} here~~
}
  \item{breaks}{end points of class intervals that  partition [0,1] after transformation
%%     ~~Describe \code{t} here~~
}
  \item{m0}{specify a starting point for choosing an optimal model degree
%%     ~~Describe \code{m0} here~~
}
  \item{m1}{specify an end point for choosing an optimal model degree
%%     ~~Describe \code{m1} here~~
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
\details{An optimal model degree is selected as the change-point of the increments of log-likelihood, log likelihood ratios, for \eqn{m\in\{m_0,m_0+1,\ldots,m_1\}}.
For each \eqn{m},  the MABLE of the coefficients \eqn{p}, the mixture proportions, are obtained using EM algorithm.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Llik }{log-likelihoods evaluated at \eqn{m\in\{m_0, \ldots, m_1\}}}
  \item{phat }{MABLE of \eqn{p}, the coefficients of Bernstein polynomial of the selected optimal degree \eqn{m}}
  \item{llik }{log-likelihood at optimal degree \eqn{m}}
  \item{opti.m }{the selected optimal degree \eqn{m}}
  \item{lr }{likelihood ratios for change-points evaluated at \eqn{m\in\{m_0, \ldots, m_1\}}}
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
%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ }
\examples{
\donttest{
# Chicken Embryo Data
data(chicken.embryo)
m0<-2; m1<-30; a<-0; b<-21
t<-((a:b)-a)/(b-a);
day<-chicken.embryo$day; nT<-chicken.embryo$nT;
Day<-rep(day,nT)
res<-mable.optim.group(x=nT, t, m0, m1)
p<-res$phat; m<-res$opti.m
op<-par(mfrow=c(1,2), lwd=2)
layout(rbind(c(1, 2), c(3, 3)))
plot(m0:m1, Lk<-res$Llik, type="b", xlab="m",
      ylab=expression(paste("Loglikelihook ",  "\u2113", "(m)")))
segments(m,-130, m, Lk[m], lty=2, col=2)
text(m+5,-119, expression(paste(hat(m)==13)),pos=4, cex=1.5)
arrows(m+5.5, -119.5, m,-121.4, length =.07, angle =20, code=2)
plot((m0+1):m1, lr<-res$lr, type="b", xlab="m", ylim=c(0, max(lr)),
     ylab="Likelihood Ratio: R(m)")
segments(m,-10, m, lr[m-2], lty=2, col=2)
text(m+3,5, expression(hat(m)==13),pos=4, cex=1.5)
arrows(m+3.45, 4.3, m,-.7, length = .07, angle = 20, code = 2)
u<-seq(0,1,len=100); y<-a+(b-a)*u
fb<-bern.approx(u, p=res$phat)/(b-a)
fk<-density(x=rep((0:20)+.5, nT), bw="sj", n=101, from=a, to=b)
hist(Day, breaks=seq(a,b,  length=12), freq=FALSE,
         main="Histogram and Density Estimates")
lines(y, fb, lty=2, col=2)
lines(fk, lty=3, col=3)
legend(14, .2, lty=1:3, c("Histogram", "MABLE", "Kernel"),
               bty="n", col=1:3)
par(op)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~models }
\keyword{ ~nonparametric }
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
