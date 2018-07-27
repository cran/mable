\name{mable.optim}
\alias{mable.optim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit raw data with the Bernstein polynomial model with optimal degree m
%%  ~~function to do ... ~~
}
\description{It also returns the MABLE of \eqn{p} and the log-likelihood.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mable.optim(x, m0, m1, a = 0, b = 1, maxit = 50000,
                  eps = 1e-09, tini = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{raw sample data
%%     ~~Describe \code{x} here~~
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
  \item{tini}{a tiny positive number to keep initial \eqn{p} in the interior of the simplex in each iteration
%%     ~~Describe \code{tini} here~~
}
}
\details{An optimal model degree is selected as the change-point of the increments of log-likelihood, log likelihood ratios, for \eqn{m\in\{m_0,m_0+1,\ldots,m_1\}}.
For each \eqn{m},  the MABLE of the coefficients \eqn{p}, the mixture proportions, are obtained using EM algorithm.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Llik }{log-likelihoods evaluated at \eqn{m\in\{m_0,\ldots, m_1\}}}
  \item{phat }{MABLE of \eqn{p}, the coefficients of Bernstein polynomial of the selected optimal degree \eqn{m}}
  \item{llik }{log-likelihood at optimal degree \eqn{m}}
  \item{opti.m }{the selected optimal degree \eqn{m}}
  \item{lr }{likelihood ratios for change-points evaluated at \eqn{m\in\{m_0,\ldots, m_1\}}}
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
\note{Since the Bernstein polynomial model of degree \eqn{m} is nested in the model of degree \eqn{m+1}, the maximum likelihood is increasing in \eqn{m}. The change-point method is used to choose an optimal degree \eqn{m}.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mable.em}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
# Vaal Rive Flow Data
data(Vaal.Flow)
x<-Vaal.Flow
a<-0
b<-3000
res<-mable.optim(x, m0=2, m1=40, a, b)
m<-res$opti.m
p<-res$phat
op<-par(mfrow=c(1,2))
layout(rbind(c(1, 2), c(3, 3)))
plot(2:40, Lk<-res$Llik, type="b", pch=20,  xlab="m",
                      ylab="Loglikelihook  \u2113(m)")
segments(m,-10, m, Lk[m], lty=2, col=2)
text(m+5.3,45, expression(paste(hat(m)==19)),pos=4)
arrows(m+5, 44, m,41, length = .07, angle = 20, code = 2)
plot(3:40, lr<-res$lr, type="b", pch=20, xlab="m",
                    ylab="LR for Change-Point: R(m)")
segments(m,-1, m, lr[m-2], lty=2, col=2)
text(m+5.1,3.5, expression(hat(m)==19),pos=4)
arrows(m+5, 3.0, m,-1, length = .07, angle = 20, code = 2)
op<-par(lwd=2)
hist(x, prob=TRUE, xlim=c(a,b), ylim=c(0,.0022), breaks=100*(0:30),
  main="Histogram and Densities of the Annual Flow of Vaal River",
border="dark grey",lwd=1,xlab="x", ylab="f(x)", col ="light grey")
lines(density(x, bw = "nrd0", adjust = 1), lty=4, col = 4)
lines(y<-seq(a, b, length=100), dlnorm(y, mean(log(x)),
                   sqrt(var(log(x)))), lty=2, col=2)
lines(y, bern.approx((y-a)/(b-a), p)/(b-a), col=1)
legend(1500, .0015, lty=c(1,2, 4), col=c(1,2, 4), bty="n",
c(expression(paste("MABLE:",hat(f)[B](x))),
        expression(paste("Log-Normal:",hat(f)[P](x))),
               expression(paste("KDE:",hat(f)[K](x)))))
par(op)
}
\dontrun{
# Old Faithful Data
  library(mixtools)
  x<-faithful$eruptions
  a<-0; b<-7
  u<-seq(0,1,len=512); v<-(b-a)*u+a
  mu<-c(2,4.5); sig<-c(1,1)
  pmixem<-normalmixEM(x,.5,mu, sig)
  y1<-pmixem$lambda[1]*dnorm(v,pmixem$mu[1], pmixem$sigma[1])
      +pmixem$lambda[2]*dnorm(v,pmixem$mu[2],pmixem$sigma[2])
  x0<-(x-a)/(b-a); m0<-2; m1<-300
  res<-mable.optim(x0, m0, m1, maxit=2000, eps=1.0e-4)
  m<-res$opti.m; p<-res$phat
  op<-par(mfrow=c(1,2))
  layout(rbind(c(1, 2), c(3, 3)))
  plot(m0:m1, res$Llik, type="b", pch=20,  xlab="m",
       ylab="Loglikelihook  \u2113(m)")
  segments(m,-10, m, res$Llik[m], lty=2, col=2)
  text(m+20,100, expression(paste(hat(m)==94)),pos=4)
  arrows(m+24, 95, m,70, length = .07, angle = 20, code = 2)
  plot((m0+1):m1, res$lr, type="b", pch=20, xlab="m",
                           ylab="Likelihood Ratio:  R(m)")
  segments(m-2,-10, m-2, res$lr[m-2], lty=2, col=2)
  text(m+20,40, expression(hat(m)==94),pos=4)
  arrows(m+20, 30, m,-7, length = .07, angle = 20, code = 2)
  hist(x, breaks=seq(0,7.5,len=20), xlim=c(0,7), ylim=c(0,.7),
     prob =TRUE,xlab="t", ylab="f(t)", col ="light grey",
     main="Histogram and Density of
               Duration of Eruptions of Old Faithful")
  lines(density(x, bw = "nrd0", adjust = 1), lty=4, col = 4, lwd=2)
  lines(v, bern.approx(u, p)/(b-a), lty=1, col = 1, lwd=2)
  lines(v, y1, lty=5, col=3, lwd=2)
  legend(6,.4, lty=c(5,1,4), col=c(3,1,4), lwd=2, bty="n",
       c(expression(hat(f)[P](t)),expression(hat(f)[B](t)),
                  expression(hat(f)[K](t))))
  par(op)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~models }
\keyword{ ~nonparametric }
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
