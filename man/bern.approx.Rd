\name{bern.approx}
\alias{bern.approx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bernstein Polynomial Approximation}
\description{Returns Bernstein polynomial approximation \eqn{f_m(x; p)} or \eqn{F_m(x; p)}}
\usage{
bern.approx(x, p, a=0, b=1, density = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector of values in [0,1] at which the Bernstein polynomial approximation \eqn{f_m} of \eqn{F_m} are evaluated.}
  \item{p}{the coefficients of the Bernstein polynomial approximation.}
  \item{a}{left endpoint of supporting interval}
  \item{b}{right endpoint of supporting interval}
  \item{density}{logical; whether a density or a distribution function to be calculated. Default is TRUE.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The \eqn{m+1} components of \eqn{p} must be nonnegative.
}
\value{A vector of \eqn{f_m(x; p)} or \eqn{F_m(x; p)} values at \eqn{x}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

\seealso{\code{\link{mable.em}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\donttest{
# classical Bernstein polynomial
a<--4; b<-4; m<-200
x<-seq(a,b,len=512)
u<-(0:m)/m
p<-dnorm(a+(b-a)*u)
plot(x, dnorm(x), type="l")
lines(x, (b-a)*bern.approx(x, p, a, b)/(m+1), lty=2, col=2)
legend(a, dnorm(0), lty=1:2, col=1:2, c(expression(f(x)==phi(x)),
                expression(B^{f}*(x))))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~distribution }
\keyword{ ~models }
\keyword{ ~nonparametric }
\keyword{ ~smooth }% __ONLY ONE__ keyword per line
