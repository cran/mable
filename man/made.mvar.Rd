% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/made.r
\name{made.mvar}
\alias{made.mvar}
\title{Minimum Approximate Distance Estimate
 of Multivariate Density Function}
\usage{
made.mvar(
  x,
  M0 = 1L,
  M,
  search = TRUE,
  interval = NULL,
  mar.deg = TRUE,
  method = c("cd", "quadprog"),
  controls = mable.ctrl(),
  progress = TRUE
)
}
\arguments{
\item{x}{an \code{n x d} matrix or \code{data.frame} of multivariate sample of size \code{n}}

\item{M0}{a positive integer or a vector of \code{d} positive integers specify
starting candidate degrees for searching optimal degrees.}

\item{M}{a positive integer or a vector of \code{d} positive integers specify
the maximum candidate or the given model degrees for the joint density.}

\item{search}{logical, whether to search optimal degrees between \code{M0} and \code{M} 
or not but use \code{M} as the given model degrees for the joint density.}

\item{interval}{a vector of two endpoints or a \code{2 x d} matrix, each column containing 
the endpoints of support/truncation interval for each marginal density.
If missing, the i-th column is assigned as \code{c(min(x[,i]), max(x[,i]))}.}

\item{mar.deg}{logical, if TRUE, the optimal degrees are selected based
on marginal data, otherwise, the optimal degrees are chosen the joint data. See details.}

\item{method}{method for finding minimum distance estimate. "cd": coordinate-descent;}

\item{controls}{Object of class \code{mable.ctrl()} specifying iteration limit
and the convergence criterion \code{eps}. Default is \code{\link{mable.ctrl}}. See Details.}

\item{progress}{if TRUE a text progressbar is displayed}
}
\value{
A list with components
\itemize{
 \item \code{m} a vector of the selected optimal degrees by the method of change-point
 \item \code{p} a vector of the mixture proportions \eqn{p(j_1, \ldots, j_d)}, arranged in the 
  column-major order of \eqn{j = (j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}.
 \item \code{minD}  the minimum distance at an optimal degree \code{m}
 \item \code{pval}  the p-values of change-points for choosing the optimal degrees for the 
   marginal densities
 \item \code{M} the vector \code{(m1, m2, ... , md)}, where \code{mi} is the largest candidate 
   degree when the search stoped for the \code{i}-th marginal density
 \item \code{interval} support hyperrectangle \eqn{[a, b]=[a_1, b_1] \times \cdots \times [a_d, b_d]}
 \item \code{convergence} An integer code. 0 indicates successful completion(the EM iteration is   
   convergent). 1 indicates that the iteration limit \code{maxit} had been reached in the EM iteration;
}
}
\description{
Minimum Approximate Distance Estimate
 of Multivariate Density Function
}
\details{
A \eqn{d}-variate density \eqn{f} on a hyperrectangle \eqn{[a, b]
  =[a_1, b_1] \times \cdots \times [a_d, b_d]} can be approximated 
  by a mixture of \eqn{d}-variate beta densities on \eqn{[a, b]}, 
  \eqn{\beta_{mj}(x) = \prod_{i=1}^d\beta_{m_i,j_i}[(x_i-a_i)/(b_i-a_i)]/(b_i-a_i)},
  with proportion \eqn{p(j_1, \ldots, j_d)}, \eqn{0 \le j_i \le m_i, i = 1, \ldots, d}. 
  If \code{search=TRUE} then the model degrees are chosen using a method of change-point based on 
  the marginal data if \code{mar.deg=TRUE} or the joint data if \code{mar.deg=FALSE}. 
  If \code{search=FALSE}, then the model degree is specified by \eqn{M}.
  For large data and multimodal density, the search for the model degrees is 
  very time-consuming. In this case, it is suggested that use \code{method="cd"} 
  and select the degrees based on marginal data using \code{\link{mable}} or 
  \code{\link{optimable}}.
}
\examples{
## Old Faithful Data
\donttest{  
 library(mable)
 a<-c(0, 40); b<-c(7, 110)
 ans<- made.mvar(faithful, M = c(46,19), search =FALSE, method="quadprog", 
 interval = rbind(a,b), progress=FALSE)
 plot(ans, which="density") 
 plot(ans, which="cumulative")
}
}
\references{
Guan, Z. (2016) Efficient and robust density estimation using Bernstein type polynomials. 
   \emph{Journal of Nonparametric Statistics}, 28(2):250-271.

Wang, T. and Guan, Z.,(2019) Bernstein Polynomial Model for Nonparametric Multivariate Density,    
   \emph{Statistics}, Vol. 53, no. 2, 321-338
}
\seealso{
\code{\link{mable}}, \code{\link{optimable}}
}
\author{
Zhong Guan <zguan@iu.edu>
}
\concept{density estimation}
\concept{multivariate Bernstein polynomial model}
\keyword{distribution}
\keyword{multivariate}
\keyword{nonparametric}
