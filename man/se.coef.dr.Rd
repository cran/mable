% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mable-dr.r
\name{se.coef.dr}
\alias{se.coef.dr}
\title{Standard errors of coefficients in density ratio model}
\usage{
se.coef.dr(
  obj,
  grouped = FALSE,
  B = 500L,
  parallel = FALSE,
  ncore = NULL,
  controls = mable.ctrl()
)
}
\arguments{
\item{obj}{Class 'mable_dr' object return by \code{mable.dr} or \code{mable.dr.group}  functions}

\item{grouped}{logical: are data grouped or not.}

\item{B}{number of bootstrap runs.}

\item{parallel}{logical: do parallel or not.}

\item{ncore}{number of cores used for parallel computing. Default is half of availables.}

\item{controls}{Object of class \code{mable.ctrl()} specifying iteration limit
and the convergence criterion for EM and Newton iterations. Default is 
\code{\link{mable.ctrl}}. See Details.}
}
\value{
the estimated standard errors
}
\description{
Bootstrap estimates of standard errors for the regression
  coefficients which are estimated by maximum approximate Bernstein/Beta 
  likelihood estimation method in a density ratio model based on two-sample 
  raw data.
}
\details{
Bootstrap method is used based on bootstrap samples generated from
the MABLE's of the densities f0 and f1. The bootstrap samples are fitted by
the Bernstein polynomial model and the \code{glm()} to obtain bootstrap 
versions of coefficient estimates.
}
