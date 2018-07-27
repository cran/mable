\name{chicken.embryo}
\alias{chicken.embryo}
\docType{data}
\title{Chicken Embryo Data
%%   ~~ data name/kind ... ~~
}
\description{The chicken embryo dataset which contains number of days and the corresponding frequencies.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(chicken.embryo)}
\format{
  The format is:
List of 2
 $ day: int [1:21] 1 2 3 4 5 6 7 8 9 10 ...
 $ nT : int [1:21] 6 5 11 2 2 3 0 0 0 0 ...
}

\source{
Jassim, E. W., Grossman, M., Koops, W. J. And Luykx, R. A. J. (1996). Multi-phasic analysis of embryonic
mortality in chickens. \emph{Poultry Sci.} 75, 464-71.
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Kuurman, W. W., Bailey, B. A., Koops, W. J. And Grossman, M. (2003). A model for failure of a chicken
embryo to survive incubation. \emph{Poultry Sci.} 82, 214-22.

Guan, Z. (2017) Bernstein polynomial model for grouped continuous data.
              \emph{Journal of Nonparametric Statistics}, 29(4):831-848.
}
\examples{
data(chicken.embryo)
%## maybe str(chicken.embryo) ; plot(chicken.embryo) ...
}
\keyword{datasets}
