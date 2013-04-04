\name{milk}
\alias{milk}
\docType{data}
\title{ Data on fresh milk expenditure.}
\description{
  Data on fresh milk expenditure, used by Arora and Lahiri (1997) and by You and Chapman (2006). 
}
\usage{data(milk)}
\format{
  A data frame with 43 observations on the following 6 variables.
  \describe{
    \item{\code{SmallArea}:}{areas of inferential interest.}
    \item{\code{ni}:}{sample sizes of small areas.}
    \item{\code{yi}:}{average expenditure on fresh milk for the year 1989 (direct estimates for the small areas).}
    \item{\code{SD}:}{estimated standard deviations of \code{yi}.}
    \item{\code{CV}:}{estimated coefficients of variation of \code{yi}.}
    \item{\code{MajorArea}:}{major areas created by You and Chapman (2006). These areas have similar direct estimates and produce a large CV reduction when using a FH model.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{

- Arora, V. and Lahiri, P. (1997). On the superiority of the Bayesian method over the BLUP in small area estimation problems. Statistica Sinica 7, 1053-1063.

- You, Y. and Chapman, B. (2006). Small area estimation using area level models and estimated sampling variances. Survey Methodology 32, 97-103.

}
\keyword{datasets}
