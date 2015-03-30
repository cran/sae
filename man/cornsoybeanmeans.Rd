\name{cornsoybeanmeans}
\alias{cornsoybeanmeans}
\docType{data}
\title{Corn and soy beans mean number of pixels per segment for 12 counties in Iowa.}
\description{ County means of number of pixels per segment of corn and soy beans, from satellite data, for 12 counties in Iowa. Population size, sample size and means of auxiliary variables in data set \code{\link{cornsoybean}}.}
\usage{data(cornsoybeanmeans)}
\format{
  A data frame with 12 observations on the following 6 variables.
  \describe{
    \item{\code{CountyIndex}:}{ numeric county code.}
    \item{\code{CountyName}:}{ name of the county.}
    \item{\code{SampSegments}:}{ number of sample segments in the county (sample size).}
    \item{\code{PopnSegments}:}{ number of population segments in the county (population size).}
    \item{\code{MeanCornPixPerSeg}:}{ mean number of corn pixels per segment in the county.}
    \item{\code{MeanSoyBeansPixPerSeg}:}{ mean number of soy beans pixels per segment in the county.}
  }
}

%%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%%}
\source{ 
- Battesse, G.E., Harter, R.M. and Fuller, W.A. (1988). An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data. Journal of the American Statistical Association 83, 28-36.
}
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}

\keyword{datasets}
