\name{cornsoybean}
\alias{cornsoybean}
\docType{data}
\title{ Corn and soy beans survey and satellite data in 12 counties in Iowa.}
\description{ Survey and satellite data for corn and soy beans in 12 Iowa counties, obtained from the 1978 June Enumerative Survey of the U.S. Department of Agriculture and from land observatory satellites (LANDSAT) during the 1978 growing season.
 }
\usage{data(cornsoybean)}
\format{
  A data frame with 37 observations on the following 5 variables.
  \describe{
    \item{\code{County}:}{ numeric county code.}
    \item{\code{CornHec}:}{ reported hectares of corn from the survey.}
    \item{\code{SoyBeansHec}:}{ reported hectares of soy beans from the survey.}
    \item{\code{CornPix}:}{ number of pixels of corn in sample segment within county, from satellite data. }
    \item{\code{SoyBeansPix}:}{ number of pixels of soy beans in sample segment within county, from satellite data.}
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
