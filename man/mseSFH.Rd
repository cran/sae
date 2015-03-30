\name{mseSFH}
\alias{mseSFH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mean squared error estimator of the spatial EBLUP under a spatial Fay-Herriot model.}
\description{
Calculates analytical mean squared error estimates of the spatial EBLUPs obtained from the fit of a spatial Fay-Herriot model, in which area effects follow a Simultaneously Autorregressive (SAR) process.
}
\usage{
mseSFH(formula, vardir, proxmat, method = "REML", MAXITER = 100, 
       PRECISION = 0.0001, data)
}

\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The variables included in \code{formula}
                  must have a length equal to the number of domains \code{D}. Details of model specification are given under Details.}
  \item{vardir}{vector containing the \code{D} sampling variances of direct estimators for each domain. 
                The values must be sorted as the variables in \code{formula}.}
  \item{proxmat}{\code{D*D} proximity matrix or data frame with values in the interval \code{[0,1]} containing the proximities between the row and column domains. The rows add up to 1. The rows and columns of this matrix must be sorted as the variables in \code{formula}.}
  \item{method}{ type of fitting method, to be chosen between \code{"REML"} or \code{"ML"}. Default value is \code{REML}.}
  \item{MAXITER}{maximum number of iterations allowed for the Fisher-scoring algorithm. Default value is \code{100}.}
  \item{PRECISION}{convergence tolerance limit for the Fisher-scoring algorithm. Default value is \code{0.0001}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{vardir}. By default the variables are taken from the environment from which \code{mseSFH} is called.}   
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  The function returns a list with the following objects:
  \item{est}{a list with the results of the estimation process: \code{eblup} and \code{fit}. For the description of these objects, see Value of  \code{\link{eblupSFH}} function.}
  \item{mse }{a vector with the analytical mean squared error estimates of the spatial EBLUPs.}
  
In case that \code{formula}, \code{vardir} or \code{proxmat} contain NA values a message is printed and no action is done.
}
\references{ 

- Small Area Methods for Poverty and Living Conditions Estimates (SAMPLE), 
funded by European Commission, Collaborative Project 217565, Call identifier FP7-SSH-2007-1.

- Molina, I., Salvati, N. and Pratesi, M. (2009). Bootstrap for estimating the MSE of
the Spatial EBLUP. Computational Statistics 24, 441-458.

- Singh, B., Shukla, G. and Kundu, D. (2005). Spatio-temporal models in small area estimation. Survey Methodology 31, 183-195.

}
\author{Isabel Molina, Monica Pratesi and Nicola Salvati.}

%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{eblupSFH}}, \code{\link{npbmseSFH}}, \code{\link{pbmseSFH}} }
\examples{
data(grapes)       # Load data set
data(grapesprox)   # Load proximity matrix 

# Calculate analytical MSE estimates using REML method
result <- mseSFH(grapehect ~ area + workdays - 1, var, grapesprox, data=grapes)
result
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
