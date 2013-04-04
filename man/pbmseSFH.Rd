\name{pbmseSFH}
\alias{pbmseSFH}
\title{Parametric bootstrap mean squared error estimators of the spatial EBLUPs under a spatial Fay-Herriot model.}
\description{
Calculates the parametric bootstrap mean squared error estimates of the spatial EBLUPs obtained by fitting the spatial Fay-Herriot model, in which area effects follow a Simultaneously Autorregressive (SAR) process.
%% ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pbmseSFH(formula, vardir, proxmat, B = 100, method = "REML", MAXITER = 100,
         PRECISION = 0.0001, data)
}
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The variables included in \code{formula}
                  must have a length equal to the number of domains \code{D}. Details of model specification are given under Details.}                  
  \item{vardir}{vector containing the \code{D} sampling variances of direct estimators for each domain. 
                The values must be sorted as the variables in \code{formula}.}
  \item{proxmat}{\code{D*D} proximity matrix or data frame with values in the interval \code{[0,1]} containing the proximities between the row and column domains. The rows add up to 1. The rows and columns of this matrix must be sorted as the variables in \code{formula}.}
  \item{B}{number of bootstrap replicates. Default value is \code{100}.}
   \item{method}{type of fitting method, to be chosen between \code{"REML"} or \code{"ML"}. Default value is \code{REML}.}
  \item{MAXITER}{maximum number of iterations allowed for the Fisher-scoring algorithm. Default value is \code{100}.}
  \item{PRECISION}{convergence tolerance limit for the Fisher-scoring algorithm. Default value is \code{0.0001}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{vardir}. By default the variables are taken from the environment from which \code{pbmseSFH} is called.}  
}
\details{
This function uses random number generation. To fix the seed, use \code{set.seed}.

A typical model has the form response ~ terms where response is the (numeric) response vector and 
terms is a series of terms which specifies a linear predictor for response. A terms specification of 
the form first + second indicates all the terms in first together with all the terms in second with 
duplicates removed. A terms specification of the form first + second indicates all the terms in first 
together with all the terms in second with any duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} 
for more details of allowed formulae. }

\value{
  The function returns a list with the following objects:
  \item{est}{a list with the results of the estimation process: \code{eblup} and \code{fit}. For the description of these objects, see Value of  \code{\link{eblupSFH}} function.}
  \item{mse}{data frame containing the naive parametric bootstrap mean squared error estimates (\code{mse}) and the bias-corrected parametric bootstrap mean squared error estimates of the spatial EBLUPs (\code{msebc}).}
  
  In case that \code{formula}, \code{vardir} or \code{proxmat} contain NA values a message is printed and no action is done.
}
\references{ 

- Small Area Methods for Poverty and Living Conditions Estimates (SAMPLE), 
funded by European Commission, Collaborative Project 217565, Call identifier FP7-SSH-2007-1.

- Molina, I., Salvati, N. and Pratesi, M. (2009). Bootstrap for estimating the MSE of the Spatial EBLUP. Computational Statistics 24, 441-458.

}

\author{Isabel Molina, Monica Pratesi and Nicola Salvati.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{eblupSFH}}, \code{\link{npbmseSFH}}, \code{\link{mseSFH}} }
\examples{
data(grapes)       # Load data set
data(grapesprox)   # Load proximity matrix

# Obtain the fitting values, naive and bias-corrected parametric bootstrap MSE estimates
# using REML method 
set.seed(123)
result <- pbmseSFH(grapehect ~ area + workdays - 1, var, grapesprox, B=2, data=grapes)
result
}
\keyword{ method }

