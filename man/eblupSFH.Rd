\name{eblupSFH}
\alias{eblupSFH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ EBLUPs based on a spatial Fay-Herriot model.}
\description{
This function gives small area estimators based on a spatial Fay-Herriot model, where area effects follow a SAR(1) process. Fitting method can be chosen between REML and ML.
}
\usage{
eblupSFH(formula, vardir, proxmat, method = "REML", MAXITER = 100, 
         PRECISION = 0.0001, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The variables included in \code{formula}
                  must have a length equal to the number of domains \code{D}. Details of model specification are given under Details.}
  \item{vardir}{vector containing the \code{D} sampling variances of direct estimators for each domain. 
                The values must be sorted as the variables in \code{formula}.}
  \item{proxmat}{\code{D*D} proximity matrix or data frame with values in the interval \code{[0,1]} containing the proximities between the row and column domains. The rows add up to 1. The rows and columns of this matrix must be sorted as the elements in \code{formula}.}
  \item{method}{type of fitting method, to be chosen between \code{"REML"} or \code{"ML"}. Default value is \code{REML}.}
  \item{MAXITER}{maximum number of iterations allowed for the Fisher-scoring algorithm. Default value is \code{100}.}
  \item{PRECISION}{convergence tolerance limit for the Fisher-scoring algorithm. Default value is \code{0.0001}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{vardir}. By default the variables are taken from the environment from which \code{eblupSHF} is called.}                 
}
\details{
A typical model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} 
for more details of allowed formulae.  
}
\value{
 The function returns a list with the following objects:
  \item{eblup}{vector with the values of the estimators for the domains.}
  \item{fit}{a list containing the following objects:
  \itemize{
    \item{\code{method}: }{type of fitting method applied (\code{"REML"} or \code{"ML"}).}
    \item{\code{convergence}: }{a logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.}
    \item{\code{iterations}: }{number of iterations performed by the Fisher-scoring algorithm.}
    \item{\code{estcoef}: }{a data frame with the estimated model coefficients in the first column (\code{beta}), 
              their asymptotic standard errors in the second column (\code{std.error}), 
              the \eqn{t} statistics in the third column \code{(tvalue)} and the p-values of the significance 
              of each coefficient in last column \code{(pvalue)}.}    
    \item{\code{refvar}: }{estimated random effects variance.}
    \item{\code{spatialcorr}: }{estimated spatial correlation parameter.}    
    \item{\code{goodness}: }{vector containing three goodness-of-fit measures: loglikehood, AIC and BIC.}
    }
  }

In case that \code{formula}, \code{vardir} or \code{proxmat} contain NA values a message is printed and no action is done.
}
\references{ 

- Small Area Methods for Poverty and Living Conditions Estimates (SAMPLE), funded by European Commission, Collaborative Project 217565, Call identifier FP7-SSH-2007-1.

- Molina, I., Salvati, N. and Pratesi, M. (2009). Bootstrap for estimating the MSE of the Spatial EBLUP. Computational Statistics 24, 441-458.

- Petrucci, A. and Salvati, N. (2006). Small area estimation for spatial correlation in watershed erosion assessment.
Journal of Agricultural, Biological and Environmental Statistics 11, 169-182.

- Pratesi, M. and Salvati, N. (2008). Small area estimation: the EBLUP estimator based on spatially correlated random area effects. 
Statistical Methods & Applications 17, 113-141.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}
\author{Isabel Molina, Monica Pratesi and Nicola Salvati.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mseSFH}}, \code{\link{npbmseSFH}}, \code{\link{pbmseSFH}} }
\examples{
data(grapes)       # Load data set
data(grapesprox)   # Load proximity matrix 

# Fit Spatial Fay-Herriot model using ML method
resultML <- eblupSFH(grapehect ~ area + workdays - 1, var, grapesprox,
                     method="ML", data=grapes)
resultML

# Fit Spatial Fay-Herriot model using REML method
resultREML <- eblupSFH(grapehect ~ area + workdays - 1, var, grapesprox,
                       data=grapes)
resultREML
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }

