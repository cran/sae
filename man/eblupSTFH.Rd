\name{eblupSTFH}
\alias{eblupSTFH}
\title{EBLUPs based on a spatio-temporal Fay-Herriot model.}
\description{
Fits a spatio-temporal Fay-Herriot model with area effects following a SAR(1) process
and with either uncorrelated or AR(1) time effects.
}
\usage{
eblupSTFH(formula, D, T, vardir, proxmat, model = "ST", MAXITER = 100, 
          PRECISION = 0.0001, theta_iter = FALSE,
          sigma21_start = 0.5 * median(vardir), rho1_start = 0.5, 
          sigma22_start = 0.5 * median(vardir), rho2_start = 0.5, 
          data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The variables included in \code{formula}
                  must have a length equal to \code{D*T} and sorted in the ascending order by the time instant within each domain. 
                  Details of model specification are given under Details.}
  \item{D}{total number of domains.}
  \item{T}{total number of time instants (constant for all domains).}
  \item{vardir}{vector containing the \code{D*T} sampling variances of direct estimators for each domain and time instant. 
                  The values must be sorted as the variables in \code{formula}.}
  \item{proxmat}{\code{D*D} proximity matrix or data frame with values in the interval \code{[0,1]} containing the proximities between the row and                       column domains. The rows add up to 1. The rows and columns of this matrix must be sorted by domain as the variables in  
                 \code{formula}.}
  \item{model}{type of model to be chosen between \code{"ST"} (AR(1) time-effects within each domain) or \code{"S"} (uncorrelated time effects within each domain). Default model is "ST".  }
  \item{MAXITER}{maximum number of iterations allowed for the Fisher-scoring algorithm. Default value is \code{100}.}
  \item{PRECISION}{convergence tolerance limit for the Fisher-scoring algorithm. Default value is \code{0.0001}.}
  \item{theta_iter}{If \code{TRUE} the estimated values of area effects variance, area effects spatial autocorrelation, area-time effects variance and time autocorrelation parameter of the area-time effects of each iteration of the fitting algorithm are returned in \code{estvarcomp_iterations} variable.}
  \item{sigma21_start}{Starting value of the area effects variance in the fitting algorithm. Default value is \code{0.5*median(vardir)}.}
  \item{rho1_start}{Starting value of the area effects spatial autocorrelation parameter in the fitting algorithm. Default value is \code{0.5}.}
  \item{sigma22_start}{Starting value of the area-time effects variance in the fitting algorithm. Default value is \code{0.5*median(vardir)}.}
  \item{rho2_start}{Starting value of the time autocorrelation parameter of the area-time effects in the fitting algorithm. Default value is \code{0.5}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{vardir}. By default the variables are taken from the environment from which \code{eblupSTFH} is called.}  
}
\details{
A typical model has the form response ~ terms where response is the (numeric) response vector and 
terms is a series of terms which specifies a linear predictor for response. A terms specification of 
the form first + second indicates all the terms in first together with all the terms in second with 
duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} 
for more details of allowed formulae.  
}

\value{
The function returns a list with the following objects:
 \item{eblup}{a column vector with length \code{D*T} with the values of the estimators 
               for the \code{D} domains and \code{T} time instants.}
 \item{fit}{a list containing the following objects:
  \itemize{               
    \item{\code{model}: }{type of model \code{"S"} or \code{"ST"}.}
    \item{\code{convergence}: }{a logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.}
    \item{\code{iterations}: }{number of iterations performed by the Fisher-scoring algorithm.} 
    \item{\code{estcoef}: }{a data frame with the estimated model coefficients in the first column (\code{beta}), 
              their asymptotic standard errors in the second column (\code{std.error}), 
              the \eqn{t} statistics in the third column \code{(tvalue)} and the p-values of the significance 
              of each coefficient in last column \code{(pvalue)}.}     
    \item{\code{estvarcomp}: }{a data frame with the estimated values of the variances and correlation coefficients in the first column (\code{estimate}) and their asymptotic standard errors in the second column (\code{std.error}).}
    \item{\code{estvarcomp_iterations}: }{if \code{theta_iter=TRUE}, this component contains a data frame with the estimated values of the variances and correlation coefficients obtained for each iteration of the fitting algorithm.}    
     \item{\code{goodness}: }{vector containing three goodness-of-fit measures: loglikehood, AIC and BIC.}
    }
  }   
  In case that \code{formula}, \code{vardir} or \code{proxmat} contain NA values a message is printed and no action is done.
}
\references{  

- Small Area Methods for Poverty and Living Conditions Estimates (SAMPLE), 
funded by European Commission, Collaborative Project 217565, Call identifier FP7-SSH-2007-1.

- Marhuenda, Y., Molina, I. and Morales, D. (2013). Small area estimation with spatio-temporal
Fay-Herriot models. Computational Statistics and Data Analysis 58, 308-325. 

}

\author{Yolanda Marhuenda, Isabel Molina and Domingo Morales.}
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{pbmseSTFH}}}

\examples{
data(spacetime)      # Load data set 
data(spacetimeprox)  # Load proximity matrix 

D <- nrow(spacetimeprox)            # number of domains
T <- length(unique(spacetime$Time)) # number of time instant

# Fit model S with uncorrelated time effects for each domain
resultS <- eblupSTFH(Y ~ X1 + X2, D, T, Var, spacetimeprox, "S",
                     theta_iter=TRUE, data=spacetime)
rowsT <- seq(T, T*D, by=T)
data.frame(Domain=spacetime$Area[rowsT], EBLUP_S=resultS$eblup[rowsT])
resultS$fit

# Fit model ST with AR(1) time effects for each domain
resultST <- eblupSTFH(Y ~ X1 + X2, D, T, Var, spacetimeprox,
                      theta_iter=TRUE, data=spacetime)
data.frame(Domain=spacetime$Area[rowsT], EBLUP_ST=resultS$eblup[rowsT])
resultST$fit
}
\keyword{ method }

