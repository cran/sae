\name{pbmseSTFH}
\alias{pbmseSTFH}
\title{Parametric bootstrap mean squared error estimator of a spatio-temporal Fay-Herriot model.}
\description{
Calculates parametric bootstrap mean squared error estimates of the EBLUPs based on a spatio-temporal Fay-Herriot model with area effects following a SAR(1) process and with either uncorrelated or correlated time effects within each domain following an AR(1) process.
}
\usage{
pbmseSTFH(formula, D, T, vardir, proxmat, B = 100, model = "ST", 
          MAXITER = 100, PRECISION = 0.0001, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The variables included in \code{formula}
                  must have a length equal to \code{D*T} and sorted in the ascending order by the time instant within each domain. 
                  Details of model specification are given under Details.}                  
  \item{D}{total number of domains.}
  \item{T}{total number of time instants (constant for each domain).}
  \item{vardir}{vector containing the \code{n=D*T} sampling variances for each domain and time instant. The values must be sorted as the variables in \code{formula}.}
  \item{proxmat}{\code{D*D} proximity matrix or data frame with values in the interval \code{[0,1]} containing the proximities between the row and                       column domains. The rows add up to 1. The rows and columns of this matrix must be sorted by domain as the variables in 
                  \code{formula}.}
  \item{B}{number of bootstrap replicates. Default value is \code{100}.}
  \item{model}{type of model to be chosen between \code{"ST"} (correlated time-effects within domains) or \code{"S"} (uncorrelated time-effects 
               within domains).}
  \item{MAXITER}{maximum number of iterations allowed for the Fisher-scoring algorithm. Default value is \code{100}.}
  \item{PRECISION}{convergence tolerance limit for the Fisher-scoring algorithm. Default value is \code{0.0001}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{vardir}. By default the variables are taken from the environment from which \code{pbmseSTFH} is called.}  
}

\details{
This function uses random number generation. To fix the seed, use \code{set.seed}.

A typical model has the form response ~ terms where response is the (numeric) response vector and 
terms is a series of terms which specifies a linear predictor for response. A terms specification of 
the form first + second indicates all the terms in first together with all the terms in second with 
duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} 
for more details of allowed formulae.  
}
\value{ 
  The function returns a list with the following objects:
  \item{est}{a list with the results of the estimation process: \code{eblup} and \code{fit}. For the description of these objects, see Value of  
             \code{\link{eblupSTFH}} function.}
  \item{mse}{a vector of length \code{D*T} containing the parametric bootstrap mean squared error estimates for the \code{D} domains and \code{T} time instants.}

In case that \code{formula}, \code{vardir} or \code{proxmat} contain NA values a message is printed and no action is done.
}
\references{ 

- Small Area Methods for Poverty and Living Conditions Estimates (SAMPLE), 
funded by European Commission, Collaborative Project 217565, Call identifier FP7-SSH-2007-1. 

- Marhuenda, Y., Molina, I. and Morales, D. (2013). Small area estimation with spatio-temporal
Fay-Herriot models. Computational Statistics and Data Analysis 58, 308-325. 

}

\author{Yolanda Marhuenda, Isabel Molina and Domingo Morales.}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{eblupSTFH}} }
\examples{
data(spacetime)      # Load data set 
data(spacetimeprox)  # Load proximity matrix 

D <- nrow(spacetimeprox)            # number of domains
T <- length(unique(spacetime$Time)) # number of time instant

# Calculate MSEs of EBLUPs under the spatio-temporal Fay-Herriot model
# with uncorrelated time effects nested within domains (model S)
set.seed(123)
resultS <- pbmseSTFH(Y ~ X1 + X2, D, T, Var, spacetimeprox, B=40, 
                     model="S", data=spacetime)

# Print direct estimates, variance, "S" model estimates, mse and 
# residuals of the last time instant.
output <- data.frame(Domain=spacetime$Area, Period=spacetime$Time, 
                     Direct=spacetime$Y, EBLUP_S=resultS$est$eblup, 
                     VarDirect=spacetime$Var, MSE_S=resultS$mse, 
                     Residuals=spacetime$Y-resultS$est$eblup) 
periods <- unique(spacetime$Time)                  
lastperiod <- periods[length(periods)]   
print(output[output[,"Period"]==lastperiod,], row.names=FALSE)


# Calculate MSEs of the EBLUPs based on the spatio-temporal Fay-Herriot model
# with AR(1) time effects nested within each area
attach(spacetime)
set.seed(123)
resultST <- pbmseSTFH(Y ~ X1 + X2, D, T, vardir=Var, spacetimeprox, B=40)  

# Print direct estimates, variance, "ST" model estimates, mse and 
# residuals of the last time instant.
output <- data.frame(Domain=Area, Period=Time, Direct=Y, 
                     EBLUP_ST=resultST$est$eblup, VarDirect=Var, 
                     MSE_ST=resultST$mse, 
                     Residuals=Y-resultST$est$eblup) 
periods <- unique(Time)                  
lastperiod <- periods[length(periods)]                     
print(output[output[,"Period"]==lastperiod,], row.names=FALSE)

detach(spacetime)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }

