\name{mseFH}
\alias{mseFH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mean squared error estimator of the EBLUP under a Fay-Herriot model.}
\description{
 Calculates the mean squared error estimator of the EBLUP under a Fay-Herriot model.
 The EBLUP might have been obtained by either ML, REML or by FH 
 fitting methods.
}
\usage{
mseFH(formula, vardir, method = "REML", MAXITER = 100, PRECISION = 0.0001, B = 0, 
      data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The variables included in \code{formula}
                  must have a length equal to the number of domains \code{D}. Details of model 
                  specification are given under Details.}                  
  \item{vardir}{vector containing the \code{D} sampling variances of direct estimators for each domain. The values must be sorted as the 
                variables in \code{formula}.}  
  \item{method}{method used to fit the Fay-Herriot model, which can be either \code{"ML"}, \code{"REML"} or \code{"FH"} methods.
                 Default is \code{"REML"} method.}
 \item{MAXITER}{maximum number of iterations allowed in the Fisher-scoring algorithm.
                  Default is 100 iterations.}
  \item{PRECISION}{convergence tolerance limit for the Fisher-scoring algorithm. Default value is \code{0.0001}.}   
\item{B}{number of bootstrap replicates to calculate the goodness-of-fit measures proposed by Marhuenda et al. (2014). Default value is \code{0} indicating that these measures are not calculated.}  
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{vardir}. By default the variables are taken from the environment from which \code{mseFH} is called.}                 
}
\details{

A typical model has the form response ~ terms where response is the (numeric) response vector
and terms is a series of terms which specifies a linear predictor for response. A terms specification
of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed.

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} for more details of allowed formulae

}
\value{
  The function returns a list with the following objects:
  \item{est}{a list with the results of the estimation process: \code{eblup} and \code{fit}. For the description of these objects, see Value of  \code{\link{eblupFH}} function.}
  \item{mse}{a vector with the estimated mean squared errors of the EBLUPs for the small domains.}
  
  In case that \code{formula} or \code{vardir} contain NA values a message is printed and no action is done.
}
\references{

- Datta, G.S. and Lahiri, P. (2000). A unified measure of uncertainty of estimated
best linear unbiased predictors in small area estimation problems. Statistica Sinica 10, 613-627.

- Datta, G.S., Rao, J.N.K. and Smith D.D. (2005). On measuring the variability of small area estimators under a basic area level model. Biometrika 92, 183-196.

- Fay, R.E. and Herriot, R.A. (1979). Estimation of income from small places: An application of James-Stein procedures to census data. Journal of the American Statistical Association 74, 269-277.

- Jiang, J. (1996). REML estimation: asymptotic behavior and related topics.
Annals of Statistics 24, 255-286.

- Marhuenda, Y., Morales, D. and Pardo, M.C. (2014). Information criteria for Fay-Herriot model selection. Computational Statistics and Data Analysis 70, 268-280.

- Prasad, N. and Rao, J. (1990). The estimation of the mean squared error of small-area estimators. Journal of the American Statistical Association 85, 163-171.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{eblupFH}} }
\examples{
# Load data set  
data(milk)   
attach(milk)

# Fit Fay-Herriot model using ML method with indicators 
# of 4 Major Areas as explanatory variables and compute
# estimated MSEs of EB estimators
resultML <- mseFH(yi ~ as.factor(MajorArea), SD^2, method="ML")
resultML

# Fit Fay-Herriot model using REML method and compute 
# estimated MSEs of EB estimators
resultREML <- mseFH(yi ~ as.factor(MajorArea), SD^2)
resultREML

# Fit Fay-Herriot model using FH method and compute 
# estimated MSEs of EB estimators
resultFH <- mseFH(yi ~ as.factor(MajorArea), SD^2, method="FH")
resultFH

detach(milk)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
