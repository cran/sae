\name{pbmseebBHF}
\alias{pbmseebBHF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric bootstrap mean squared error estimators of EB estimators.
}
\description{This function obtains estimators of the mean squared errors of the EB estimators of domain parameters by a parametric bootstrap method. Population values of auxiliary variables are required.
}
\usage{
pbmseebBHF(formula, dom, selectdom, Xnonsample, B = 100, MC = 100, data, 
           transform = "BoxCox", lambda = 0, constant = 0, indicator)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The details of model 
                  specification are given under Details.}
  \item{dom}{\code{n*1} vector or factor (same size as \code{y} in \code{formula}) with domain codes.}
  \item{selectdom}{\code{I*1} optional vector or factor with the domain codes for which we 
                   want to estimate the indicators. It must be a subset of the domain codes in \code{dom}. 
                   If this parameter is not included, the unique domain codes included in \code{dom} 
                   are considered.}                   
  \item{Xnonsample}{matrix or data frame containing in the first column the domain codes and in the rest of columns the values of each of 
                     \code{p} auxiliary variables for the out-of-sample units in each selected domain.}
  \item{B}{number of bootstrap replicates. Default value is \code{100}.}  
  \item{MC}{number of Monte Carlo replicates for the empirical approximation of the EB estimator. Default value is \code{100}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{dom}. By default the variables are taken from the 
              environment from which \code{pbmseebBHF} is called.} 
   \item{transform}{type of transformation for the dependent variable to be chosen between the \code{"BoxCox"} and \code{"power"} families so that the dependent variable in \code{formula} follows approximately a Normal distribution. Default value is \code{"BoxCox"}.}
   \item{lambda}{value for the parameter of the family of transformations specified in \code{transform}. Default value is \code{0}, which gives the log transformation for the two possible families.}
   \item{constant}{constant added to the dependent variable before doing the transformation, to achieve a distribution close to Normal. Default value is \code{0}.}
    \item{indicator}{function of the (untransformed) variable on the left hand side of \code{formula} that we want to estimate in each domain.}
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
  \item{est}{a list with the results of the estimation process: \code{eb} and \code{fit}. For the description of these objects, see Value of  \code{\link{ebBHF}} function.}
  \item{mse}{data frame with number of rows equal to number of selected domains, containing in its columns the domain codes (\code{domain}) and the parametric bootstrap mean squared error estimates of \code{indicator} (\code{mse}).}

Cases with NA values in \code{formula} or \code{dom} are ignored. 
}
\references{ 

- Small Area Methods for Poverty and Living Conditions Estimates (SAMPLE), 
funded by European Commission, Collaborative Project 217565, Call identifier FP7-SSH-2007-1.

- Molina, I. and Rao, J.N.K. (2010). Small Area Estimation of Poverty Indicators. The Canadian Journal of Statistics 38, 369-385.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ebBHF}} }
\examples{
data(incomedata)         # Load data set
attach(incomedata)

# Construct design matrix for sample elements
Xs<-cbind(age2,age3,age4,age5,nat1,educ1,educ3,labor1,labor2)

# Select the domains to compute EB estimators
data(Xoutsamp)
domains <- c(5)

# Poverty incidence indicator
povertyline <- 0.6*median(incomedata$income)
povertyline                         # 6477.484
povinc <- function(y)    
{
   z <- 6477.484
   result <- mean(y<z)
   return (result)
}

# Compute parametric bootstrap MSE estimators of the EB 
# predictors of poverty incidence. Take constant=3600 to achieve 
# approximately symmetric residuals.
set.seed(123)
result <- pbmseebBHF(income~Xs, dom=prov, selectdom=domains,
                     Xnonsample=Xoutsamp, B=2, MC=2, constant=3600,
                     indicator=povinc)
result$est$eb
result$mse
result$est$fit$refvar

detach(incomedata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
