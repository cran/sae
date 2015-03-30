\name{ebBHF}
\alias{ebBHF}
\title{EB estimators of an indicator with non-sample values of auxiliary variables.}
\description{
Fits by REML method the unit level model of Battese, Harter and Fuller (1988) to a transformation of the specified dependent variable by a Box-Cox family or power family and obtains Monte Carlo approximations 
of EB estimators of the specified small area indicators, when the values of auxiliary variables for out-of-sample units are available.
}
\usage{
ebBHF(formula, dom, selectdom, Xnonsample, MC = 100, data,
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
   \item{Xnonsample}{matrix or data frame containing in the first column the domain codes and in the rest of columns the values of each of \code{p} auxiliary variables for the out-of-sample units in each selected domain. The domains considered in \code{Xnonsample} must contain at least those specified in \code{selectdom}.}
   \item{MC}{number of Monte Carlo replicates for the empirical approximation of the EB estimator. Default value is \code{MC=100}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{dom}. By default the variables are taken from the environment from which \code{ebBHF} is called.}   
   \item{transform}{type of transformation for the dependent variable to be chosen between the \code{"BoxCox"} and \code{"power"} families so that the dependent variable in \code{formula} follows approximately a Normal distribution. Default value is \code{"BoxCox"}.}
   \item{lambda}{value for the parameter of the family of transformations specified in \code{transform}. Default value is \code{0}, which gives the log transformation for the two possible families.}
   \item{constant}{constant added to the dependent variable before doing the transformation, to achieve a distribution close to Normal. Default value is \code{0}.}
   \item{indicator}{function of the (untransformed) variable on the left hand side of \code{formula} that we want to estimate in each domain.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function uses random number generation. To fix the seed, use \code{set.seed}.

A typical model has the form response ~ terms where response is the (numeric) response vector and 
terms is a series of terms which specifies a linear predictor for response. A terms specification of 
the form first + second indicates all the terms in first together with all the terms in second with 
duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} 
for more details of allowed formulae. 

}
\value{ The function returns a list with the following objects:
  \item{eb}{data frame with number of rows equal to number of selected domains,
                          containing in its columns the domain codes (\code{domain}), the EB 
                          estimators of \code{indicator} (\code{eb}) and the sample sizes (\code{sampsize}). 
For domains with zero sample size, the EB estimators are based on the synthetic regression. For domains in \code{selectdom} not included in \code{Xnonsample} the EB estimators are NA.}
  \item{fit}{a list containing the following objects: 
  \itemize{
    \item{\code{summary}: }{summary of the unit level model fitting.}
    \item{\code{fixed}: }{vector with the estimated values of the fixed regression coefficient.}
    \item{\code{random}: }{vector with the predicted random effects.}
    \item{\code{errorvar}: }{estimated model error variance.}
    \item{\code{refvar}: }{estimated random effects variance.}
    \item{\code{loglike}: }{log-likelihood.}
    \item{\code{residuals}: }{vector with raw residuals from the model fit.}
    }  
  } 

Cases with NA values in \code{formula} or \code{dom} are ignored. 
}
\references{ 

- Molina, I. and Rao, J.N.K. (2010). Small Area Estimation of Poverty Indicators. The Canadian Journal of Statistics 38, 369-385.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{  \code{\link{pbmseebBHF}} }
\examples{
data(incomedata)         # Load data set
attach(incomedata)

# Construct design matrix for sample elements
Xs <- cbind(age2, age3, age4, age5, nat1, educ1, educ3, labor1, labor2)

# Select the domains to compute EB estimators. 
data(Xoutsamp)
domains <- unique(Xoutsamp[,"domain"])

# Poverty gap indicator
povertyline <- 0.6*median(income)
povertyline                         # 6477.484
povgap <- function(y)     
{
   z <- 6477.484
   result <- mean((y<z) * (z-y) / z) 
   return (result)
}

# Compute EB predictors of poverty gap. The value constant=3600 is selected
# to achieve approximately symmetric residuals.
set.seed(123)
result <- ebBHF(income ~ Xs, dom=prov, selectdom=domains,
                Xnonsample=Xoutsamp, MC=10, constant=3600, indicator=povgap)
result$eb
result$fit$summary
result$fit$fixed
result$fit$random[,1]
result$fit$errorvar
result$fit$refvar
result$fit$loglike
result$fit$residuals[1:10]

detach(incomedata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
