\name{pbmseBHF}
\alias{pbmseBHF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric bootstrap mean squared error estimators of the EBLUPs of means obtained under a nested error linear regression model.}
\description{Calculates, for selected domains, parametric bootstrap mean squared error  
 estimators of the EBLUPs of means, when EBLUPs are obtained from a nested error linear regression model.}

\usage{
pbmseBHF(formula, dom, selectdom, meanxpop, popnsize, B = 200, method = "REML", 
         data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The details of model 
                  specification are given under Details.}
  \item{dom}{\code{n*1} vector or factor (same size as \code{y} in \code{formula}) with domain codes.}
  \item{selectdom}{\code{I*1} optional vector or factor with the domain codes for which we 
                   want to estimate the means. It must be a subset of the domain codes in \code{dom}. 
                   If this parameter is not included all the domain codes included in \code{dom} 
                   are considered.}  
  \item{meanxpop}{\code{D*(p+1)} data frame with domain codes in the first column. Each remaining column contains the  
                  population means of each of the \code{p} auxiliary variables for the D domains. 
                  The domains considered in \code{meanxpop} must contain those specified in \code{selectdom} (D>=I).}
  \item{popnsize}{\code{D*2} data frame with domain codes in the first column
                  and the corresponding domain population sizes in the second column. 
                  The domains considered in \code{popnsize} must contain those specified in \code{selectdom} (D>=I).}
  \item{B}{number of bootstrap replicates. Default is \code{50}.}
  \item{method}{a character string. If \code{"REML"} the model is fitted by maximizing the restricted 
                log-likelihood. If \code{"ML"} the log-likelihood is maximized. Defaults to \code{"REML"}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{dom}. By default the variables are taken from the environment from which \code{pbmseBHF} is called.}                
                
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
  \item{est}{a list with the results of the estimation process: \code{eblup} and \code{fit}. For the description of these objects, see Value of  \code{\link{eblupBHF}} function.}
  \item{mse}{data frame with number of rows equal to number of selected domains, containing in its columns the domain codes (\code{domain}) and the parametric bootstrap mean squared error estimators (\code{mse}).}

Cases with NA values in \code{formula} or \code{dom} are ignored. 
}
\references{

- Gonzalez-Manteiga, W., Lombardia, M., Molina, I., Morales, D. and Santamaria, L. (2008). Analytic and bootstrap approximations of prediction errors under a multivariate Fay-Herriot model. Computational Statistics and Data Analysis 52, 5242-5252.

- Molina, I. and Rao, J.N.K. (2010). Small Area Estimation of Poverty Indicators. The Canadian Journal of Statistics 38, 369-385.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
%%\section{Warning }{MIRAR} 

\seealso{ \code{\link{eblupBHF}} }
\examples{
# Load data set for segments (units within domains)
data(cornsoybean)

# Load data set for counties
data(cornsoybeanmeans)
attach(cornsoybeanmeans)

# Construct data frame with county means of auxiliary variables for 
# domains. First column must include the county code
Xmean <- data.frame(CountyIndex, MeanCornPixPerSeg, MeanSoyBeansPixPerSeg)
Popn  <- data.frame(CountyIndex, PopnSegments)

# Compute parametric bootstrap MSEs of the EBLUPs of means of crop areas 
# for each county.
set.seed(123) 
result <- pbmseBHF(CornHec ~ CornPix + SoyBeansPix, dom=County,
                   selectdom=c(10,1,5), meanxpop=Xmean, popnsize=Popn,
                   B=50, data=cornsoybean)
result

detach(cornsoybeanmeans)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
