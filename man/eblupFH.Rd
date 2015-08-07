\name{eblupFH}
\alias{eblupFH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ EBLUPs based on a Fay-Herriot model.}
\description{ This function gives the EBLUP (or EB predictor under normality) based on a Fay-Herriot model. Fitting method can be chosen between ML, REML and FH methods.}
\usage{
eblupFH(formula, vardir, method = "REML", MAXITER = 100, PRECISION = 0.0001, 
        B = 0, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ an object of class \code{\link{formula}} (or one that can be coerced to that class): 
                  a symbolic description of the model to be fitted. The variables included in \code{formula}
                  must have a length equal to the number of domains \code{D}. Details of model 
                  specification are given under Details.}
  \item{vardir}{ vector containing the \code{D} sampling variances of direct estimators for each domain. The values must be sorted as the variables in \code{formula}.}
  \item{method}{ type of fitting method, to be chosen between \code{"ML"}, \code{"REML"} or \code{"FH"} methods.}
  \item{MAXITER}{ maximum number of iterations allowed in the Fisher-scoring algorithm.
                  Default is 100 iterations.}
  \item{PRECISION}{ convergence tolerance limit for the Fisher-scoring algorithm. Default value is \code{0.0001}.}
  \item{B}{number of bootstrap replicates to calculate the goodness-of-fit measures proposed by Marhuenda et al. (2014). Default value is \code{0} indicating that these measures are not calculated.}  
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{vardir}. By default the variables are taken from the environment from which \code{eblupFH} is called.}                  
}
\details{
A typical model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}}  for more details of allowed formulae. 
}
\value{ 
  The function returns a list with the following objects:
  \item{eblup}{vector with the values of the estimators for the domains.}
  \item{fit}{a list containing the following objects:
  \itemize{
    \item{\code{method}: }{type of fitting method applied (\code{"REML"}, \code{"ML"}or \code{"FH"}).}
    \item{\code{convergence}: }{a logical value equal to \code{TRUE} if Fisher-scoring algorithm converges in less than \code{MAXITER} iterations.}
    \item{\code{iterations}: }{number of iterations performed by the Fisher-scoring algorithm.}
   \item{\code{estcoef}: }{a data frame with the estimated model coefficients in the first column (\code{beta}), 
              their asymptotic standard errors in the second column (\code{std.error}), 
              the \eqn{t} statistics in the third column \code{(tvalue)} and the p-values of the significance 
              of each coefficient in last column \code{(pvalue)}.}
    \item{\code{refvar}: }{estimated random effects variance.}
    \item{\code{goodness}: }{vector containing several goodness-of-fit measures: loglikehood, AIC, BIC, KIC and the measures proposed by Marhuenda et al. (2014): AICc, AICb1, AICb2, KICc, KICb1, KICb2. \code{B} must be must be greater than 0 to obtain these last measures.}
    }
  }

In case that \code{formula} or \code{vardir} contain NA values a message is printed and no action is done.
}
\references{

- Fay, R.E. and Herriot, R.A. (1979). Estimation of income from small places: An application of James-Stein procedures to census data. Journal of the American Statistical Association 74, 269-277.

- Marhuenda, Y., Morales, D. and Pardo, M.C. (2014). Information criteria for Fay-Herriot model selection. Computational Statistics and Data Analysis 70, 268-280.

- Rao, J.N.K. (2003). Small Area Estimation. Wiley, London.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{ 
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{mseFH}} }
\examples{
# Load data set  
data(milk)   
attach(milk)

# Fit FH model using REML method with indicators of 4 Major Areas as 
# explanatory variables.
resultREML <- eblupFH(yi ~ as.factor(MajorArea), SD^2)
resultREML

#Fit FH model using FH method
resultFH <- eblupFH(yi ~ as.factor(MajorArea), SD^2, method="FH")
resultFH

detach(milk)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
