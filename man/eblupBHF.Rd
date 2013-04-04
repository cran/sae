\name{eblupBHF}
\alias{eblupBHF}
\title{EBLUPs of domain means based on a nested error linear regression model.}
\description{ This function calculates, for selected domains, EBLUPs of domain means based on the 
nested error linear regression model of Battese, Harter and Fuller (1988).}
\usage{
eblupBHF(formula, dom, selectdom, meanxpop, popnsize, method = "REML", data)
}

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
  \item{method}{a character string. If \code{"REML"}, the model is fitted by maximizing the restricted 
                log-likelihood. If \code{"ML"} the log-likelihood is maximized. Defaults to \code{"REML"}.}
  \item{data}{optional data frame containing the variables named in \code{formula} and \code{dom}. By default the variables are taken from the environment from which \code{eblupBHF} is called.}                
            
}
\details{
A typical model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} 
for more details of allowed formulae.
}

  
\value{ The function returns a list with the following objects:
  \item{eblup }{data frame with number of rows equal to number of selected domains (\code{selectdom}),
                          containing in its columns the domain codes (\code{domain}) and the EBLUPs 
                          of the means of selected domains based on the nested error linear regression 
                          model (\code{eblup}). For domains with zero sample size, the EBLUPs are the synthetic regression estimators.}
  \item{fit}{a list containing the following objects: 
  \itemize{
    \item{\code{summary}: }{summary of the unit level model fitting.}
    \item{\code{fixed}: }{vector with the estimated values of the fixed regression coefficient.}
    \item{\code{random}: }{vector with the predicted random effects.}
    \item{\code{errorvar}: }{estimated model error variance.}
    \item{\code{refvar}: }{estimated random effects variance.}
    \item{\code{loglike}: }{log-likelihood.}
    \item{\code{residuals}: }{vector with raw residuals.}
    }  
  }    

Cases with NA values in \code{formula} or \code{dom} are ignored. 
}

\references{ 

- Battese, G.E., Harter, R.M. and Fuller, W.A. (1988). An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data, Journal of the American Statistical Association 83, 28-36. 

- Rao, J.N.K. (2003). Small Area Estimation. New York: John Wiley and Sons.
 }
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}
%%\section{Warning }{ MIRAR } 
\seealso{ \code{\link{pbmseBHF}} }
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

# Compute EBLUPs of county means of corn crop areas for all counties
resultCorn <- eblupBHF(CornHec ~ CornPix + SoyBeansPix, dom=County, 
                       meanxpop=Xmean, popnsize=Popn, data=cornsoybean)
resultCorn$eblup

# Compute EBLUPs of county means of soy beans crop areas for
# a subset of counties using ML method
domains <- c(10,1,5)
resultBean <- eblupBHF(SoyBeansHec ~ CornPix + SoyBeansPix, dom=County, 
                       selectdom=domains, meanxpop=Xmean, popnsize=Popn,
                       method="ML", data=cornsoybean)
resultBean$eblup
resultBean$fit

detach(cornsoybeanmeans)
}

\keyword{ method }
