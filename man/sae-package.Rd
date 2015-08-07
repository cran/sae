\name{sae-package}
\alias{sae-package}
\alias{sae}
\docType{package}
\title{Small area estimation}
\description{
This package provides a variety of functions for small area estimation, including functions for mean squared error estimation. Basic estimators include direct, poststratified synthetic and sample size dependent. Model-based estimators include the EBLUP based on a Fay-Herriot model and the EBLUP based on a unit level nested error model. Estimators obtained from spatial and spatio-temporal Fay-Herriot models and the EB method based on the unit level nested error model for estimation of general non linear parameters are also included.}
\details{
 This package provides functions for estimation in domains with small sample sizes. For a complete list of functions, see library(help=sae).
\tabular{ll}{
Package: \tab sae\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-08-07\cr
License: \tab GPL-2\cr
Depends: \tab stats, nlme, MASS\cr
}

%%~~ An overview of how to use the package, including the most important ~~
%%~~ functions ~~
}
\author{
	Isabel Molina <isabel.molina@uc3m.es> and 
	Yolanda Marhuenda <y.marhuenda@umh.es>	
}
\references{

- Arora, V. and Lahiri, P. (1997). On the superiority of the Bayesian method over the BLUP in small area estimation problems. Statistica Sinica 7, 1053-1063.

- Battesse, G.E., Harter, R.M. and Fuller, W.A. (1988). An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data. Journal of the American Statistical Association 83, 28-36.

- Box, G.E.P. and Cox, D.R. (1964). An analysis of transformations. Journal of Royal Statistical Society  Series B 26, 211-246. 

- Cochran, W.G. (1977). Sampling techniques. Wiley, New York.

- Datta, G.S. and Lahiri, P. (2000). A unified measure of uncertainty of estimated best linear unbiased predictors in small area estimation problems. Statistica Sinica 10, 613-627.

- Datta, G.S., Rao, J.N.K. and Smith D.D. (2005). On measuring the variability of small area estimators under a basic area level model. Biometrika 92, 183-196.

- Drew, D., Singh, M.P. and Choudhry, G.H. (1982). Evaluation of small area estimation techniques for the Canadian Labour Force Survey. Survey Methodology 8, 17-47.

- Fay, R.E. and Herriot, R.A. (1979). Estimation of income from small places: An application of James-Stein procedures to census data. Journal of the American Statistical Association 74, 269-277.

- Gonzalez-Manteiga, W., Lombardia, M., Molina, I., Morales, D. and Santamaria, L. (2008). Analytic and bootstrap approximations of prediction errors under a multivariate Fay-Herriot model. Computational Statistics and Data Analysis 52, 5242-5252.

- Jiang, J. (1996). REML estimation: asymptotic behavior and related topics.
Annals of Statistics 24, 255-286.

- Marhuenda, Y., Molina, I. and Morales, D. (2013). Small area estimation with spatio-temporal Fay-Herriot models. Computational Statistics and Data Analysis 58, 308-325. 

- Marhuenda, Y., Morales, D. and Pardo, M.C. (2014). Information criteria for Fay-Herriot model selection. Computational Statistics and Data Analysis 70, 268-280. 

- Molina, I., Salvati, N. and Pratesi, M. (2009). Bootstrap for estimating the MSE of the Spatial EBLUP. Computational Statistics 24, 441-458.

- Molina, I. and Rao, J.N.K. (2010). Small Area Estimation of Poverty Indicators. The Canadian Journal of Statistics 38, 369-385.

- Petrucci, A. and Salvati, N. (2006). Small area estimation for spatial correlation in watershed erosion assessment.
Journal of Agricultural, Biological and Environmental Statistics 11, 169-182.

- Prasad, N. and Rao, J. (1990). The estimation of the mean squared error of small-area estimators. Journal of the American Statistical Association 85, 163-171.

- Pratesi, M. and Salvati, N. (2008). Small area estimation: the EBLUP estimator based on spatially correlated random area effects. 
Statistical Methods & Applications 17, 113-141.

- Rao, J.N.K. (2003). Small Area Estimation. Wiley, London.

- Sarndal, C.E., Swensson, B. and Wretman, J. (1992). Model Assisted Survey Sampling. Springer-Verlag.

- Singh, B., Shukla, G. and Kundu, D. (2005). Spatio-temporal models in small area estimation. Survey Methodology 31, 183-195.

- Small Area Methods for Poverty and Living Conditions Estimates (SAMPLE), funded by European Commission, Collaborative Project 217565, Call identifier FP7-SSH-2007-1.

- You, Y. and Chapman, B. (2006). Small area estimation using area level models and estimated sampling variances. Survey Methodology 32, 97-103.
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%%~~ the R documentation directory ~~
%%\keyword{ package }
%%\seealso{
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%%}
%%\examples{
%%## simple examples of the most important functions 
%%}
