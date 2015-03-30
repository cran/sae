\name{direct}
\alias{direct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Direct estimators.}
\description{This function calculates direct estimators of domain means.}
\usage{
direct(y, dom, sweight, domsize, data, replace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ vector specifying the individual values of the variable for which we want to estimate the domain means. }
  \item{dom}{ vector or factor (same size as \code{y}) with domain codes. }
  \item{sweight}{ optional vector (same size as \code{y}) with sampling weights. When this argument is not included, by default estimators 
                  are obtained under simple random sampling (SRS). }
  \item{domsize}{ \code{D*2} data frame with domain codes in the first column and the corresponding domain population sizes in the second column. This argument is not required when \code{sweight} is not included and \code{replace=TRUE} (SRS with replacement).}
  \item{data}{ optional data frame containing the variables named in \code{y}, \code{dom} and \code{sweight}. By default the variables are taken from the environment from which \code{direct} is called.}
  \item{replace}{ logical variable with default value \code{FALSE} for random sampling without replacement within each domain is considered 
                  and \code{TRUE} for random sampling with replacement within each domain.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{ 
  The function returns a data frame of size \code{D*5} with the following columns:
  \item{Domain}{ domain codes in ascending order.}
  \item{SampSize}{ domain sample sizes.}
  \item{Direct}{ direct estimators of domain means of variable \code{y}.}
  \item{SD}{ estimated standard deviations of domain direct estimators. If sampling design is SRS or Poisson sampling, estimated variances are unbiased. Otherwise, estimated variances are obtained under the approximation that second order inclusion probabilities are the product of first order inclusion probabilities.}
  \item{CV}{ absolute value of percent coefficients of variation of domain direct estimators.}

Cases with NA values in \code{y}, \code{dom} or \code{sweight} are ignored.
%% ...
}
\references{

- Cochran, W.G. (1977). Sampling techniques. Wiley, New York.

- Rao, J.N.K. (2003). Small Area Estimation. Wiley, London.

- Sarndal, C.E., Swensson, B. and Wretman, J. (1992). Model Assisted Survey Sampling. Springer-Verlag.

}
%%\author{}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pssynt}} for post-stratified synthetic estimator, \code{\link{ssd}} for sample size dependent estimator. 

In case that the sampling design is known, see packages \code{survey} or \code{sampling} for more exact variance estimation.
}
\examples{
# Load data set with synthetic income data for provinces (domains)
data(incomedata)

# Load population sizes of provinces
data(sizeprov)   

# Compute Horvitz-Thompson direct estimator of mean income for each 
# province under random sampling without replacement within each province.
result1 <- direct(y=income, dom=prov, sweight=weight,
                   domsize=sizeprov[,2:3], data=incomedata)
result1

# The same but using province labels as domain codes
result2 <- direct(y=incomedata$income, dom=incomedata$provlab,
                  sweight=incomedata$weight, domsize=sizeprov[,c(1,3)])
result2

# The same, under SRS without replacement within each province.
result3 <- direct(y=income ,dom=provlab, domsize=sizeprov[,c(1,3)],
                  data=incomedata)
result3

# Compute direct estimator of mean income for each province
# under SRS with replacement within each province
result4 <- direct(y=income, dom=provlab, data=incomedata, replace=TRUE)
result4
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
