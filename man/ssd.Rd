\name{ssd}
\alias{ssd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample size dependent estimator.}
\description{
Calculates sample size dependent estimators of domain means, 
as composition of direct and synthetic estimators. The estimators
 involved in the composition must be given as function arguments.
}
\usage{
ssd(dom, sweight, domsize, direct, synthetic, delta = 1, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dom}{ vector or factor (same size as \code{y}) with domain codes. }
  \item{sweight}{ vector (same size as \code{dom}) with sampling weights of the units.}
  \item{domsize}{ matrix or data frame with domain codes in the first column and the corresponding domain population sizes in the second column.}
  \item{direct}{ matrix or data frame with domain codes in the first column and the corresponding direct estimators of domain means in the second column.}
  \item{synthetic}{ matrix or data frame with domain codes in the first column and the corresponding synthetic estimators of domain means in the second column.}
  \item{delta}{ constant involved in sample size dependent estimator,
              controlling how much strength to borrow. Default value is 1. }
  \item{data}{ optional data frame containing the variables named in \code{dom} and \code{sweight}. By default the variables are taken from the environment from which \code{ssd} is called.}
}

\value{
The function returns a data frame of size \code{D*2} with the following columns:
  \item{Domain}{ domain codes in ascending order.}
  \item{ssd }{ sample size dependent estimators of domain means.}
  \item{CompWeight}{ weights attached to direct estimators in the composition.}

Cases with NA values in \code{dom} or \code{sweight} are ignored.  
}
\references{

- Drew, D., Singh, M.P. and Choudhry, G.H. (1982). Evaluation of small area estimation techniques for the Canadian Labour Force Survey. Survey Methodology 8, 17-47.

- Rao, J. N. K. (2003). Small Area Estimation. Wiley, London.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{direct}, \link{pssynt}}
}
\examples{
# We compute sample size dependent estimators of mean income by  
# composition of the Horvitz-Thompson direct estimator and the  
# post-stratified synthetic estimator with age groups as post-strata.

# Load data set
data(incomedata)

# Load population sizes of provinces (domains)
data(sizeprov)  

# First we compute Horvitz-Thompson direct estimators
dir <- direct(y=income, dom=provlab, sweight=weight,
              domsize=sizeprov[,c(1,3)], data=incomedata)

# Now we compute post-stratified synthetic estimators with education 
# levels as post-strata
# Load province sizes by education levels
data(sizeprovedu)

# Compute post-stratified synthetic estimators
colnames(sizeprovedu) <- c("provlab", "prov", "0", "1", "2", "3")
synth <- pssynt(y=income, sweight=weight, ps=educ,
                domsizebyps=sizeprovedu[,-2], data=incomedata)
                
# Compute sample size dependent estimators of province mean income
# by composition of Horvitz-Thompson direct estimators and 
# post-stratified estimators for delta=1
comp <- ssd(dom=provlab, sweight=weight, domsize=sizeprov[,c(1,3)],
            direct=dir[,c("Domain","Direct")], synthetic=synth, data=incomedata)
comp
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
