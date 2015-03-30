\name{pssynt}
\alias{pssynt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Post-stratified synthetic estimators of domain means.}
\description{
Calculates post-stratified synthetic estimators of domain means 
using the categories of a cualitative variable as post-strata.
}
\usage{
pssynt(y, sweight, ps, domsizebyps, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ vector specifying the individual values of the variable for which we want to estimate the domain means. }
  \item{sweight}{ vector (same size as \code{y}) with the sampling weights of the units.}
  \item{ps}{ vector (same size as \code{y}) of factor with post-strata codes.}
  \item{domsizebyps}{ data frame with domain codes in the first column. Each remaining column contains the domain population sizes for each post-strata. Names of these columns must be the post-strata identifiers specified in \code{ps}.}
   \item{data}{ optional data frame containing the variables named in \code{y}, \code{sweight} and \code{ps}. By default the variables are taken from the environment from which \code{pssynt} is called.}  
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{ The function returns a data frame of size \code{D*2} with the following columns:
  \item{Domain}{ domain codes in ascending order.}
  \item{PsSynthetic }{ post-stratified synthetic estimators of domain means of variable \code{y}.}
  
  Cases with NA values in \code{y}, \code{sweight} or \code{ps} are ignored.
}
\references{ 

- Rao, J.N.K. (2003). "Small Area Estimation". Wiley, London.

}
%%\author{Isabel Molina, Yolanda Marhuenda}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{direct}, \link{ssd}}
}
\examples{
# Compute post-stratified synthetic estimators of mean income 
# for provinces considering the education levels codes 
# (variable educ) as post-strata.

# Load data set
data(incomedata)   

# Load province sizes by education levels
data(sizeprovedu)

# Compute post-stratified synthetic estimators with province labels 
# as domain codes
colnames(sizeprovedu) <- c("provlab", "prov", "0", "1", "2", "3")
result1 <- pssynt(y=income, sweight=weight, ps=educ,
                  domsizebyps=sizeprovedu[,-2], data=incomedata)
result1

# Now with province codes as domain codes
colnames(sizeprovedu) <- c("provlab", "prov", "0", "1", "2", "3")
result2 <- pssynt(y=income, sweight=weight, ps=educ,
                  domsizebyps=sizeprovedu[,-1], data=incomedata)
result2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
