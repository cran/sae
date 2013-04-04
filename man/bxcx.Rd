\name{bxcx}
\alias{bxcx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Box-Cox Transformation and its Inverse 
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Box-Cox or power transformation or its inverse. For lambda!=0, the Box-Cox transformation of x is (x^lambda-1)/lambda, whereas the regular power transformation is simply x^lambda. When lambda=0, it is log in both cases. The inverse of the Box-Cox and the power transform can also be obtained. 
}
\usage{
bxcx(x, lambda, InverseQ = FALSE, type = "BoxCox")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector or time series}
  \item{lambda}{power transformation parameter}
  \item{InverseQ}{if TRUE, the inverse transformation is done}
  \item{type}{either "BoxCox" or "power"}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A vector or time series of the transformed data 
}
\references{
%% ~put references to the literature/web site here ~
- Box, G.E.P. and Cox, D.R. (1964). An analysis of transformations. Journal of Royal Statistical Society  Series B 26, 211-246.
}
\author{
%%  ~~who you are~~
A.I. McLeod. R package \code{FitAR}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{

#lambda=0.5
z<-AirPassengers; lambda<-0.5
y<-bxcx(z, lambda)
z2<-bxcx(y, lambda, InverseQ=TRUE)
sum(abs(z2-z))
#
z<-AirPassengers; lambda<-0.0
y<-bxcx(z, lambda)
z2<-bxcx(y, lambda, InverseQ=TRUE)
sum(abs(z2-z))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
