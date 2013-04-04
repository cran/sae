\name{diagonalizematrix}
\alias{diagonalizematrix}
\title{It constructs a block-diagonal matrix.}
\description{
Using a \code{n*m} matrix \code{A}, this function constructs a block-diagonal matrix with dimension \code{(n*ntimes) *} \code{(m*ntimes)}, with all blocks equal to matrix \code{A} and the rest of entries equal to \code{0}.
}
\usage{
diagonalizematrix(A, ntimes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{\code{n*m} matrix with the values.}
  \item{ntimes}{number of times.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
X <- matrix(data=c(1,2,3,4,5,6), nrow=3, ncol=2)
diagonalizematrix(X,3)
}

\keyword{ method }
