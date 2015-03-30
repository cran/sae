\name{incomedata}
\alias{incomedata}
\docType{data}
\title{ Synthetic income data.}
\description{ Synthetic data on income and other related variables for Spanish provinces.}
\usage{data(incomedata)}
\format{
  A data frame with 17199 observations on the following 21 variables.
  \describe{
    \item{\code{provlab}:}{ province name.}
    \item{\code{prov}:}{ province code.}
    \item{\code{ac}:}{ region of the province.}
    \item{\code{gen}:}{ gender: 1:male, 2:female.}
    \item{\code{age}:}{ age group: 0:<=13, 1:14-15, 2:16-24, 3:25-49, 4:50-64, 5: >=65.}
    \item{\code{nat}:}{ nationality: 1:Spanish, 2:other.}
    \item{\code{educ}:}{ education level: 0:age<16, 1:primary education (compulsory educ.), 2:secondary education, 3:post-secondary education.}
    \item{\code{labor}:}{ labor force status: 0:age<16, 1:employed, 2:unemployed, 3:inactive.}
    \item{\code{age2}:}{ indicator of age group 16-24.}
    \item{\code{age3}:}{ indicator of age group 25-49.}
    \item{\code{age4}:}{ indicator of age group 50-64.}
    \item{\code{age5}:}{ indicator of age group >=65.}
    \item{\code{educ1}:}{ indicator of education level 1 (primary education).}
    \item{\code{educ2}:}{ indicator of education level 2 (secondary education.}
    \item{\code{educ3}:}{ indicator of education level 3 (post-secondary education).}
    \item{\code{nat1}:}{ indicator of Spanish nationality.}
    \item{\code{labor1}:}{ indicator of being employed.}
    \item{\code{labor2}:}{ indicator of being unemployed.}
    \item{\code{labor3}:}{ indicator of being inactive.}
    \item{\code{income}:}{ normalized income.}
    \item{\code{weight}:}{ sampling weight.}
  }
}
%\details{ There are not missing values}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\keyword{datasets}
