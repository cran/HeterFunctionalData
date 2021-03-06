% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{NorRanGen}
\alias{NorRanGen}
\title{Generate a vector of random effects with specific correlation structure and given variance}
\usage{
NorRanGen(m, sigmaj = stats::runif(m, 1.2, 1.4))
}
\arguments{
\item{m}{the length of the vector of the random effects}

\item{sigmaj}{standard deviations for the random effect vector. Default is a vector from U(1.2, 1.4).}
}
\value{
the random effect vector of length m
}
\description{
Generate a vector of random effects that follow AR(1)
correlation structure with rho=exp(-1/m) and variances sigmaj being
given or andomly generated from uniform distribution on (1.2, 1.4).
}
\examples{
   m=50;  raneff=NorRanGen(m)
# Note: If X ~ N(0, I), then tran X ~ N(0, A) with
# A being the cov matrix of AR(1), which contains the standard deviations sigma and the
# correlation coeff rho=exp(-1/m).
# i.e. corr= (1  rho  rho^2 ... rho^(m-1)
#             rho 1   rho   ... rho^(m-2)
#                 ...................
#             rho^(m-1) rho^(m-2) ... rho )
#
# To see the correlation values, run the following example
# j1=seq(25); cv=numeric()
# for (j in 1:25){
#  lag=abs(j1-j)/25; cv=rbind(cv, exp(-lag))
#}
# row.names(cv)=j1;  colnames(cv)=j1;  cv[1,]
}
