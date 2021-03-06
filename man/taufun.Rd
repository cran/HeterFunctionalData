% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{taufun}
\alias{taufun}
\title{An intermediate function to calculate the partial sums}
\usage{
taufun(u, ranks, d1, d2, d3, a, b, mn, mcon, coln)
}
\arguments{
\item{u}{a vector}

\item{ranks}{a vector corresponds to either the fourth column of data or rank of the observations.}

\item{d1}{a vector corresponds to the first column of data in the long format.}

\item{d2}{a vector corresponds to the second column of data in the long format.}

\item{d3}{a vector corresponds to the third column of data in the long format.}

\item{a}{The number of treatments.}

\item{b}{The number of time points or repeated measurements per subject.}

\item{mn}{The vector of sample sizes in treatments.}

\item{mcon}{a vector corresponds to $b^h$ in Theorem 3.2 of Wand and Akritas (2010a)}

\item{coln}{takes value 4 or 5. Value 4 is for the tests based on original data
and value 5 is for the tests based on ranks.}
}
\value{
Asymptotic variances to be used in Heter.test().
}
\description{
This function is internally used to calculate the partial sums which will
then be used in calculating the asymptotic variances of the test statistics in Heter.test()
}
