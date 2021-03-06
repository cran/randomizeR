% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doublyT.R
\name{doublyTValues}
\alias{doublyTValues}
\title{Calculation of the biased type-one-error probability (resp. power) of Student`s t-test}
\usage{
doublyTValues(randSeq, bias, endp)
}
\arguments{
\item{randSeq}{object of the class randSeq.}

\item{bias}{object of the class bias.}

\item{endp}{object of the class endpoint.}
}
\value{
the biased type-one-error probability (resp. power) of all randomization sequences.
}
\description{
Computes the biased type-one-error probability (resp. power) of Student`s t-test due to shifts in the expectation vectors
in both treatment groups.
}
\keyword{internal}
