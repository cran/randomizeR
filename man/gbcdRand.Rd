% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gbcdPar.R
\name{gbcdRand}
\alias{gbcdRand}
\title{Sampling algorithm for gbcd}
\usage{
gbcdRand(N, rho, K = 2)
}
\arguments{
\item{N}{integer for the total sample size of the trial.}

\item{rho}{nonnegative parameter which my be adjusted according to how strongly it is
desired to balance the experiment. If \code{rho = 1}, we have Wei's urn design with
\code{alpha = 0}. If \code{rho = 0}, we have complete randomization.}

\item{K}{number of treatment groups (e.g. K=2 if we compare one experimental
against one control treatment).}
}
\value{
A vector with the allocation sequence for a clinical trial. 
It will contain a zero (resp. 1) at position \code{i}, when patient \code{i}
is allocated to treatment A (resp. B).
}
\description{
Sampling algorithm for gbcd
}
\references{
R. L. Smith (1984) Sequential treatment allocation using biased coin designs. 
\emph{Journal of the Royal Statistical Society B},
\strong{46}, 519-543. \cr
W. F. Rosenberger and J. M. Lachin (2002) Randomization in Clinical Trials. \emph{Wiley},
64-65
}
\keyword{internal}
