% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survEndp.R
\name{survEndp}
\alias{survEndp}
\title{Representation of survival endpoints}
\usage{
survEndp(
  cenRate,
  accrualTime = 0,
  cenTime,
  shape,
  scale,
  weights = c(0, 0),
  maxcombo = FALSE
)
}
\arguments{
\item{cenRate}{exponential censoring rate in a survival study.}

\item{accrualTime}{duration of the accrual period in a survival study.}

\item{cenTime}{total duration of a survival study (maximum length of followup).}

\item{shape}{parameter of the Weibull distribution (must be positive)}

\item{scale}{parameter of the Weibull distribution (must be positive)}

\item{weights}{numeric specifying the weights used for the test. Unless specified an unweighted test is conducted.}

\item{maxcombo}{logical specifying if the maxcombo test is used}
}
\value{
A \code{S4} object representing
a survival endpoint in a clinical trial.
}
\description{
Represents survival endpoints in clinical trials.
}
\details{
The \code{survEnd} function is a constructor function
for an S4 object of the class \code{survEnd} representing
a survival endpoint in a clinical trial.
}
\seealso{
Compute exact or simulated type-I-error: \code{\link{assess}}.

Other endpoint types: 
\code{\link{expEndp}},
\code{\link{normEndp}}
}
\concept{endpoint types}
