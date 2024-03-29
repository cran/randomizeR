% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randPar.R, R/pbrPar.R, R/rarPar.R, R/abcdPar.R,
%   R/bbcdPar.R, R/ebcPar.R, R/bsdPar.R, R/chenPar.R, R/crPar.R, R/gbcdPar.R,
%   R/hadaPar.R, R/mpPar.R, R/tbdPar.R, R/udPar.R
\name{generateAllSequences}
\alias{generateAllSequences}
\alias{getAllSeq}
\alias{getAllSeq,pbrPar-method}
\alias{getAllSeq,rarPar-method}
\alias{getAllSeq,abcdPar-method}
\alias{getAllSeq,bbcdPar-method}
\alias{getAllSeq,ebcPar-method}
\alias{getAllSeq,bsdPar-method}
\alias{getAllSeq,chenPar-method}
\alias{getAllSeq,crPar-method}
\alias{getAllSeq,gbcdPar-method}
\alias{getAllSeq,hadaPar-method}
\alias{getAllSeq,mpPar-method}
\alias{getAllSeq,tbdPar-method}
\alias{getAllSeq,udPar-method}
\title{Complete set of randomization sequences}
\usage{
getAllSeq(obj)

\S4method{getAllSeq}{pbrPar}(obj)

\S4method{getAllSeq}{rarPar}(obj)

\S4method{getAllSeq}{abcdPar}(obj)

\S4method{getAllSeq}{bbcdPar}(obj)

\S4method{getAllSeq}{ebcPar}(obj)

\S4method{getAllSeq}{bsdPar}(obj)

\S4method{getAllSeq}{chenPar}(obj)

\S4method{getAllSeq}{crPar}(obj)

\S4method{getAllSeq}{gbcdPar}(obj)

\S4method{getAllSeq}{hadaPar}(obj)

\S4method{getAllSeq}{mpPar}(obj)

\S4method{getAllSeq}{tbdPar}(obj)

\S4method{getAllSeq}{udPar}(obj)
}
\arguments{
\item{obj}{object specifying the randomization procedure, see \code{\link{randPar}}
or \code{\link{createParam}}.}
}
\value{
An object inheriting from \linkS4class{randSeq}, representing the set
of randomization sequences for the given parameters.
The output consists of the parameters used for the generation of the
randomization sequences (see \code{\link{createParam}}) and the matrix \code{M}
that stores the randomization sequences in its rows.
}
\description{
Computes all randomization sequences for the given randomization procedure,
and stores them in an object along with the parameters belonging to the
randomization procedure.
}
\details{
\code{getAllSeq} is a generic function which dispatches different
methods depending on the type of input. The set of sequences of a procedure
is computed by enumerating all possible sequences and eliminating those that
are not possible in the randomization procedure specified by \code{obj}. The
parameters of the randomization procedure are saved along with the sequences
to ensure reproducibility of the results.
}
\examples{
# all randomization sequences of Efron's Biased Coin Design with p = 0.667 for N = 6
myPar <- ebcPar(6, 0.667)
getAllSeq(myPar)

# all randomization sequences of Big Stick Design with mti = 2 for N = 6
myPar <- bsdPar(6, 2)
getAllSeq(myPar)

# all randomization sequences of Permuted Block Randomization with block sizes 4 and 2
myPar <- pbrPar(c(4, 2))
getAllSeq(myPar)



}
\seealso{
\code{\link{createParam}}
}
