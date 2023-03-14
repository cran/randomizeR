# --------------------------------------------
# Generic function for Expectations
# --------------------------------------------

#' Get expectations of a randomization list
#'
#' Generates a matrix of the expectations of the included patients in the
#' clinical trial.
#'
#' @param randSeq object of the class randSeq.
#' @param issue object of the class issue (optional).
#' @param endp object of the class endpoint (optional).
#'
#' @details
#' It is assumed that the expectations of the included patients in a clinical trial
#' can be influenced in three different ways:
#' \itemize{
#' \item The strength of selection bias and the guessing strategy
#' of the investigator (see \code{\link{selBias}}).
#' \item The strength of a linear time trend, which is described by an object
#'  of the class \code{\link{chronBias}}.
#' \item The expectations of the investigated treatment groups can be different
#' (see e.g. \code{\link{normEndp}}).
#' }
#'
#' @examples
#' # get Expectation for a normal endpoint
#' myPar <- bsdPar(10, 2)
#' M <- genSeq(myPar, 2)
#' cs <- selBias("CS", 2, "sim")
#' endp <- normEndp(mu = c(2, 2), sigma = c(1, 1))
#' getExpectation(M, cs, endp)
#'
#' # get Expectation for an exponential endpoint
#' cs <- selBias("CS", 0.1 , "sim")
#' endp <- expEndp(lambda = c(0.5, 1), cenTime = 10, cenRate = 0.01)
#' getExpectation(M, cs, endp)
#'
#' @name getExpectation
NULL


#' @rdname getExpectation
#' @returns A matrix of the expectations of the included patients in the clinical trial.
#' @export
setGeneric("getExpectation", function(randSeq, issue, endp) standardGeneric("getExpectation"))


