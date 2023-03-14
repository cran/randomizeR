# --------------------------------------------
# Generic function for distribution parameters
# --------------------------------------------

#' Get distribution parameters of a randomization list
#'
#' Generates a matrix of the distribution parameters of the included patients
#' in the clinical trial.
#'
#' @param randSeq object of the class randSeq.
#' @param issue object of the class issue (optional).
#' @param endp object of the class endpoint.
#'
#' @examples
#' # return the shape and scale parameters of a Weibull distribution
#' endp <- survEndp(shape = c(1,1), scale = c(0.5,1), cenTime = 10, cenRate = 0.01)
#' biasSB <- selBias("CS", log(2), "exact")
#' randSeq <- genSeq(rpbrPar(rb = 2, N = 12))
#' getDistributionPars(randSeq,biasSB,endp)
#' @name getDistributionPars
NULL

#' @rdname getDistributionPars
#' @returns a matrix of the distribution parameters of the included
#'  patients in the clinical trial.
#'
#' @export
setGeneric("getDistributionPars", function(randSeq, issue, endp) standardGeneric("getDistributionPars"))
