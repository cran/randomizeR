#' @include getExpectation.R
#' @include getStat.R
#' @include issue.R
#' @include randSeq.R
#' @include util.R
#' @include endpoint.R
NULL

###############################################
# --------------------------------------------#
# Class combinedBias                          #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Class definition for combinedBias
# --------------------------------------------

# The combinedBias class
setClass("combinedBias", slots = c(eta = "numeric", typeSB = "character",
                   theta = "numeric", typeCB = "character",
                   method = "character", alpha = "numeric"))
setClass("combinedBiasStepTrend", slots = c(saltus = "numeric"),
         contains = "combinedBias")


#' Combined bias criterion
#'
#' This class combines a \code{selBias} object and a \code{chronBias} object
#' to a new object. In the analysis within the new object the
#' two types of bias are treated as additive effect for normal endpoints
#' and as multiplicative effect for exponential endpoints.
#'
#' @param selBias object of class \code{selBias}
#' @param chronBias object of class \code{chronBias}
#'
#' @family issues
#'
#' @examples
#' chronBias <- chronBias(type="linT", theta=1, method="sim")
#' selBias <- selBias(type="CS", eta=1, method="sim")
#' combineBias(selBias, chronBias)
#' @returns A combined bias object that combines a \code{selBias} and
#' a \code{chronBias} object
#' @export
combineBias <- function(selBias, chronBias) {
  stopifnot(is(selBias, "selBias"), is(chronBias, "chronBias"))
  if (selBias@alpha != chronBias@alpha) {
    warning("Parameter alpha taken from object selBias.")
  }
  if (selBias@method != chronBias@method) {
    warning("Parameter method taken from object selBias.")
  }
  if (!.hasSlot(chronBias, "saltus")) {
    new("combinedBias", eta = selBias@eta, typeSB = selBias@type, theta = chronBias@theta,
        method = selBias@method, alpha = selBias@alpha, typeCB = chronBias@type)
  } else {
    new("combinedBiasStepTrend", eta = selBias@eta, typeSB = selBias@type,
        theta = chronBias@theta,  method = selBias@method,
        alpha = selBias@alpha, saltus = chronBias@saltus, typeCB = chronBias@type)
  }
}


# --------------------------------------------
# Class union of bias and combinedBias
# --------------------------------------------

setClassUnion("bias", c("combinedBias", "combinedBiasStepTrend"))
setClassUnion("issue", c("combinedBias", "combinedBiasStepTrend"))


# --------------------------------------------
# Methods for combinedBias
# --------------------------------------------

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBias",
                                      endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            if(randSeq@K == 2)
              expectationSB <- getExpectation(randSeq, selBias, endp)
            else
              expectationSB <- makeBiasedExpectation(randSeq, endp@mu, selBias)
            expectationCB <- getExpectation(randSeq, chronBias)
            expectationSB + expectationCB
})

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBias",
                                      endp = "expEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@alpha)
            selBias   <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            expectationSB <- getExpectation(randSeq, selBias, endp)
            expectationCB <- getExpectation(randSeq, chronBias, endp)
            # Combined expectation
            expectationCombined <- expectationSB * expectationCB
            expectationCombined[randSeq@M == 0] <- expectationCombined[randSeq@M == 0] * endp@lambda[1]
            expectationCombined[randSeq@M == 1] <- expectationCombined[randSeq@M == 1] * endp@lambda[2]
            expectationCombined
})

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBias",
                                      endp = "survEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            expectationSB <- getExpectation(randSeq, selBias, endp)
            expectationCB <- getExpectation(randSeq, chronBias, endp)
            expectationNB <- getExpectation(randSeq, issue = "missing", endp)
            # Combined expectation
            expectationCombined <- expectationSB * expectationCB / expectationNB
            expectationCombined
})


#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBiasStepTrend",
                                      endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@saltus,
                                   issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            if(randSeq@K == 2)
              expectationSB <- getExpectation(randSeq, selBias, endp)
            else {
              if(dim(randSeq@M)[1]){
                R_ <- genSeq(crPar(randSeq@N, randSeq@K))
                expectationSB <- t(apply(randSeq@M, 1, function(x){
                  R_@M <- matrix(x, ncol = N(randSeq))
                  makeBiasedExpectation(R_, endp@mu, selBias)
                }))
              }else{
                expectationSB <- makeBiasedExpectation(randSeq, endp@mu, selBias)
              }
            }
            expectationCB <- getExpectation(randSeq, chronBias, endp)
            expectationSB + expectationCB
})

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBiasStepTrend",
                                      endp = "expEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2, randSeq@K == length(endp@lambda))
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@saltus,
                                   issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            expectationSB <- getExpectation(randSeq, selBias, endp)
            expectationCB <- getExpectation(randSeq, chronBias, endp)
            # Combined expectation
            expectationCombined <- expectationSB * expectationCB
            expectationCombined[randSeq@M == 0] <- expectationCombined[randSeq@M == 0] * endp@lambda[1]
            expectationCombined[randSeq@M == 1] <- expectationCombined[randSeq@M == 1] * endp@lambda[2]
            expectationCombined
})

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "combinedBiasStepTrend",
                                      endp = "survEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@saltus,
                                   issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            expectationSB <- getExpectation(randSeq, selBias, endp)
            expectationCB <- getExpectation(randSeq, chronBias, endp)
            expectationNB <- getExpectation(randSeq, issue = "missing", endp)
            # Combined expectation
            expectationCombined <- expectationSB * expectationCB / expectationNB
            expectationCombined
})

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "combinedBias",
                               endp = "endpoint"),
          function(randSeq, issue, endp) {
            stopifnot(validObject(randSeq), validObject(endp))
            if (issue@method == "sim") {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("testDec", " ", issue@method, "(combined)", sep = "")
              D
            } else {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("rejection prob.", " ", issue@method, "(combined)",
                                   sep = "")
              D
            }
          }
)



# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "combinedBiasStepTrend",
                               endp = "endpoint"),
          function(randSeq, issue, endp) {
            stopifnot(validObject(randSeq), validObject(endp))
            if (issue@method == "sim") {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("testDec", " ", issue@method, "(combined)", sep = "")
              D
            } else {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("rejection prob.", " ", issue@method, "(combined)", sep = "")
              D
            }
          }
)

# --------------------------------------------
# Get Parameters for combined Bias
# --------------------------------------------

#' @rdname getDistributionPars
setMethod("getDistributionPars", signature(randSeq = "randSeq", issue = "combinedBias",
                                           endp = "survEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            scaleSB <- getDistributionPars(randSeq, selBias, endp)$scale
            scaleCB <- getDistributionPars(randSeq, chronBias, endp)$scale

            # Combined distribution parameters
            shape <- matrix(numeric(0), ncol = ncol(randSeq@M),
                            nrow = nrow(randSeq@M))
            scale <- matrix(numeric(0), ncol = ncol(randSeq@M),
                            nrow = nrow(randSeq@M))
            for(i in 0:(randSeq@K-1)) {
              shape[randSeq@M == i] <- endp@shape[i+1]
              scale[randSeq@M == i] <- scaleSB[randSeq@M == i]*scaleCB[randSeq@M == i]/endp@scale[i+1]
            }
            list(shape = shape, scale = scale)
          }
)

#' @rdname getDistributionPars
setMethod("getDistributionPars", signature(randSeq = "randSeq", issue = "combinedBiasStepTrend",
                                           endp = "survEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(endp)
            chronBias <- chronBias(issue@typeCB, issue@theta, issue@method, issue@saltus,
                                   issue@alpha)
            selBias <- selBias(issue@typeSB, issue@eta, issue@method, issue@alpha)
            scaleSB <- getDistributionPars(randSeq, selBias, endp)$scale
            scaleCB <- getDistributionPars(randSeq, chronBias, endp)$scale

            # Combined distribution parameters
            shape <- matrix(numeric(0), ncol = ncol(randSeq@M),
                            nrow = nrow(randSeq@M))
            scale <- matrix(numeric(0), ncol = ncol(randSeq@M),
                            nrow = nrow(randSeq@M))
            for(i in 0:(randSeq@K-1)) {
              shape[randSeq@M == i] <- endp@shape[i+1]
              scale[randSeq@M == i] <- scaleSB[randSeq@M == i]*scaleCB[randSeq@M == i]/endp@scale[i+1]
            }
            list(shape = shape, scale = scale)
          }
)
