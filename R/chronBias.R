#' @include doublyT.R
#' @include getStat.R
#' @include testDec.R
#' @include getExpectation.R
#' @include endpoint.R
NULL

###############################################
# --------------------------------------------#
# Class chronBias                             #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the chronBias class
#
# @inheritParams overview
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateChronBias <- function(object) {
  errors <- character()
  lengthType <- length(object@type)
  if (lengthType != 1) {
    msg <- paste("type is length ", lengthType, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }

  type <- object@type[1]
  if (!(type %in% c("linT", "stepT", "logT"))) {
    msg <- paste("(First) Argument of type is ", type, ". Should be in linT,
                 logT, or stepT.", sep = "")
    errors <- c(errors, msg)
  }

  lengthMethod <- length(object@method)
  if (lengthMethod != 1) {
    msg <- paste("method is length ", lengthMethod, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }

  method <- object@method[1]
  if (!(method %in% c("exact", "sim"))) {
    msg <- paste("(First) Argument of method is ", method, ". Should be exact or sim."
                 , sep = "")
    errors <- c(errors, msg)
  }

  lengthAlpha<- length(object@alpha)
  if (lengthAlpha != 1) {
    msg <- paste("alpha is length ", lengthAlpha, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }

  alpha <- object@alpha[1]
  if (!(alpha >= 0 && alpha <= 1 )) {
    msg <- paste("(First) Argument of alpha is ", alpha, ". Should be in [0,1]."
                 , sep = "")
    errors <- c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for chronBias
# --------------------------------------------

# Randomization parameters generic
setClass("chronBias",
         slots = c("type" = "character", "theta" = "numeric",
                   method = "character", alpha = "numeric"),
         validity = validateChronBias)

# --------------------------------------------
# Constructor function for chronBias
# --------------------------------------------

#' Representing chronological bias
#'
#' Represents the issue of chronological bias in a clinical trial.
#'
#' @inheritParams overview
#' @param method character string, should be one of \code{"sim"} or \code{"exact"}, see Description.
#' @param type character string, should be one of "\code{linT}", "\code{logT}", or "\code{stepT}",
#' see Details.
#' @param alpha  significance level
#'
#'
#' @details
#' Chronological bias can be an issue in the design of a clinical trial. The
#' \code{chronBias} function is a constructor function
#' for an S4 object of the class \code{chronBias} representing the issue of
#' chronological bias, s.a. time trends, in a clinical trial. It supports two possible modes,
#' \code{method="sim"} and \code{method="exact"}, and three different types of trend.
#'
#' If \code{method="sim"}, the object represents the simulated type-I-error rate given
#'  the level \code{alpha}, the selection effect \code{eta} and the biasing
#'  strategy \code{type}. When calling \code{assess} for a \code{chronBias} object
#'  with \code{method="sim"}, one test decision is computed for each sequence of
#' \code{randSeq}. The type-I-error rate (power) is the proportion of falsely
#' (correctly) rejected null hypotheses.
#'
#' If \code{method="exact"}, the object represents the exact type-I-error probability
#'  given the level \code{alpha}, the selection effect \code{eta} and the
#'  biasing strategy \code{type}. When calling \code{assess} for a \code{chronBias}
#'  object with \code{method="exact"}, the \emph{p}-value of each randomization
#'  sequence is computed. For normal endpoints and two treatment groups these p-values
#'  are exact values which can be calculated from the sum of the corresponding quantiles
#'  of the doubly noncentral t-distribution. For more than two treatment groups, exact
#'  p-values are computed using a doubly noncentral F distribution. For exponential
#'  endpoints the p-values are obtained using an approximation formula.
#'
#' \subsection{Types of chronological bias}{
#' \describe{
#' 	 \item{\code{type = "linT"}}{
#'    Represents linear time trend. Linear time trend means that the time trend function of the patients,
#'    i.e. expected response for normal endpoints, increases evenly by \code{theta/(N-1)} with
#'    every patient included in the study, until reaching \code{theta} after \code{N} patients.
#'    Linear time trend may occur as a result of gradually relaxing in- or exclusion criteria
#'    throughout the trial.
#'    It can be represented by the formula:
#'    \deqn{f(i) = (i-1)/(N-1) \theta}{f(i) = (i-1)/(N-1) \theta}
#'   }
#' 	 \item{\code{type = "logT"}}{
#'    Represents logarithmic time trend. Logarithmic time trend means that the time trend function of
#'    the patients, i.e. expected response for normal endpoints, increases logarithmically in the
#'    patient index by \code{theta/log(N)} with every patient included in the study, until reaching
#'    \code{theta} after \code{N} patients. Logarithmic time trend may occur as a result of a learning
#'    curve, i.e. in a surgical trial.
#'    It can be represented by the formula:
#'    \deqn{\log(i)/\log(N) \theta}{f(i) = log(i)/log(N) \theta}
#'   }
#' 	 \item{\code{type = "stepT"}}{
#'    Represents step trend. Step trend means that the expected response of the patients increases
#'    by \code{theta} after a given point (\code{"saltus"}) in the allocation process.
#'    Step trend may occur if a new device is used after the point \eqn{c} = \code{"saltus"}, or if
#'    the medical personal changes after this point.
#'    Step time trend can be represented by the formula:
#'   \deqn{f(i) = 1_{c < i \leq N} \theta}{f(i) = 1_{c < i \le N} \theta}
#'   }
#' }
#' }
#'
#' @return
#' \code{S4} object of class \code{chronBias}, a formal representation of the
#' issue of chronological bias in a clinical trial.
#'
#' @references
#' G. K. Rosenkranz (2011) The impact of randomization on the analysis of
#' clinical trials. \emph{Statistics in Medicine}, \strong{30}, 3475-87.
#'
#' M. Tamm and R.-D. Hilgers (2014) Chronological bias in randomized clinical
#' trials under different types of unobserved time trends.
#' \emph{Methods of Information in Medicine}, \strong{53}, 501-10.
#'
#' @family issues
#'
#' @examples
#' # create a linear time trend with theta = 0.5 for which the exact rejection probabilities
#' # are calculated
#' cbias <- chronBias("linT", 0.5, "exact")
#'
#' # create a stepwise time trend with theta = 1 after 10 allocations for which the test
#' # decision is simulated
#' cbias <- chronBias("stepT", 1, "sim", 10)
#'
#' @export
chronBias <- function(type, theta, method, saltus, alpha = 0.05) {
  if(missing(saltus) && type %in% c("linT", "logT")) {
    new("chronBias", type = type, theta = theta, method = method,  alpha = alpha)
  } else {
    new("chronBiasStepT", type = type, theta = theta, method = method, saltus = saltus, alpha = alpha)
  }
}


# --------------------------------------------
# Methods for chronBias
# --------------------------------------------

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "chronBias", endp = "missing"),
          function(randSeq, issue, endp) stop("Need an object of endpoint class."))

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "chronBias", endp = "endpoint"),
          function(randSeq, issue, endp) {
            validObject(randSeq); validObject(issue); validObject(endp)
            if (issue@method == "sim") {
              D <- data.frame(testDec = testDec(randSeq, issue, endp))
              colnames(D) <- paste("testDec(", issue@type, ")", sep = "")
              D
            } else {
              D <- data.frame(testDec(randSeq, issue, endp))
              colnames(D) <- paste("P(rej)(", issue@type, ")", sep = "")
              D
            }
          }
)

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "chronBias", endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(issue); validObject(endp)
            n <- N(randSeq)
            # linear time trend
            if (issue@type == "linT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                0:(n-1)/(n-1) * issue@theta
              }))
            }
			# logarithmic time trend
			else if (issue@type == "logT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                log(1:n)/log(n) * issue@theta
              }))
            }
			# step time trend
			else if (issue@type == "stepT") {
              stopifnot(randSeq@N > issue@saltus)
              issue <- t(apply(randSeq@M, 1, function(x) {
                c(rep(0, issue@saltus), rep(issue@theta, n - issue@saltus))
              }))
            }

            issue[randSeq@M == 0] <- issue[randSeq@M == 0] + endp@mu[1]
            issue[randSeq@M == 1] <- issue[randSeq@M == 1] + endp@mu[2]
            issue
          }
)

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "chronBias", endp = "expEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2, randSeq@K == length(endp@lambda))
            validObject(randSeq); validObject(issue); validObject(endp)
            n <- N(randSeq)
            # linear time trend
            if (issue@type == "linT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                0:(n-1)/(n-1) * issue@theta
              }))
            }
            # logarithmic time trend
            else if (issue@type == "logT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                log(1:n)/log(n) * issue@theta
              }))
            }
            # step time trend
            else if (issue@type == "stepT") {
              stopifnot(randSeq@N > issue@saltus)
              issue <- t(apply(randSeq@M, 1, function(x) {
                c(rep(0, issue@saltus), rep(issue@theta, n - issue@saltus))
              }))
            }
            issue[randSeq@M == 0] <- exp(issue[randSeq@M == 0]) * endp@lambda[1]
            issue[randSeq@M == 1] <- exp(issue[randSeq@M == 1]) * endp@lambda[2]
            1/issue
          }
)

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "chronBias", endp = "survEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(issue); validObject(endp)
            n <- N(randSeq)
            # linear time trend
            if (issue@type == "linT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                0:(n-1)/(n-1) * issue@theta
              }))
            }
            # logarithmic time trend
            else if (issue@type == "logT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                log(1:n)/log(n) * issue@theta
              }))
            }
            # step time trend
            else if (issue@type == "stepT") {
              stopifnot(randSeq@N > issue@saltus)
              issue <- t(apply(randSeq@M, 1, function(x) {
                c(rep(0, issue@saltus), rep(issue@theta, n - issue@saltus))
              }))
            }
            issue[randSeq@M == 0] <- exp(issue[randSeq@M == 0])^{-1/endp@shape[1]} *
              endp@scale[1] * gamma(1+1/endp@shape[1])
            issue[randSeq@M == 1] <- exp(issue[randSeq@M == 1])^{-1/endp@shape[2]} *
              endp@scale[2] * gamma(1+1/endp@shape[2])
            issue
          }
)




#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "chronBias", endp = "missing"),
          function(randSeq, issue) {
            #stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(issue)
            n <- N(randSeq)
            # linear time trend
            if (issue@type == "linT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                0:(n-1)/(n-1) * issue@theta
              }))
            # logarithmic time trend
            } else if (issue@type == "logT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                log(1:n)/log(n) * issue@theta
              }))
            # step time trend
            } else if (issue@type == "stepT") {
              stopifnot(randSeq@N > issue@saltus, is(issue,"chronBiasStepT"))
              issue <- t(apply(randSeq@M, 1, function(x) {
                c(rep(0, issue@saltus), rep(issue@theta, n - issue@saltus))
              }))
            }
            issue
          }
)
# --------------------------------------------
# Get Parameters for chronBias
# --------------------------------------------

#' @rdname getDistributionPars
setMethod("getDistributionPars", signature(randSeq = "randSeq", issue = "chronBias",
                                           endp = "survEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == 2)
            validObject(randSeq); validObject(issue); validObject(endp)
            n <- N(randSeq)
            # linear time trend
            if (issue@type == "linT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                0:(n-1)/(n-1) * issue@theta
              }))
            }
            # logarithmic time trend
            else if (issue@type == "logT") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                log(1:n)/log(n) * issue@theta
              }))
            }
            # step time trend
            else if (issue@type == "stepT") {
              stopifnot(randSeq@N > issue@saltus)
              issue <- t(apply(randSeq@M, 1, function(x) {
                c(rep(0, issue@saltus), rep(issue@theta, n - issue@saltus))
              }))
            }
            shape <- matrix(numeric(0), ncol = ncol(randSeq@M),
                            nrow = nrow(randSeq@M))
            scale <- matrix(numeric(0), ncol = ncol(randSeq@M),
                            nrow = nrow(randSeq@M))
            for(i in 0:(randSeq@K-1)) {
              shape[randSeq@M == i] <- endp@shape[i+1]
              scale[randSeq@M == i] <- exp(issue[randSeq@M == i])^{-1/endp@shape[i+1]} * endp@scale[i+1]
            }
            list(shape = shape, scale = scale)
          }
)
