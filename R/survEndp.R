#' @include getExpectation.R
#' @include getParameters.R
###############################################
# --------------------------------------------#
# Class survEndp                              #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the endpoint class
#
# @inheritParams overview
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateSurvEndp <- function(object) {
  errors <- character()

  if(!(length(object@cenRate) == 1)) {
    msg <- paste("Censoring rate should have length 1. Has length ", length(object@cenRate),
                 ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }
  if (!(all(object@cenRate > 0))) {
    msg <- ("The common exponential censoring rate must be positive.")
    errors <- c(errors, msg)
  }

  if(!(length(object@accrualTime) == 1)) {
    msg <- paste("Accrual time should have length 1. Has length ", length(object@accrualTime),
                 ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }
  if (!(all(object@accrualTime >= 0))) {
    msg <- ("The accrual time must non-negative.")
    errors <- c(errors, msg)
  }

  if(!(length(object@cenTime) == 1)) {
    msg <- paste("Censoring time should have length 1. Has length ", length(object@cenTime),
                 ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if (!(all(object@cenTime >= 0)) ){
    msg <- ("Censoring time must non-negative.")
    errors <- c(errors, msg)
  }

  if(length(object@cenTime) == 1 && length(object@accrualTime) == 1){
    if (!(object@cenTime > object@accrualTime)) {
      msg <- ("The censoring time must be greater than the accrual time.")
      errors <- c(errors, msg)
    }
  }

  if (!(all(object@weights >= 0))) {
    msg <- ("The Fleming-Harrington weight parameters must be non-negative")
    errors <- c(errors, msg)
  }

  if (!(all(object@shape > 0))) {
    msg <- ("The shape parameter must be positive.")
    errors <- c(errors, msg)
  }

  if (!(all(object@scale > 0))) {
    msg <- ("The scale parameter must be positive.")
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for survEndp
# --------------------------------------------

# Representation of the exponential endpoints
setClass("survEndp",
         slots = c(cenRate="numeric", accrualTime="numeric", cenTime="numeric", weights = "numeric",shape = "numeric", scale = "numeric", maxcombo = "logical"),
         validity = validateSurvEndp)



# --------------------------------------------
# Constructor function for survEndp
# --------------------------------------------

#' Representation of survival endpoints
#'
#' Represents survival endpoints in clinical trials.
#'
#' @inheritParams overview
#' @param shape parameter of the Weibull distribution (must be positive)
#' @param scale parameter of the Weibull distribution (must be positive)
#'
#' @details
#' The \code{survEnd} function is a constructor function
#' for an S4 object of the class \code{survEnd} representing
#' a survival endpoint in a clinical trial.
#'
#' @family endpoint types
#'
#' @seealso Compute exact or simulated type-I-error: \code{\link{assess}}.
#'
#' @returns A \code{S4} object representing
#' a survival endpoint in a clinical trial.
#'
#' @export
survEndp <- function(cenRate, accrualTime = 0, cenTime, shape, scale, weights = c(0,0), maxcombo = FALSE) {
  new("survEndp", cenRate = cenRate, accrualTime = accrualTime, cenTime = cenTime, shape = shape, scale = scale, weights = weights, maxcombo = maxcombo)
}

# --------------------------------------------
# Generic function for survEndp
# --------------------------------------------

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "missing",
                                      endp = "survEndp"),
          function(randSeq, endp) {
            stopifnot(randSeq@K == length(endp@shape))
            stopifnot(randSeq@K == length(endp@scale))
            validObject(randSeq); validObject(endp)
            expectation <- matrix(numeric(0), ncol = ncol(randSeq@M),
                                  nrow = nrow(randSeq@M))
            for(i in 0:(randSeq@K-1)) {
             expectation[randSeq@M == i] <- endp@scale[i+1]*gamma(1+1/endp@shape[i+1])
           }
           expectation
          }
)

# --------------------------------------------
# Get shape and scale parameters of survEndp
# --------------------------------------------

#' @rdname getDistributionPars
setMethod("getDistributionPars", signature(randSeq = "randSeq", issue = "missing",
                                      endp = "survEndp"),
          function(randSeq, endp) {
            stopifnot(randSeq@K == length(endp@shape))
            stopifnot(randSeq@K == length(endp@scale))
            validObject(randSeq); validObject(endp)

            shape <- matrix(numeric(0), ncol = ncol(randSeq@M),
                                  nrow = nrow(randSeq@M))
            scale <- matrix(numeric(0), ncol = ncol(randSeq@M),
                            nrow = nrow(randSeq@M))
            for(i in 0:(randSeq@K-1)) {
              shape[randSeq@M == i] <- endp@shape[i+1]
              scale[randSeq@M == i] <- endp@scale[i+1]
            }
            list(shape = shape, scale = scale)
          }
)
