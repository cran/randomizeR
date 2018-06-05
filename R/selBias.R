#' @include doublyT.R
#' @include getStat.R
#' @include testDec.R
#' @include getExpectation.R
#' @include endpoint.R
NULL

###############################################
# --------------------------------------------#
# Class selBias                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the selBias class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateSelBias <- function(object) {
  errors <- character()
  lengthType <- length(object@type)
  if (lengthType != 1) {
    msg <- paste("Type is length ", lengthType, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  type <- object@type[1]
  if (!(type %in% c("CS", "DS", "CS2"))) {
    msg <- paste("(First) Argument of type is ", type, ". Should be \"CS\", \"DS\" or \"CS2\" ."
                 , sep = "")
    errors <- c(errors, msg)
  }
  
  lengthEta <- length(object@eta)
  if (lengthEta != 1) {
    msg <- paste("eta is length ", lengthEta, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  lengthMethod <- length(object@method)
  if (lengthMethod != 1) {
    msg <- paste("Method is length ", lengthMethod, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  method <- object@method[1]
  if (!(method %in% c("exact", "sim"))) {
    msg <- paste("(First) Argument of method is ", method, ". Should be \"exact\" or \"sim\"."
                 , sep = "")
    errors <- c(errors, msg)
  }
  
  lengthAlpha<- length(object@alpha)
  if (lengthAlpha != 1) {
    msg <- paste("Alpha is length ", lengthAlpha, ". Should be 1.", sep = "")
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
# Class definition for selBias
# --------------------------------------------

# The selBias class
# 
setClass("selBias", slots = c(eta = "numeric", type = "character",
                              method = "character", alpha = "numeric"), 
         validity = validateSelBias)


# --------------------------------------------
# Constructor function for selBias
# --------------------------------------------

#' Representing selection bias
#' 
#' Represents the issue of selection bias in a clinical trial. 
#' 
#' @family issues
#' 
#' @inheritParams overview
#' 
#' @importFrom stats anova lm
#' 
#' @param 
#' method character string, should be one of \code{"sim"} or \code{"exact"}, see 
#' Details.
#' @param
#' type character string, should be one of \code{"CS"}, \code{"CS2"} or \code{"DS"}, see 
#' Details.
#' @param alpha significance level.
#'
#' @details
#' Selection bias can be an issue in the design of a clinical trial. The 
#' \code{selBias} function is a constructor function
#' for an S4 object of the class \code{selBias} representing the issue of
#' third order selection bias in a clinical trial. It supports two possible modes,
#' \code{method="sim"} and \code{method="exact"}. This representation is 
#' particularly useful in interaction with the \code{\link{assess}} function. 
#' 
#' \describe{
#'  \item{\code{method="sim"}}{Represents the simulated type-I-error rate given 
#'  the level \code{alpha}, the selection effect \code{eta} and the biasing 
#'  strategy \code{type}. When calling \code{assess} for a \code{selBias} object 
#'  with \code{method="sim"}, one test decision is computed for each sequence of
#' \code{randSeq}. The type-I-error rate (power) is the proportion of falsely
#' (correctly) rejected null hypotheses.
#'  }
#'  \item{\code{method="exact"}}{Represents the exact type-I-error proabability 
#'  given the level \code{alpha}, the selection effect \code{eta} and the 
#'  biasing strategy \code{type}. When calling \code{assess} for a \code{selBias} 
#'  object with \code{method="exact"}, the exact \emph{p}-value of each 
#'  randomization sequence is computed. So far, this is only supported for
#'  normal endpoints. Then the type-I-error probability is
#'  the sum of the corresponding quantiles of the doubly noncentral t-distribution.
#'  }
#' }
#' 
#' @return
#' \code{S4} object of class \code{selBias}, a formal representation of the
#' issue of selection bias in a clinical trial.
#'
#' @seealso Compute exact or simulated type-I-error: \code{\link{assess}}.
#'
#' @references
#' D. Blackwell and J.L. Hodges Jr. (1957) Design for the control of 
#' selection bias. \emph{Annals of Mathematical Statistics}, \strong{25}, 449-60. 
#' 
#' M. Proschan (1994) Influence of selection bias on the type-I-error rate  
#' under random permuted block designs. \emph{Statistica Sinica}, \strong{4}, 219-31. 
#' 
#' @export
selBias <- function(type, eta, method, alpha = 0.05) { 
  new("selBias", type = type, eta = eta, method = method, alpha = alpha)
}

# --------------------------------------------
# Methods for selBias
# --------------------------------------------

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "selBias", endp = "missing"),
          function(randSeq, issue, endp) stop("Need an object of endpoint class."))

# @rdname getStat
setMethod("getStat", signature(randSeq = "randSeq", issue = "selBias", endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(validObject(randSeq), validObject(issue), validObject(endp), randSeq@K == length(endp@mu))
            if (issue@method == "sim") {
              D <- data.frame(testDec(randSeq, issue, endp))
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
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "selBias",
                                     endp = "normEndp"),
          function(randSeq, issue, endp) {
            stopifnot(randSeq@K == length(endp@mu))
            validObject(randSeq); validObject(issue); validObject(endp)
            # for the second convergence strategy, use the already implemented function
            # in the doublyF file
            if(issue@type == "CS2"){
              R_ <- genSeq(crPar(N(randSeq), K(randSeq)))
              issue_ <- t(apply(randSeq@M, 1, function(x) {
                R_@M <- matrix(x, nrow = 1)
                makeBiasedExpectation(R_, endp@mu, issue)
              }))
              return(issue_)
            } else if (issue@type == "CS") {
              # convergence strategy
              issue <- t(apply(randSeq@M, 1, function(x) {
                issue <- sign(cumsum(2*x - 1)) * issue@eta
                issue <- issue[-length(issue)]
                issue <- c(0, issue)
                issue
              }))
            } else if (issue@type == "DS") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                issue <- sign(cumsum(2*x - 1)) * issue@eta * (-1)
                issue <- issue[-length(issue)]
                issue <- c(0, issue)
                issue
              }))  
            }
            issue[randSeq@M == 0] <- issue[randSeq@M == 0] + endp@mu[1]
            issue[randSeq@M == 1] <- issue[randSeq@M == 1] + endp@mu[2]
            issue
          }
)

#' @rdname getExpectation
setMethod("getExpectation", signature(randSeq = "randSeq", issue = "selBias", 
                                     endp = "missing"), 
          function(randSeq, issue) {
            stopifnot(validObject(randSeq), validObject(issue))
            # if no endpoint is specified, create a vector with 0s for the
            # second convergence strategy
            if(issue@type == "CS2"){
              mu <- rep(0, randSeq@K)
              R_ <- genSeq(crPar(N(randSeq), K(randSeq)))
              issue_ <- t(apply(randSeq@M, 1, function(x) {
                R_@M <- matrix(x, nrow = 1)
                makeBiasedExpectation(R_, mu, issue)
              }))
              return(issue_)
            }else if (issue@type == "CS") {
              # convergence strategy
              issue <- t(apply(randSeq@M, 1, function(x) {
                issue <- sign(cumsum(2*x - 1)) * issue@eta
                issue <- issue[-length(issue)]
                issue <- c(0, issue)
                issue
              }))
            } else if (issue@type == "DS") {
              issue <- t(apply(randSeq@M, 1, function(x) {
                issue <- sign(cumsum(2*x - 1)) * issue@eta
                issue <- issue[-length(issue)]
                issue <- c(0, issue)
                issue
              }))  
            }

            issue
          }
)