#' @include issue.R
#' @include randSeq.R
#' @include util.R
#' @include endpoint.R
#' @include assess.R
#' @include desScores.R
NULL

###############################################
# --------------------------------------------#
# Class ProbUnDes                             #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the desirability class
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateProbUnDes <- function(object) {
  errors <- character()
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for probUnDesirable
# --------------------------------------------

# probUnDesirable parameters generic
setClass("probUnDesirable",
         slots = c(D = "data.frame", design = "character", N = "numeric", K = "numeric", 
                   desFuncs = "character", weights = "numeric"),
         validity = validateProbUnDes)


# --------------------------------------------
# Accesssor functions for probUnDesirable
# --------------------------------------------

#' Method defining the $ operator for the probUnDesirable class
#' @keywords internal
#' @inheritParams overview
setMethod("$", "probUnDesirable",
          function(x, name) slot(x, name))


# --------------------------------------------
# Show function for probUnDesirable
# --------------------------------------------

setMethod("show", "probUnDesirable", function(object) {
  # headline
  cat("\nObject of class \"", class(object)[1],"\"\n\n", sep="")
  # iterate through all slots of the object
  names <- slotNames(object)
  names <- names[!(names == "D")] # without D
  for(name in names) {
    if(is.numeric(slot(object, name))){
      cat(name, "=", round(slot(object, name), digits = 3), "\n")
    } else{
      cat(name, "=", slot(object, name), "\n")
    }
  }
  cat("\n") 
  # Note that this show function is different from show function for assess, evaluate etc.
  # since the matrix has only one row
  object@D <- round(object@D, digits = 3)
  print(object@D)
  cat("\n") 
}  
)


# --------------------------------------------
# Generic functions for using objects of type desScores
# --------------------------------------------

#' Computing the probability of having desirability scores of zero
#' 
#' Computing the probability of having desirability scores of zero for each desirability
#' function applied to an issue. 
#' 
#' @family desirability topics
#' 
#' @param desScore an object of the class \code{desScores}, i.e. an object resulting from
#' applying the function \code{\link{getDesScores}}
#'
#' @details
#'
#' The function \code{probUnDes} expects an object that results from the \code{\link{getDesScores}}
#' function. For each issue it computes the probability that it achieves an undesirable score, 
#' i.e. a desirability score of 0. In doing so, it weights the zero desirability scores 
#' with the probability that the sequence occurs. 
#'
#' @examples 
#' # compare Random Allocation Rule to Big Stick Design with respect to different issues
#' # and their corresponding desirability functions
#' RAR <- getAllSeq(rarPar(4))
#' issue1 <- corGuess("CS")
#' issue2 <- corGuess("DS")
#' A1 <- assess(RAR, issue1, issue2)
#' 
#' d1 <- derFunc(TV = 0.1, 0.7, 2)
#' d2 <- derFunc(0.5, c(0.3, 0.8), c(1, 1))
#' DesScore <- getDesScores(A1, d1, d2, weights = c(5/6, 1/6))
#' 
#' probUnDes(DesScore)
#' 
#'
#' @return
#' \code{S4} object of class \code{probUnDesirable} computing the probability of getting 
#' undesirable scores, i.e. desirability scores of 0. 
#'
#' @seealso Representation of randomization procedures: \code{\link{randPar}}
#' @seealso Generation of randomization sequences: \code{\link{genSeq}}
#' @seealso \code{\link{issues}} for the desirability of randomization sequences
#' 
#' @name probUnDes
NULL

#' @rdname probUnDes
#'
#' @export
setGeneric("probUnDes", function(desScore) standardGeneric("probUnDes"))

# --------------------------------------------
# Methods for probUnDesirable
# --------------------------------------------

#' @rdname probUnDes
setMethod("probUnDes", signature(desScore = "desScores"),
          function(desScore) {
            stopifnot(is(desScore, "desScores"))
            D <- desScore$D[, -c(1,2), drop = FALSE]
            W <- apply(D, 2, function(x) which(x==0))
            if(is.matrix(W)){
              dColnames <- colnames(W)
            } else{
              dColnames <-  names(W)
              max.length <- max(sapply(W, length))
              W <- data.frame(lapply(W, function(v) { c(v, rep(0, max.length-length(v)))}))
            }
            P <- desScore$D[2]
            Z <- data.frame(t(apply(W, 2 , function(x) sum(P[x,]))))
            rownames(Z) <- "1"
            colnames(Z) <- unlist(lapply(dColnames, function(x) paste("P(", x, "=0)", sep = "")))
            new("probUnDesirable", 
                D = Z, design = desScore$design, N = desScore$N, K = desScore$K, 
                desFuncs = desScore$desFuncs, weights = desScore$weights)   
          }
)