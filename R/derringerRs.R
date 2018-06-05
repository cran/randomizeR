#' @include getDesFunc.R
#' @include derFunc.R
NULL

###############################################
# --------------------------------------------#
# Class derringerRs                           #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the derringerRs class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateDerringerRs <- function(object) {
  errors <- character()
  SLs <- object@SLs
  TV <- object@TV
  b <- object@b
  
  if(any(is.na(SLs))){
    msg <- "SLs contains NA values. Should be numeric."
    errors <- c(errors, msg)
  }
  
  if(any(is.na(b))){
    msg <- "b contains NA values. Should be numeric."
    errors <- c(errors, msg)
  }
  
  if(b[1] < 0) {
    msg <- paste("First element of b is ", b, ". Should be greater or equal than zero. ", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(SLs) != 1) {
    msg <- paste("SLs has length  ", length(SLs), ". Should be one. ", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(b) != 1) {
    msg <- paste("b has length  ", length(b), ". Should be one. ", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(SLs <= TV){
    msg <- paste("SLs = ", SLs, " is smaller than or equal to TV = ", TV, ". SLs should be 
                 greater than TV. ", sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for derringerRs
# --------------------------------------------

# Right-sided desirabilty parameters generic
setClass("derringerRs",
         slots = c(SLs = "numeric", b = "numeric"),
         contains = "derFunc",
         validity = validateDerringerRs)

# --------------------------------------------
# Constructor function for derringerRs
# --------------------------------------------

# Representing desirability functions after Derringer and Suich (1980)
# 
# Represents the (left-) one-sided Derringer-Suich desirability function.
#
# @details
# derringerLs represents the framework for right-sided desirability functions introduced 
# by Derringer and Suich (1980). By applying the \code{derringer} function to an object
# of this class the right-sided desirability function is computed. One property of this 
# class is that all \code{x} greater or equal than the upper specified border \code{SLs} 
# are mapped to zero and all \code{x} smaller or equal to the target value \code{TV} are
# mapped to one. For \code{b = 1} the function for punishment deviations from the target
# value is linear and for \code{b = 0} every \code{x} in \code{[TV,SLs]} is mapped to one. 
# 
# @family derFunc
# 
# @inheritParams overview
# 
# @return 
# \code{S4} object of the class \code{derringerRs}.
#
# @references 
# TBD
# 
# @export
derringerRs <- function(TV, SLs, b = 1) {
  new("derringerRs", TV = TV, SLs = SLs, b = b)
}

# --------------------------------------------
# Definition of (right-) one-sided Derringer-Suich desirability function
# --------------------------------------------
rightDerFunc <- function(TV, USL, bR, x) {
  stopifnot(is.numeric(TV), is.numeric(USL), is.numeric(bR), is.numeric(x), bR >= 0)
  d <- 0
 if(x >= USL){
   d <- 0
 } else if(TV < x && x < USL){
    d <- ((USL - x)/(USL - TV))^bR  
 } else{
    d <- 1
 }
 d
}

# --------------------------------------------
# Methods for derringerRs
# --------------------------------------------
setMethod("derringer", signature(obj = "derringerRs"),
          function(obj, ...) {
            L <- list(...)
            TV <- TV(obj); USL <- obj@SLs; bR <- obj@b
            stopifnot(is.numeric(TV), is.numeric(USL), is.numeric(bR), bR >= 0)
			      stopifnot(all(sapply(unlist(L), function(x)  is(x, "numeric"))))
            d <- sapply(unlist(L), function(x) rightDerFunc(TV, USL, bR, x))
            d
          }
)

#' @rdname getDesFunc
setMethod("getDesFunc", 
          signature(obj = "derringerRs"),
          function(obj) {
            paste("derringerRs(", round(obj@TV, digits = 2), ", ", round(obj@SLs, digits = 2), ", ", 
                  round(obj@b, digits = 2), ")", sep = "")
          }
)



