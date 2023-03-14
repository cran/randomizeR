#' @include getDesFunc.R
#' @include derFunc.R

NULL

###############################################
# --------------------------------------------#
# Class derringerLs                           #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the derringerLs class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateDerringerLs <- function(object) {
  errors <- character()
  TV <- object@TV
  SLs <- object@SLs
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
    msg <- paste("First element of b is ", b[1], ". Should be greater or equal than zero. ", 
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
  
  if(SLs >= TV){
    msg <- paste("SLs = ", SLs, " is greater than or equal to TV = ", TV, ". SLs should be 
                 smaller than TV. ", sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for derringerLs
# --------------------------------------------

# Left-sided desirabilty parameters generic
setClass("derringerLs",
         slots = c(SLs = "numeric", b = "numeric"), 
         contains = "derFunc",
         validity = validateDerringerLs)

# --------------------------------------------
# Constructor function for derringerLs
# --------------------------------------------

# Representing desirability functions after Derringer and Suich (1980)
# 
# Represents the (left-) one-sided Derringer-Suich desirability function.
#
# @details
# derringerLs represents the framework for left-sided desirability functions introduced 
# by Derringer and Suich (1980). By applying the \code{derringer} function to an object
# of this class the left-sided desirability function is computed. One property of this 
# class is that all \code{x} smaller than the lower specified border \code{LSL} are 
# mapped to zero and all \code{x} greater or equal to the target value \code{TV} are
# mapped to one. For \code{bL = 1} the function for punishment deviations from the target
# value is linear and for \code{bL = 0} every \code{x} in \code{[LSL,TV]} is mapped to one. 
# 
# @family derFunc
#
# @inheritParams overview
# 
# @return 
# \code{S4} object of the class \code{derringerLs}.
#
# @references 
# TBD
# 
# @export
derringerLs <- function(TV, SLs, b = 1) {
  new("derringerLs", TV = TV, SLs = SLs, b = b)
}


# --------------------------------------------
# Definition of (left-) one-sided Derringer-Suich desirability function
# --------------------------------------------
leftDerFunc <- function(TV, LSL, bL, x) {
  stopifnot(is.numeric(TV), is.numeric(LSL), is.numeric(bL), is.numeric(x), bL >= 0)
  d <- 0
  if(x <= LSL){
    d <- 0
  } else if( LSL < x && x < TV ){
    d <- ((x-LSL)/(TV-LSL))^bL  
  } else{
    d <- 1
  }
  d
}


# --------------------------------------------
# Methods for derringerLs
# --------------------------------------------
setMethod("derringer", signature(obj = "derringerLs"),
          function(obj, ...) {
            L <- list(...)
            TV <- TV(obj); LSL <- obj@SLs; bL <- obj@b
            stopifnot(is.numeric(TV), is.numeric(LSL), is.numeric(bL), bL >= 0)
			      stopifnot(all(sapply(unlist(L), function(x)  is(x, "numeric"))))
            d <- sapply(unlist(L), function(x) leftDerFunc(TV, LSL, bL, x))
            d
          }
)


#' @rdname getDesFunc
setMethod("getDesFunc", 
          signature(obj = "derringerLs"),
          function(obj) {
            paste("derringerLs(", round(obj@TV, digits = 2), ", ", round(obj@SLs, digits = 2), ", ", 
                  round(obj@b, digits = 2), ")", sep = "")
          }
)



