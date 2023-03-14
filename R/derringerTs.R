#' @include getDesFunc.R
#' @include derFunc.R
NULL

###############################################
# --------------------------------------------#
# Class derringerTs                           #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the derringerTs class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateDerringerTs <- function(object) {
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
  
  if(length(SLs) != length(b)){
    msg <- paste("SLs and b have different lengths. Should be equal. ")
    errors <- c(errors, msg)
  }
  
  if(any(b < 0, na.rm = TRUE)) {
    msg <- paste("b is ", b, ". Each element of b should be greater or equal than zero. ", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(length(SLs) != 2) {
    msg <- paste("SLs has length  ", length(SLs), ". Should be two. ", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(any(SLs == TV, na.rm = TRUE)){
    msg <- sapply(which(SLs == TV), 
                  function(x) paste("SLs[", x, "] = ", SLs[x], " is equal to TV = ", TV, 
                                    ". Should be smaller or greater than TV.", sep = ""))
    errors <- c(errors, msg)
  }
  
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for derringerTs
# --------------------------------------------

# Two-sided desirabilty parameters generic
setClass("derringerTs",
         slots = c(SLs = "numeric", b = "numeric"),
         contains = "derFunc",
         validity = validateDerringerTs)

# --------------------------------------------
# Constructor function for derringerTs
# --------------------------------------------

# Representing desirability functions after Derringer and Suich (1980)
# 
# Represents the (left-) one-sided Derringer-Suich desirability function.
#
# @details
# derringerLs represents the framework for two-sided desirability functions introduced 
# by Derringer and Suich (1980). By applying the \code{derringer} function to an object
# of this class the two-sided desirability function is computed. One property of this 
# class is that all \code{x} greater than the upper specified border \code{SLs} or 
# smaller than the lower specified border \code{LSL} are mapped to zero. For 
# \code{bR = b = 1} the function for punishment deviations from the target value is 
# linear and for \code{bR = b = 0} every \code{x} in \code{[LSL,SLs]} is mapped to one. 
# 
# @family derFunc
# 
# @inheritParams overview
# 
# @return 
# \code{S4} object of the class \code{derringerTs}.
#
# @references 
# TBD
# 
# @export
derringerTs <- function(TV, SLs, b = c(1,1)) {
  new("derringerTs", TV = TV, SLs = sort(SLs), b = c(b[which(SLs <= TV)], b[which(SLs > TV)]))
}

# --------------------------------------------
# Definition of two-sided Derringer-Suich desirability function
# --------------------------------------------
derringTsFunc <- function(TV, LSL, USL, bL, bR, x) {
  stopifnot(is.numeric(TV), is.numeric(LSL), is.numeric(USL), is.numeric(bL),
            is.numeric(bR), is.numeric(x), bL >= 0, bR >= 0)
  d <- 0
  if(x > USL || x < LSL){
    d <- 0
  } else if(TV < x && x <= USL){
    d <- ((USL - x)/(USL - TV))^bR  
  } else{
    d <- ((x - LSL)/(TV - LSL))^bL  
  }
  d
}


# --------------------------------------------
# Methods for derringerTs
# --------------------------------------------
setMethod("derringer", signature(obj = "derringerTs"),
          function(obj, ...) {
            L <- list(...)
            TV <- TV(obj)
            SLs <- obj@SLs
            # ind1 <- which(SLs < TV)
            # ind2 <- which(SLs > TV)
            LSL <- SLs[1]; USL <- SLs[2]; bL <- obj@b[1]; bR <- obj@b[2]
            stopifnot(is.numeric(TV), is.numeric(USL), is.numeric(bL), is.numeric(bR),
                      bL >= 0, bR >= 0)
			      stopifnot(all(sapply(unlist(L), function(x)  is(x, "numeric"))))
            d <- sapply(unlist(L), function(x) derringTsFunc(TV, LSL, USL, bL, bR, x))
            d
          }
)

#' @rdname getDesFunc
setMethod("getDesFunc", 
          signature(obj = "derringerTs"),
          function(obj) {
            paste("derringerTs(", round(obj@TV, digits = 2), ", ", 
                  round(obj@SLs[1], digits = 2), ", ", 
                  round(obj@SLs[2], digits = 2), ", ",
                  round(obj@b[1], digits = 2), ", ", 
                  round(obj@b[2], digits = 2), ")",
                  sep = "")
          }
)



