# @include desFunc.R
# @include derringerLs.R
# @include derringerRs.R
# @include derringerTs.R
NULL
###############################################
# --------------------------------------------#
# Class derFunc                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check for objects of the desirability functions class
# 
# @inheritParams overview
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateDerFunc <- function(object) {
  errors <- character()
  TV <- object@TV
  
  if(length(TV) != 1) {
    msg <- paste("TV has length ", length(TV), ". Should be one.", 
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }
  
  if(any(is.na(TV))){
    msg <- "Tv contains NA values. Should be numeric."
    errors <- c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for derFunc
# --------------------------------------------

# Parameters generic
setClass("derFunc", 
         slots = c(TV = "numeric"),
         validity = validateDerFunc)

# --------------------------------------------
# Constructor function for derFunc
# --------------------------------------------
#' Representing Derringer-Suich desirability functions
#' 
#' Represents the Derringer-Suich desirability approach.
#' 
#' @family desirability topics
#' 
#' @inheritParams overview
#' 
#' @param TV numeric specifying the optimal desired value called the target value.
#' 
#' @details 
#' derFunc represents the framework for left, right and two-sided desirability functions 
#' introduced by Derringer and Suich (1980). For all three different kinds of desirability
#' functions the parameter \code{TV} must be specified. If the parameter \code{SLs} has 
#' length 1, either the left- or right-sided desirability function is created depending
#' from whether the value is smaller (left-sided) or greater (right-sided) than the target 
#' value. By specifying \code{SLs} as a vector of length 2 a two-sided  desirability 
#' function is created where the lower specified border is determined as the smaller value 
#' of \code{SLs} and thus the upper specified border is determined as the greater value. 
#' If there are no values specified for the weights, then they are automatically set to 1 
#' (linear loss). \cr
#' 
#' @examples 
#' # create an object of a left-sided desirability function
#' dLeft <- derFunc(0.5, 0.3, 2)
#' 
#' # create an object of a right-sided desirability function
#' dRight <- derFunc(0.5, 0.8, 1)
#' 
#' # create an object of a two-sided desirability function
#' dLR <- derFunc(0.5, c(0.3, 0.9), c(3, 1))
#' 
#' @return 
#' \code{S4} object of class \code{derFunc}, a formal representation of desirability 
#' functions introduced by Derringer and Suich.
#' 
#' @references 
#' Derringer, G., and Suich, R., (1980) Simultaneous Optimization of Several Response 
#' Variables. \emph{Journal of Quality Technology}, \strong{12}, 214-219. 
#' 
#' 
#' @export
derFunc <- function(TV, SLs, b){
  bMiss <- missing(b)
  if(length(SLs) == 1 && SLs <= TV){
    new("derringerLs", TV = TV, SLs = SLs, b = if(bMiss) 1 else b)
  }
  else if(length(SLs) == 1 && SLs > TV){
    new("derringerRs", TV = TV, SLs = SLs, b = if(bMiss) 1 else b)
  }else{
    if(bMiss) b <- c(1,1)
    else b <- c(b[which(SLs <= TV)], b[which(SLs > TV)])
    new("derringerTs", TV = TV, SLs = sort(SLs, na.last = TRUE), b = b[!is.na(b)])
  }
}



# --------------------------------------------
# Show function for derFunc
# --------------------------------------------

setMethod("show", "derFunc", function(object) {
  validObject(object)
  # headline
  cat("\nObject of class \"", class(object)[1], "\"\n\n", sep = "")
  # crop the method from the class name of the derFunc object
  cat("desirability function =", getDesFunc(object), "\n") 
  # iterate through all slots of the randPar object
  names <- slotNames(object)
  for(name in names) {
    cat(name, "=", slot(object, name),"\n")
  }
  cat("\n") 
})


# --------------------------------------------
# Accesssor functions for derFunc
# --------------------------------------------

#' Method defining the $ operator for the derFunc class
#' @keywords internal
#' @inheritParams overview
setMethod("$", "derFunc",
          function(x, name) slot(x, name))

#' Function returning the target value slot of an S4 object
#'
#' @param obj object inheriting from derFunc 
#' 
#' @export
TV <- function(obj) {
  if (.hasSlot(obj, "TV")) {
    obj@TV
  } else {
    stop("object has no slot named TV.")
  }
}

# --------------------------------------------
# Generic functions for derFunc
# --------------------------------------------

# Derringer-Suich desirability function
# 
# @param obj object of the class \code{derFunc}
# @param ... at least one argument of class \code{numeric} or just a list of objects
# of the class \code{numeric}
#
# @details 
# \code{derringer} computes the corresponding Derringer-Suich desirability function for
# an argument \code{x}. If the value of \code{x} corresponds to its optimal desired value 
# the desirability function is one. 
# One property of this class is that e.g. for two-sided desirabilty functions all \code{x} 
# greater than the upper specified border \code{USL} or smaller than the lower specified 
# border \code{LSL} are mapped to zero. For \code{bR = bL = 1} the function for punishment 
# deviations from the target value is linear and for \code{bR = bL = 0} every \code{x} in 
# \code{[LSL,USL]} is mapped to one. 
# 
# @inheritParams overview
#  
# @return The result of the applied Derringer-Suich desirability function.
# 
# @name derringer
# NULL
setGeneric("derringer", function(obj, ...) standardGeneric("derringer"))
