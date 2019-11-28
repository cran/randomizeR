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
  if (length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for survEndp
# --------------------------------------------

# Representation of the exponential endpoints
setClass("survEndp", 
         slots = c(cenRate="numeric", accrualTime="numeric", cenTime="numeric"),
         validity = validateSurvEndp)



# --------------------------------------------
# Constructor function for survEndp
# --------------------------------------------

#' Representation of survival endpoints
#' 
#' Represents survival endpoints in clinical trials.
#'
#' @inheritParams overview
#'
#' @details
#' The \code{survEnd} function is a constructor function
#' for an S4 object of the class \code{survEnd} representing 
#' a survival endpoint in a clinical trial.
#'
#' @family endpoint types
#' 
#' @export
survEndp <- function(cenRate, accrualTime, cenTime) {
  new("survEndp", cenRate = cenRate, accrualTime = accrualTime, cenTime = cenTime)
}
