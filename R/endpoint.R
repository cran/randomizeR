#' @include normEndp.R
#' @include expEndp.R
NULL

###############################################
# --------------------------------------------#
# Class endpoint                              #
# --------------------------------------------#
###############################################

# All endpoints should be added to this class union

# Common representation of the endpoints.
#
# @name endpoint
setClassUnion("endpoint", c("normEndp","expEndp"))


# --------------------------------------------
# Accessor functions for endpoints
# --------------------------------------------

#' Access the expectation value slot of a normEndp S4 object
#' 
#' @param obj object of class normEndp
mu <- function(obj) {
  if (.hasSlot(obj, "mu")) obj@mu else stop("object has no slot named mu.")  
}

#' Function returning the standard deviation slot of a normEndp S4 object
#' 
#' @param obj object of class normEndp
sigma <- function(obj) {
  if (.hasSlot(obj, "sigma")) obj@sigma else stop("object has no slot named sigma.") 
}

#' Method returning the rate parameter of an expEndp S4 object
#' 
#' @param obj object of class expEndp
lambda <- function(obj) {
  if (.hasSlot(obj, "lambda")) obj@lambda else stop("object has no slot named lambda.") 
}

#' Method defining the $ operator for the endpoint class
#' @keywords  internal
#' @inheritParams overview
setMethod("$", "endpoint",
          function(x, name) slot(x, name))


# --------------------------------------------
# Show function for endopoints
# --------------------------------------------

setMethod("show", signature = "endpoint", definition = function(object){
  validObject(object)
  # headline
  cat("\n Object of class \"", is(object)[1], "\"\n\n", sep="")
  # iterate through all slots of the object
  names <- slotNames(object)
  for(name in names){
    cat("\t", name, "=", slot(object, name), "\n")
  }
  cat("\n")
})
