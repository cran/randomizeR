#' @include derFunc.R
#' @include desFunc.R
NULL
# --------------------------------------------
# Generic function for Desirability Functions
# --------------------------------------------

#' Type of Desirability function
#' 
#' Generates a \code{character} vector which specifies the used desirability function and
#' its parameters
#' 
#' @param obj object of the class \code{desFunc}.
#'
#' @name getDesFunc
NULL


#' @rdname getDesFunc
#' 
#' @export
setGeneric("getDesFunc", function(obj) standardGeneric("getDesFunc"))