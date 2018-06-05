#' @include selBias.R
#' @include chronBias.R
#' @include bias.R
#' @include corGuess.R
#' @include imbalance.R
#' @include power.R
NULL

###############################################
# --------------------------------------------#
# Class issue                                 #
# --------------------------------------------#
###############################################

#' Assessment criteria for clinical trials
#' 
#' Summarizes the criteria for the assessment of randomization procedures.
#' 
#' @details
#' Randomization in clinical trials is supposed to 
#' control certain properties in clinical trials. In the randomizeR package, 
#' these properties are called \code{issues}.
#' It is crucial to decide which of the issues is relevant in the present 
#' clinical trial, because a randomization procedure that mananges well one issue
#' might behave very badly for another. The issues include
#' \itemize{
#' \item \strong{Selection bias} 
#' 		can occur if future treatment allocations are predictable due to 
#' 		restricted randomization and unmasking of past treatment assigments.
#'		The influence of selection bias on the test decision is represented by 
#'		the \code{\link{selBias}} class. The measure for the predictability of
#'		a randomization procedure is impemented in the \code{\link{corGuess}} class
#'		representing the expected number of correct guesses.
#' \item \strong{Chronological bias} 
#' 		can occur if a time trend is present in the data. Time trends occur
#' 		due to learning curves, relaxed inclusion/ exclusion criteria or
#'    new co-medication.
#'		Chronological bias is represented by the \code{\link{chronBias}} class.
#' \item \strong{Additive combination of chronological and selection bias} 
#' 		may occur if a time trend and selection bias are present in the data.
#'		The combined bias is represented by the \code{\link{combineBias}} class.
#' \item \strong{Balance}
#' 		is important in order to ensure proper power estimation properties of
#'    the treatments.
#' 		However, a high degree of balance favours selection bias. Depending on the
#' 		clinical context, a randomization procedure should be chosen that admits 
#' 		a suitable imbalance.
#'		Imbalance bias is represented by the \code{\link{imbal}} class. The power
#'		loss due to imbalance can be assessed directly via the \code{\link{setPower}}
#'		class
#' }
#' @name issue
#'
#' @seealso Representation of randomization procedures: \code{\link{randPar}}
#' @seealso Generation of randomization sequences: \code{\link{genSeq}}
#' @seealso Assessment of randomization sequences: \code{\link{assess}}
#' @seealso Comparison of randomization sequences: \code{\link{compare}}
#'
#' @family issues
#' 
#' @aliases issues
NULL

# Issue class
# 
# @name issue
setClassUnion("issue", c("selBias", "chronBias", "corGuess", "imbal", "power"))
#setClassUnion("issue", c("bias", "corGuess", "imbal", "power"))


# --------------------------------------------
# Accesssor functions for issue
# --------------------------------------------

#' Method defining the $ operator for the issue class
#' 
#' @inheritParams overview
setMethod("$", "issue",
          function(x, name) slot(x, name))


# --------------------------------------------
# Show function for issue
# --------------------------------------------

setMethod("show", signature = "issue", definition = function(object) {
  validObject(object)
  # headline
  cat("\n Object of class \"", class(object)[1], "\"\n\n", sep = "")
  # iterate through all slots of the object
  names <- slotNames(object) 
  for(name in names){
    cat("\t", toupper(name), "=", slot(object, name), "\n")
  }
  cat("\n") 
})








