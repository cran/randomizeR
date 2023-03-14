#' @include derFunc.R
NULL

###############################################
# --------------------------------------------#
# Class desirability                          #
# --------------------------------------------#
###############################################

#' Desirability functions within the scope of clinical trials
#' 
#' Illustrates the interplay between functions related to desirability indices.
#' 
#' @details 
#' Currently, \code{randomizeR} encompasses the class of desirability functions introduced
#' by Derringer and Suich (1980) and corresponding functions to evaluate and compare 
#' randomization sequences which have been assessed on the basis of desirability indices 
#' of specific issues:
#' \itemize{
#' \item \strong{\link{derFunc}}
#'    represents the class of desirability functions according to Derringer-Suich (1980).
#' \item \strong{\link{getDesScores}}
#'    can be applied to an object of class \code{assessment} together with prespecified
#'    desirability functions to compare the behavior of randomization sequences (on a 
#'    common scale \[0,1\]).
#' \item \strong{\link{plotDes}}
#'    plots a \code{desScores} object on a radar chart. 
#' \item \strong{\link{evaluate}}
#'    performs a comparison of sequences from different randomization sequences on the 
#'    basis of object of the class \code{desScores}.
#' \item \strong{\link{plotEv}}
#'    plots an \code{evaluation} object on a radar chart. 
#' \item \strong{\link{probUnDes}}
#'    computes the probability of undesired randomization sequences with respect to 
#'    certain issues and desirability functions.
#' }
#' 
#' @examples 
#' # perform a comparison of randomization sequences from different randomization procedures 
#' # with the help of desirability functions
#' 
#' issue1 <- corGuess("CS")
#' issue2 <- chronBias(type = "linT", theta = 1/4, method = "exact")
#' RAR <- getAllSeq(rarPar(4))
#' BSD <- getAllSeq(bsdPar(4, mti = 2))
#' A1 <- assess(RAR, issue1, issue2, endp = normEndp(c(0,0), c(1,1)))
#' A2 <- assess(BSD, issue1, issue2, endp = normEndp(c(0,0), c(1,1)))
#' 
#' d1 <- derFunc(TV = 0.5, 0.75, 2)
#' d2 <- derFunc(0.05, c(0, 0.1), c(1, 1))
#'
#' # apply the getDesScores function to the assessment output with the specified desirability
#' # functions to evaluate the behaviour of randomization sequences on a [0,1] scale
#' 
#' DesScore <- getDesScores(A1, d1, d2, weights = c(5/6, 1/6))
#' DesScore2 <- getDesScores(A2, d1, d2, weights = c(5/6, 1/6))
#' 
#' # plotting the desScores objects
#' plotDes(DesScore, quantiles = TRUE)
#' plotDes(DesScore2, quantiles = TRUE)
#' 
#' # summarize the results of getDesScore with respect to the statistic "mean"
#' evaluate(DesScore, DesScore2)
#' 
#' # plot the evaluation objects for a visualized comparison
#' plotEv(evaluate(DesScore, DesScore2))
#' 
#' # display which randomzation procedure produces more undesired randomization sequences 
#' # with respect to certain issues and desirability functions
#' probUnDes(DesScore)
#' probUnDes(DesScore2)
#' 
#' @name desirability
NULL
