#' Randomization for Clinical Trials
#'
#' This tool enables the user to choose a randomization procedure based on
#' sound scientific criteria. It comprises the generation of randomization
#' sequences as well the assessment of randomization procedures based on
#' carefully selected criteria. Furthermore, randomizeR provides a function for
#' the comparison of randomization procedures.
#'
#' @references
#' W. F. Rosenberger and J. M. Lachin (2002) \emph{Randomization in Clinical Trials}.
#' Wiley.
#'
#' #' D. Uschner, D. Schindler, R. D. Hilgers and N. Heussen (2018). "randomizeR: An R Package
#' for the Assessment and Implementation of Randomization in Clinical Trials." \emph{Journal
#' of Statistical Software}, \strong{85}(8), pp. 1-22. doi: 10.18637/jss.v085.i08
#' .
#'
#' D. Schindler (2016) \emph{Assessment of Randomization Procedures in the Presence of
#' Selection and Chronological Bias}. PhD Thesis.
#'
#' D. Uschner, R. D. Hilgers, N. Heussen (2018). "The impact of selection bias in randomized
#' multi-arm parallel group clinical trials." PLOS ONE, 13(1): e0192065.
#' doi: 10.1371/journal.pone.0192065.
#'
#' M. V. Rueckbeil, R. D. Hilgers, N. Heussen (2019). "Randomization in survival studies: An
#' evaluation method that takes into account selection and chronological bias." PLOS ONE, 14(6):
#' e0217946. doi: 10.1371/journal.pone.0217946.
#'
#'
#' @section Acknowledgement:
#' This research is embedded in the
#' \href{https://www.ideal.rwth-aachen.de/}{IDeAl project}, which has received
#' funding from the European Union's Seventh Framework Programme for
#' research, technological development and demonstration under Grant Agreement
#' no 602552.
#'
#' @seealso
#' For functionality for randomization procedures, see \code{\link{randPar}} and
#' \code{\link{genSeq}}.
#' For the criteria for the assessment of randomization procedures, see
#' \code{\link{issues}}.
#' For the assessment and comparison of randomization procedures, see
#' \code{\link{assess}} and \code{\link{compare}}.
#'
#' @docType package
#' @name randomizeR-package
#' @aliases randomizeR
#' @author David Schindler \email{dv.schindler@@gmail.com}, Diane Uschner \email{Diane.Uschner@@gmail.com},
#' Ralf-Dieter Hilgers, Nicole Heussen, Marcia Viviane Rueckbeil \email{marcia.rueckbeil@@rwth-aachen.de}
#'
#' @import methods
#' @import ggplot2
#' @import plotrix
#' @import survival
#' @importFrom stats dpois pt qpois qt rbinom rnorm t.test
#' @importFrom utils capture.output packageVersion sessionInfo write.table
#' @importFrom graphics abline axis box lines plot.new plot.window title
#' @importFrom grDevices col2rgb rainbow
#' @importFrom graphics legend
#' @importFrom stats pbeta
#' @importFrom stats integrate dexp pchisq pexp pnorm punif qexp qf qnorm rexp runif
NULL
