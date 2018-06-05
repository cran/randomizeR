#' @include issue.R
#' @include randSeq.R
#' @include util.R
#' @include endpoint.R
#' @include assess.R
#' @include desScores.R
NULL

###############################################
# --------------------------------------------#
# Class Evaluation                            #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validity check function for objects of the desirability class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateEvaluation <- function(object) {
  errors <- character()
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for evaluation
# --------------------------------------------

# Evaluation paramters generic
setClass("evaluation",
         slots = c(D = "data.frame", desFuncs = "character", weights = "numeric", 
                   statistic = "character"),
         validity = validateEvaluation)


# --------------------------------------------
# Accesssor functions for evaluation
# --------------------------------------------

#' Method defining the $ operator for the evaluation class
#' 
#' @inheritParams overview
setMethod("$", "evaluation",
          function(x, name) slot(x, name))


# --------------------------------------------
# Show function for desirability
# --------------------------------------------

setMethod("show", "evaluation", function(object) {
  # headline
  cat("\nObject of class \"", class(object)[1],"\"\n\n", sep="")
  # iterate through all slots of the object
  names <- slotNames(object)
  names <- names[!(names == "D")] # without D
  for(name in names) {
    if(is.numeric(slot(object, name))){
      cat(name, "=", round(slot(object, name), digits = 3), "\n")
    } else{
      cat(name, "=", slot(object, name), "\n")
    }
  }
  cat("\n") 
  # The data.frame D is printed seperately dependent on its size.
  if (dim(object@D)[1] <= 3) {
    if (nchar(as.character(object@D[1, 1])) >= 10)
      object@D[ ,1] <- paste(substr(object@D[, 1], 1, 9), "...")
    object@D[ ,-1] <- round(object@D[ ,-1], digits = 3)
    print(object@D) 
  } else {
    cat("\nThe first 3 rows of", dim(object@D)[1], "rows of D: \n\n")
    object <- object@D[1:3, ]
    if (nchar(as.character(object[1, 1])) >= 10)
      object[ ,1] <- paste(substr(object[, 1], 1, 9), "...")
    object[ ,-1] <- round(object[ ,-1],digits = 3)
    print(object) 
    cat("...")
  }
  cat("\n") 
}  
)


# --------------------------------------------
# Generic functions for using objects of type desScores
# --------------------------------------------

#' Evaluation of several randomization procedures with respect to certain desirability
#' functions applied to specified issues.
#'
#'
#' @inheritParams overview
#' 
#' @param ... at least one object of the class \code{desScores} or a list of objects of 
#' the class \code{desScores}.
#' @param statistic character string that specifies on the basis of which statistic the 
#' \code{evaluate} function should be applied. The statistic can be chosen from "mean", 
#' "median", "min" or "max". 
#' 
#' @details
#'
#' The \code{evaluate} function allows the user to compare and evaluate different 
#' randomization procedures. It expects a number of objects that result when applying the 
#' \code{getDesScores} function to an assess object and specified desirability functions. 
#' The \code{evaluate} function summarizes the desirability scores of each randomization 
#' procedure on the basis of a prespecified statistic and encorporates them into a data 
#' frame. If no statistic is specified then it is automatically set to \code{mean}. If 
#' the function is applied to only one object it corresponds simply to 
#' \code{summary(getDesScores(...))}.
#'
#' @examples 
#' # Compare Random Allocation Rule to Big Stick Design with respect to different issues
#' # and their corresponding desirability functions
#' issue1 <- corGuess("CS")
#' issue2 <- corGuess("DS")
#' RAR <- getAllSeq(rarPar(4))
#' BSD <- getAllSeq(bsdPar(4, mti = 2))
#' A1 <- assess(RAR, issue1, issue2)
#' A2 <- assess(BSD, issue1, issue2)
#' 
#' d1 <- derFunc(TV = 0.1, 0.7, 2)
#' d2 <- derFunc(0.5, c(0.3, 0.8), c(1, 1))
#' DesScore <- getDesScores(A1, d1, d2, weights = c(5/6, 1/6))
#' DesScore2 <- getDesScores(A2, d1, d2, weights = c(5/6, 1/6))
#' 
#' evaluate(DesScore, DesScore2)
#' evaluate(DesScore, DesScore2, statistic = "max")
#' 
#'
#' @return
#' \code{S4} object of class \code{evaluation} Comparison of randomization procedures 
#' with respect to desirability functions applied to specified issues, summarized by a
#' prespecified statistic. 
#' 
#' @references 
#' D. Schindler \emph{Assessment of Randomization Procedures in the Presence of 
#' Selection and Chronological Bias}. PhD Thesis.
#'
#' @seealso Representation of randomization procedures: \code{\link{randPar}}
#' @seealso Generation of randomization sequences: \code{\link{genSeq}}
#' @seealso \code{\link{issues}} for the desirability of randomization sequences
#' 
#' @name evaluate
NULL

#' @rdname evaluate
#'
#' @family desirability topics
#'
#' @export
setGeneric("evaluate", function(..., statistic) standardGeneric("evaluate"))

# --------------------------------------------
# Methods for Evaluation
# --------------------------------------------

#' @rdname evaluate
setMethod("evaluate", signature(statistic = "missing"),
          function(...) {
            dScores <- list(...)
            if(length(dScores) == 1 && is.list(dScores[[1]])){
              dScores <- c(...)
            }
            
            stopifnot(all(sapply(dScores, function(x) is(x, "desScores"))))
            n <- ncol(dScores[[1]]$D)
            if(!all(unlist(lapply(dScores, function(x) all(n == ncol(x$D)))))){
              stop("Error: Objects have different number of columns!")
            }
            colnames <- colnames(dScores[[1]]$D)
            if(!all(unlist(lapply(dScores, function(x) all(colnames == colnames(x$D)))))){
              stop("Error: Colnames do not coincide!")
            }
            desFuncs <- dScores[[1]]$desFuncs
            if(!all(unlist(lapply(dScores, function(x) all(desFuncs == x$desFuncs))))){
              warning("The desirability functions do not coincide. The show function only 
                      displays the desirability functions of the first desScores object.")
            }
            weights <- dScores[[1]]$weights
            if(!all(unlist(lapply(dScores, function(x) all(weights == x$weights))))){
              warning("The weights do not coincide. The show function only 
                      displays the weights of the first desScores object.")
            }
            
            # Creates the first column which contains the designs of the different 
            # randomization procedures
            D <- data.frame("RandProc" = sapply(dScores, function(x) x@design))
            # Uses summary(.) to generate the means of the desirability functions and
            # puts it in one matrix
            M <- lapply(dScores, function(x) summary(x)[1,])
            M <- do.call(rbind, M)
            D <- cbind(D, M)           
            
            new("evaluation", 
                D = D, desFuncs = dScores[[1]]$desFuncs, weights = dScores[[1]]$weights, 
                statistic = "mean")   
          }
)

#' @rdname evaluate
setMethod("evaluate", signature(statistic = "character"),
          function(..., statistic) {
            dScores <- list(...)
            if(length(dScores) == 1 && is.list(dScores[[1]])){
              dScores <- c(...)
            }
            
            stats <- c("mean", "median", "min", "max")
            stopifnot(statistic %in% stats)
            stopifnot(all(sapply(dScores, function(x) is(x, "desScores"))))
            n <- ncol(dScores[[1]]$D)
            if(!all(unlist(lapply(dScores, function(x) all(n == ncol(x$D)))))){
              stop("Error: Objects have different number of columns!")
            }
            colnames <- colnames(dScores[[1]]$D)
            if(!all(unlist(lapply(dScores, function(x) all(colnames == colnames(x$D)))))){
              stop("Error: Colnames do not coincide!")
            }
            desFuncs <- dScores[[1]]$desFuncs
            if(!all(unlist(lapply(dScores, function(x) all(desFuncs == x$desFuncs))))){
              warning("The desirability functions do not coincide. The show function only 
                      displays the desirability functions of the first desScores object.")
            }
            weights <- dScores[[1]]$weights
            if(!all(unlist(lapply(dScores, function(x) all(weights == x$weights))))){
              warning("The weights do not coincide. The show function only 
                      displays the weights of the first desScores object.")
            }
            
            
            # Creates the first column which contains the designs of the different 
            # randomization procedures
            D <- data.frame("RandProc" = sapply(dScores, function(x) x@design))
            # colStats contains the row number of the corresponding statistics
            colStats <- c(1, 7, 4, 3)
            # ind determines which statistic was chosen
            ind <- which(statistic == stats)
            # Uses summary(.) to generate the statistic of the desirability functions and
            # puts it in one matrix
            M <- lapply(dScores, function(x) summary(x)[colStats[ind],])
            M <- do.call(rbind, M)
            D <- cbind(D, M)           
            
            new("evaluation", 
                D = D, desFuncs = dScores[[1]]$desFuncs, weights = dScores[[1]]$weights, 
                statistic = statistic)   
            }
)


# --------------------------------------------
# Plot function for evaluate Object
# --------------------------------------------

#' Evaluation plotting
#'
#' Plot of an \code{evaluation} object.
#' 
#' @family desirability topics
#' 
#' @param evaluation object of type \code{evaluation}.
#' @param labels labels used in the plot. Can be \code{NULL}.
#' @param cols colors of the lines representing the  desirability scores in the plot. Can be \code{NULL}.
#' 
#' @examples 
#' # Compare Random Allocation Rule to Big Stick Design with respect to different issues
#' # and their corresponding desirability functions
#' issue1 <- corGuess("CS")
#' issue2 <- chronBias(type = "linT", theta = 1/4, method = "exact")
#' RAR <- getAllSeq(rarPar(4))
#' BSD <- getAllSeq(bsdPar(4, mti = 2))
#' A1 <- assess(RAR, issue1, issue2, endp = normEndp(c(0,0), c(1,1)))
#' A2 <- assess(BSD, issue1, issue2, endp = normEndp(c(0,0), c(1,1)))
#' 
#' d1 <- derFunc(TV = 0.5, 0.75, 2)
#' d2 <- derFunc(0.05, c(0, 0.1), c(1, 1))
#' DesScore <- getDesScores(A1, d1, d2, weights = c(5/6, 1/6))
#' DesScore2 <- getDesScores(A2, d1, d2, weights = c(5/6, 1/6))
#' 
#' E <- evaluate(DesScore, DesScore2)
#' plotEv(E)
#' 
#' @export
plotEv <- function(evaluation , labels, cols) { 
  stopifnot(is(evaluation, "evaluation" ))
  
  if (missing(labels)) {
    labels <- rev(colnames(evaluation@D)[-1])
  } else {
    if (length(labels) != ncol(evaluation@D)-1) stop(paste("Length of labels must be ", ncol(evaluation@D)-1), ".", sep = "")
  }
  
  if (missing(cols)) {
    cols <- rainbow(nrow(evaluation@D))
  } else {
    if (length(cols) != nrow(evaluation@D)) stop(paste("Length of cols must be ", ncol(evaluation@D)-1), ".", sep = "")
  }
  
  addTrans <- function(color, trans){
    # This function adds transparancy to a color.
    # Define transparancy with an integer between 0 and 255
    # 0 being fully transparant and 255 being fully visable
    # Works with either color and trans a vector of equal length,
    # or one of the two of length 1.
    if (length(color) != length(trans) &! any(c(length(color), length(trans)) == 1)) stop("Vector lengths not correct")
    if (length(color) == 1 & length(trans) > 1) color <- rep(color, length(trans))
    if (length(trans) == 1 & length(color) > 1) trans <- rep(trans, length(color))
    num2hex <- function(x){
      hex <- unlist(strsplit("0123456789ABCDEF", split = ""))
      return(paste(hex[(x-x%%16)/16+1], hex[x%%16+1], sep = ""))
    }
    rgb <- rbind(col2rgb(color), trans)
    res <- paste("#", apply(apply(rgb, 2, num2hex), 2, paste, collapse = ""), sep = "")
    return(res)
  }
  
  # backround colour
  bgGrid <- addTrans("grey", 100)
  
  values <- as.numeric(evaluation@D[1,])[-1]
  radial.plot(rev(values), rp.type = "p", start = pi/2, clockwise = TRUE,
              labels = labels, 
              radial.lim = c(0,1), line.col = cols[1], lwd = 4,
              grid.bg = bgGrid, show.grid.labels = 4)
  
  if (nrow(evaluation@D) > 1) {
    for (i in 1:(nrow(evaluation@D) - 1)) {
      values <- as.numeric(evaluation@D[i+1,])[-1]
      radial.plot(rev(values), rp.type = "p", start = pi/2, clockwise = TRUE,
                  labels = labels, 
                  radial.lim = c(0,1), line.col = cols[i+1], lwd = 4,
                  grid.bg = bgGrid, show.grid.labels = 4, add = TRUE)
    }
  }
  
  legend(0.65, 1.2, legend = as.character(evaluation@D[,1]),
         bty = "n", col = cols,
         lwd = 3)
  
}  