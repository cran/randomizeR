#' @include getDesign.R
NULL

###############################################
# --------------------------------------------#
# Class randPar                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------

# Validate randomization generators
#
# @param object object
validateRandPar <- function(object) {
  errors <- character()
  N <- object@N
  K <- object@K
  groups <- object@groups
  ratio <- object@ratio

  if(N <= 1) {
    msg <- paste("N should be greater than one.",
                 sep = "", collapse = "")
    errors <- c(errors, msg)
  }

  if(!(length(N) == 1)) {
    msg <- paste("N = ", N, " should have length 1. Has length ", length(N),
               ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(N == ceiling(N))) {
    msg <- paste("N should be an integer.",
                  sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(length(K) == 1)) {
    msg <- paste("K =", K, " should have length 1. Has length ", length(K),
                 ".", sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(K == ceiling(K))) {
    msg <- paste("K should be an integer.",
                  sep = "", collapse = ",")
    errors <- c(errors, msg)
  }

  if(!(length(groups) == K) && K == ceiling(K)) {
    msg <- paste("Length of groups is ", length(groups), ". Should have length ", K,
                 "." , sep = "")
    errors <- c(errors, msg)
  }

  if(sum(duplicated(groups)) > 0) {
    msg <- paste("Duplicated group names selected, must be unique.",
                 sep = "")
    errors <- c(errors, msg)
  }

  if(!(length(ratio) == K) && K == ceiling(K)) {
    msg <- paste("Length of ratio is ", length(ratio), ". Should have length ", K,
                 "." , sep = "")
    errors <- c(errors, msg)
  }

  if(!(all(ratio >= 1))) {
    msg <- paste("All entries of ratio must be greater than 1.",
                 sep = "")
    errors <- c(errors, msg)
  }

  if(!(all(ratio == ceiling(ratio)))) {
    msg <- paste("All entries of ratio must be integers.",
                 sep = "")
    errors <- c(errors, msg)
  }

  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for randPar
# --------------------------------------------

# Randomization parameters
setClass("randPar",
         slots = c(N = "numeric", K = "numeric" , ratio = "numeric",
           groups = "character"),
         validity = validateRandPar)




#' Settings for randomization procedures
#'
#' Randomization procedures in randomizeR are represented by objects that inherit
#' from \code{randPar}. The representation can then be used in order to
#' generate randomization sequences. In order generate a representation of a
#' randomization procedure, call \code{\link{createParam}} or one of the following
#' functions.
#'
#' @section Supported randomization procedures:
#' \itemize{
#'   \item Complete Randomization (\code{\link{crPar}})
#'   \item Random Allocation Rule (\code{\link{rarPar}})
#'   \item Permuted Block Randomization (\code{\link{pbrPar}})
#'   \item Permuted Block Randomization with random block length (\code{\link{rpbrPar}})
#'   \item Truncated Binomial Design (\code{\link{tbdPar}})
#'   \item Truncated Binomial Design with random block length (\code{\link{rtbdPar}})
#'   \item Efron's Biased Coin Design (\code{\link{ebcPar}})
#'   \item Big Stick Design (\code{\link{bsdPar}})
#'   \item Maximal Procedure (\code{\link{mpPar}})
#'   \item Wei's Urn Design (\code{\link{udPar}})
#'   \item Chen's Design (\code{\link{chenPar}})
#'   \item Generalized Biased Coin Design (\code{\link{gbcdPar}})
#'   \item Accelerated Biased Coin Design (\code{\link{abcdPar}})
#'   \item Bayesian Biased Coin Design (\code{\link{bbcdPar}})
#'   \item Hadamard Randomization (\code{\link{hadaPar}})
#' }
#'
#' @seealso Generate randomization sequences \code{\link{genSeq}}.
#' Calculate the the complete set of randomization sequences of a randomization
#' procedure.
#' \code{\link{getAllSeq}}.
#'
#' @name randPar
NULL


# --------------------------------------------
# Accesssor functions for randPar
# --------------------------------------------

#' Method defining the $ operator for the randPar class
#' @keywords internal
#' @inheritParams overview
setMethod("$", "randPar",
          function(x, name) slot(x, name))

#' Function returning the sample size slot of an S4 object
#'
#' @param obj object inheriting from randPar
#' @returns the sample size slot of an \code{S4} object
#' @export
N <- function(obj) {
  if (.hasSlot(obj, "N")) {
    obj@N
  } else {
    stop("object has no slot named N.")
  }
}

#' Function returning the block slot of an S4 object
#'
#' @param obj object of class pbrPar
#' @returns a vector with the lenghts of each block of a \code{pbrPar} object
#' @export
blocks <- function(obj) {
  if (.hasSlot(obj, "bc")) {
    return(obj@bc)
  } else {
    stop("object has no slot for blocks.")
  }
}

#' Function returning the block slot of an S4 object
#'
#' @param obj object of class pbrPar
#' @returns a vector with the lengths of each random block of a \code{pbrPar} object
#' @export
randBlocks <- function(obj) {
  if (.hasSlot(obj, "rb")) {
    obj@rb
  } else {
    stop("object has no random blocks.")
  }
}

#' Function returning the MTI slot of an S4 object
#'
#' @param obj object of class bsdPar or mpPar
#'
#' @export
mti <- function(obj) {
  if (.hasSlot(obj, "mti")) {
    obj@mti
  } else {
    stop("object has no slot named mti.")
  }
}

#' Function returning the coin slot of an S4 object
#'
#' @param obj object extending class randPar or randSeq
#' @returns The success probability of the biased coin
#' @export
coin <- function(obj) {
  if (.hasSlot(obj, "p")) {
    obj@p
  } else stop("object has no slot named p.")
}

#' Function returning the number of trial arms slot of an S4 object
#'
#' @param obj object of class randPar
#' @returns  The number of trial arms
#' @export
K <- function(obj) {
  if (.hasSlot(obj, "K")) {
    obj@K
  } else {
    stop("object has no slot named K.")
  }
}

#' Function returning the allocation ratio slot of an S4 object
#'
#' @param obj object of class randPar
#' @returns A vector containng the allocation ratio of an \code{S4} object
#' @export
ratio <- function(obj) {
  if (.hasSlot(obj, "ratio")) {
    obj@ratio
  } else {
    stop("object has no slot named ratio.")
  }
}

#' Function returning the method of an S4 object
#'
#' @param obj object inheriting from randPar
#' @returns  The method of an \code{S4} object
#' @export
method <- function(obj) {
  toupper(sub("Par", "", class(obj)[1]))
}

#' Function returning the adjusting parameter rho slot of an S4 object
#'
#' @param obj object of class randPar
#' @return the value of the rho parameter of an \code{S4} object
#' @export
rho <- function(obj) {
  if (.hasSlot(obj, "rho")) {
    obj@rho
  } else {
    stop("object has no slot named rho.")
  }
}

#' Function returning the adjusting parameter a slot of an S4 object
#'
#' @param obj object of class randPar
#' @returns  the value of the adjusting parameter \code{a} of an
#' \code{S4} object
#' @export
a <- function(obj) {
  if (.hasSlot(obj, "a")) {
    obj@a
  } else {
    stop("object has no slot named a.")
  }
}

# --------------------------------------------
# Show function for randPar
# --------------------------------------------

setMethod("show", "randPar", function(object) {
  validObject(object)
  # headline
  cat("\nObject of class \"", class(object)[1], "\"\n\n", sep = "")
  # crop the method from the class name of the randPar object
  cat("design =", getDesign(object), "\n")
  # iterate through all slots of the randPar object
  names <- slotNames(object)
  if (K(object) == 2) names <- names[!(names %in% "K")]
  if (all(ratio(object) == rep(1, K(object)))) {
    names <- names[!(names %in% "ratio")]
  }

  for(name in names) {
    cat(name, "=", slot(object, name),"\n")
  }
  cat("\n")
})


# --------------------------------------------
# Generic functions for randPar
# --------------------------------------------

#' Complete set of randomization sequences
#'
#' Computes all randomization sequences for the given randomization procedure,
#' and stores them in an object along with the parameters belonging to the
#' randomization procedure.
#'
#' @details \code{getAllSeq} is a generic function which dispatches different
#' methods depending on the type of input. The set of sequences of a procedure
#' is computed by enumerating all possible sequences and eliminating those that
#' are not possible in the randomization procedure specified by \code{obj}. The
#' parameters of the randomization procedure are saved along with the sequences
#' to ensure reproducibility of the results.
#'
#' @inheritParams overview
#'
#' @return An object inheriting from \linkS4class{randSeq}, representing the set
#' of randomization sequences for the given parameters.
#' The output consists of the parameters used for the generation of the
#' randomization sequences (see \code{\link{createParam}}) and the matrix \code{M}
#' that stores the randomization sequences in its rows.
#'
#' @seealso \code{\link{createParam}}
#'
#' @examples
#' # all randomization sequences of Efron's Biased Coin Design with p = 0.667 for N = 6
#' myPar <- ebcPar(6, 0.667)
#' getAllSeq(myPar)
#'
#' # all randomization sequences of Big Stick Design with mti = 2 for N = 6
#' myPar <- bsdPar(6, 2)
#' getAllSeq(myPar)
#'
#' # all randomization sequences of Permuted Block Randomization with block sizes 4 and 2
#' myPar <- pbrPar(c(4, 2))
#' getAllSeq(myPar)
#'
#'
#'
#' @name generateAllSequences
NULL

#' @rdname generateAllSequences
#'
#' @export
setGeneric("getAllSeq", function(obj) standardGeneric("getAllSeq"))

#' Generate random sequences
#'
#' Generates randomization sequences from a given randomization procedure.
#'
#' @details
#' \code{genSeq} generates randomization sequences for a randomization
#' procedure as defined by the input parameters.
#' \code{genSeq} has two modes, according to the input.
#' \enumerate{
#'   \item \code{genSeq(obj,r)}: gives \code{r} random sequences from the
#'   design specified by \code{obj}, along with the parameters stored in \code{obj}.
#'   \item \code{genSeq(obj)}: gives one random sequences from the
#'   design specified by \code{obj}, along with the parameters stored in \code{obj}.
#' }
#' The sequences are generated by using the Monte-Carlo sampling technique to sample
#' from the true distribution of the sequences according to the randomization procedure
#' specified by \code{obj}.
#' The parameters of the randomization procedure are saved along with the sequences
#' to ensure reproducibility of the results.
#'
#' @inheritParams overview
#'
#' @return An object inheriting from \linkS4class{randSeq}, representing the \code{r}
#' randomization sequences generated at random for the specified randomization procedure.
#' The output consists of the parameters used for the generation of the
#' randomization sequences (see \code{\link{createParam}}) and the matrix \code{M}
#' that stores the randomization sequences in its \code{r} rows.
#' If \code{r} is missing, one sequence is generated by default.
#'
#' @examples
#' # generate randomization sequences using Complete Randomization for N = 10
#' myPar <- crPar(10)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#'
#' # generate randomization sequences using the Random Allocation Rule for N = 10
#' myPar <- rarPar(10)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#'
#' # generate randomization sequences using the Maximal Procedure with mti = 2 and N = 10
#' myPar <- mpPar(10, 2)
#' genSeq(myPar, 4)
#' genSeq(myPar)
#'
#' @name generateRandomSequences
NULL


#' @rdname generateRandomSequences
#'
#' @export
setGeneric("genSeq", function(obj, r, seed) standardGeneric("genSeq"))


