#' @include randSeq.R
NULL

###############################################
# --------------------------------------------#
# Class chenSeq                               #
# --------------------------------------------#
###############################################

# --------------------------------------------
# Function for validity check
# --------------------------------------------


# --------------------------------------------
# Class definition for chenSeq
# --------------------------------------------

# Representation of sequences for the Chen's Design
# 
# @description This set of classes provides functionality of storing randomization
# sequences of Chen's Design along with the parameters 
# representing the design.
# 
# @slot N total number of patients included in the trial.
# @slot p success probability of the biased coin.
# @slot mti The maximum tolerated imbalance during the trial.
# @slot M matrix containing randomization sequences of length \code{N} in its
# rows.
setClass("chenSeq", slots = c(p="numeric", mti = "numeric"), contains = "randSeq")

# --------------------------------------------
# Class definition for rChenSeq
# --------------------------------------------

# Representation of random sequences from the Chen's Design
# 
# @description This set of classes provides functionality of storing random randomization
# sequences of Chen's Design along with the parameters 
# representing the design.
setClass("rChenSeq", contains = c("rRandSeq", "chenSeq"))

# --------------------------------------------
# Methods for chenSeq
# --------------------------------------------
#' @rdname getProbabilities
setMethod("getProb", signature = c(obj = "chenSeq"),
          function(obj) {
            if(obj@K == 2) {
              apply(obj@M, 1, function(x, p, mti){
                rw <- abs(c(0, cumsum(2*x - 1)))
                # hitting zero imbalance
                origin <- sum(rw[-length(rw)] == 0)
                # hitting the mti
                mtiHit <- sum(rw[-length(rw)] == mti)
                # reducing the imbalance object
                together <- sum(rw[-length(rw)] > rw[-1]) 
                # increasing the imbalance
                apart <- sum((rw[-length(rw)] < rw[-1])*(rw[-length(rw)] > 0))
                0.5^origin * p^(together-mtiHit) * (1 - p)^apart 
              }, p = coin(obj), mti = mti(obj))
            } 
            else "Only supported for K=2."
          }  
)

#' @rdname getDesign
setMethod("getDesign", 
          signature(obj = "chenSeq"),
          function(obj) {
            paste("CHEN(", obj@mti, ", ",  round(obj@p, digits = 2), ")", sep = "")
          }
)

