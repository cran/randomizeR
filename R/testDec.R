#' @include doublyF.R
#' @include logRankDistribution.R
# --------------------------------------------
# Function for the simulated test decision
# --------------------------------------------

# Returns the exact/simulated test decision of Student's t-test or the log rank test.
#
# Calculates for every randomization sequence the exact/simulated p.value.
#
# @param randSeq object of the class randSeq.
# @param bias object of the class bias.
# @param endp object of the class endpoint.
#
# @return
# vector of the simulated/exact p.value of a randomization sequence.
testDec <- function(randSeq, bias, endp) {
  stopifnot(is(randSeq, "randSeq"),
            #is(bias, "chronBias") || is(bias, "selBias") || is(bias, "power"),
            is(endp, "normEndp") || is(endp, "expEndp") || is(endp, "survEndp") )

  if (is(endp, "normEndp"))
  {
    stopifnot(randSeq@K == length(endp@mu))
    if(is(bias, "combinedBias") && bias@typeSB == "DS"){
      stop("Error: Selection bias for K > 2 can only be calculated for convergence strategy.")
    }

    if (bias@method == "sim") {
      # calculates the bias matrix
      if(randSeq@K >= 2 && is(bias, "selBias") && bias@type != "DS"){
        R_ <- genSeq(crPar(randSeq@N, randSeq@K))
        biasM <- t(apply(randSeq@M, 1, function(x){
          R_@M <- matrix(x, ncol = N(randSeq))
          makeBiasedExpectation(R_, endp@mu, bias)
        }))
      } else{
        biasM <- getExpectation(randSeq, bias, endp)
      }

      # matrix of the standard deviations
      sdM <- matrix(numeric(0), ncol = dim(randSeq@M)[2], nrow = dim(randSeq@M)[1])
      for(i in 1:randSeq@K){
        sdM[randSeq@M == (i-1)] <- endp@sigma[i]
      }

      sapply(1:dim(randSeq@M)[1], function(i) {
        randVar <- rnorm(length(biasM[i, ]) , mean = biasM[i, ], sd = sdM[i, ] )
        if (sum(randSeq@M[i,]) == 0 || sum(randSeq@M[i,]) == length(biasM[i, ]) ) {
          return(as.numeric(FALSE))
        } else {
          as.numeric(anova(lm(randVar ~ randSeq@M[i,]))[1, "Pr(>F)"] <= bias@alpha)
        }
      })
    } else if (bias@method == "exact") {
      if(randSeq@K == 2){
        if(is(bias, "selBias") && bias@type == "CS2"){
          doublyF_values(randSeq, bias, endp)$p
        } else{
          doublyTValues(randSeq, bias, endp)
        }
      } else{
        if(is(bias, "selBias") && bias@type == "DS"){
          stop("Error: Selection bias for K > 2 can only be calculated for convergence strategy.")
        }
        doublyF_values(randSeq, bias, endp)$p
      }

    }

  }


  else if (is(endp, "expEndp"))
  {
    stopifnot(randSeq@K == 2)
    if(is(bias, "selBias") && bias@type == "CS2"){
      stop("Error: Selection bias for exponential endpoints can only be calculated for convergence and divergence strategy.")
    }
    if(is(bias, "combinedBias") && bias@typeSB == "CS2"){
      stop("Error: Selection bias for exponential endpoints can only be calculated for convergence and divergence strategy.")
    }

    if (bias@method == "sim") {
      logRankDecSim(randSeq, bias, endp)
    }
    else if (bias@method == "exact") {
      message("The rejection probabilities are calculated using an approximation formula.")
      logRankRejectionProb(randSeq, bias, endp)
    }

  }
  else if (is(endp, "survEndp"))
  {
    stopifnot(randSeq@K == 2)
    if(is(bias, "chronBias")){
      stop("Error: Survival endpoints only permit the consideration of allocation bias.")
    }
    if(is(bias, "combinedBias")){
      stop("Error: Survival endpoints only permit the consideration of selection bias.")
    }
    if(!(all(abs(bias@delta) < endp@shape))){
      stop("Error: The absolute value of the allocation bias parameter delta must be
           smaller than the minimum of the shape parameters.")
    }

    if (bias@method == "sim") {
      logRankDecSim(randSeq, bias, endp)
    }
    else if (bias@method == "exact") {
      message("The rejection probabilities are calculated using an approximation formula.")
      logRankRejectionProb(randSeq, bias, endp)
    }

  }
}
