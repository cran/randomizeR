#' @include doublyF.R
# --------------------------------------------
# Function for the simulated test decision
# --------------------------------------------

# Returns the exact/simulated test decision of Student's t-test
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
  stopifnot(is(randSeq, "randSeq"), randSeq@K == length(endp@mu),
            #is(bias, "chronBias") || is(bias, "selBias") || is(bias, "power"), 
            is(endp, "normEndp"))
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
      if(is(bias, "selBias") && bias@type == "DS")
        stop("Error: Selection bias for K > 2 can be calculated just for convergence strategy.")
      
      doublyF_values(randSeq, bias, endp)$p
    }
    
  }       
}
