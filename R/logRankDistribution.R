###############################################
# --------------------------------------------#
# Distribution of the log rank statistic      #
# --------------------------------------------#
###############################################


#' Simulation of the test decision of the log rank test
#'
#' Simulates one test decision of the log rank test in the presence of bias. 
#'
#' @param randSeq object of the class randSeq.
#' @param bias object of the class bias.
#' @param endp object of the class endpoint.
#' @keywords internal
#' @return a simulated test decision for each randomization sequence.
logRankDecSim <- function(randSeq, bias, endp){
  stopifnot(is(randSeq, "randSeq"), randSeq@K == 2,
            is(bias, "issue"), is(endp, "endpoint"))  
  # calculates the bias matrix
  biasM <- 1 / getExpectation(randSeq, bias, endp)
  followUp <- endp@cenTime - endp@accrualTime
  
  decision <- sapply(1:dim(randSeq@M)[1], function(i) {
    # time matrix
    timeVar <- rexp(length(biasM[i,]), rate = biasM[i,])
    # random censoring, common exponential censoring rate
    randCenVar <- rexp(length(biasM[i,]), rate = endp@cenRate)
    # hard censoring
    endCenVar <- runif(length(biasM[i,]), min = followUp, max = endp@cenTime )
    # observed survival time
    randVar <- pmin(timeVar, randCenVar, endCenVar)
    # censoring status
    status <- (randVar == timeVar)*1
    if (sum(randSeq@M[i,]) == 0 || sum(randSeq@M[i,]) == length(biasM[i, ])) {
      return(FALSE)
    } else {
      sdf <- survdiff(Surv(randVar, status) ~ randSeq@M[i, ])
      # Compute p-values and test decision
      p.value <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
      return(as.numeric(p.value <= bias@alpha))
    }
  })
  
  decision
}



#' Approximation of the biased type I error probability (resp. power) of the log rank test
#'
#' Computes the biased type I error probability (resp. power) of the log rank test due to shifts in the expectation vectors 
#' in both treatment groups.
#'
#' @param randSeq object of the class randSeq.
#' @param bias object of the class bias.
#' @param endp object of the class endpoint.
#' @keywords internal
#' @return the rejection probability (type I error probability resp. power) of all randomization sequences.
logRankRejectionProb <- function(randSeq, bias, endp) {
  stopifnot(is(randSeq, "randSeq"), randSeq@K == 2,
            is(bias, "issue"), is(endp, "endpoint"))  
  # calculates the bias matrix
  biasM <- 1 / getExpectation(randSeq, bias, endp)
  # variable for the defined alpha quantile
  alpha <- bias@alpha
  followUp <- endp@cenTime - endp@accrualTime
  
  # function for calculating the rejection probability of each randomization sequence
  rej.prob <- sapply(1:dim(randSeq@M)[1], function(i) {
    # Define approximation functions
    phi <- function(t){ sum( (1-randSeq@M[i,]) * dexp(t, rate = biasM[i,]) ) / sum( dexp(t, rate = biasM[i,]) )  }
    pi <- function(t){ sum( (1-randSeq@M[i,]) * (1-pexp(t, rate = biasM[i,])) ) / sum( (1-pexp(t, rate = biasM[i,])) )  }
    V <- function(t){ sum( dexp(t, rate = biasM[i,]) ) / randSeq@N *
        ( 1-pexp(t, rate = endp@cenRate) ) * ( 1-punif(t, min = followUp, max = endp@cenTime) )  }
    # Define combined functions
    f1 <- function(t){(phi(t)-pi(t))*V(t)}
    f2 <- function(t){pi(t)*(1-pi(t))*V(t)}
    # Define integrals
    up <- endp@cenTime
    int1 <- integrate(Vectorize(f1),0,up)$value
    int2 <- integrate(Vectorize(f2),0,up)$value

    # Compute expected value of the normal distribution
    Exp.approx <- int1/sqrt(1/randSeq@N * int2)
    # Compute rejection probability
    qlow <- qnorm( alpha/2 )
    qup <- qnorm( 1 - alpha/2 )
    rej.prob <- pnorm( qlow, Exp.approx, 1 ) + (1 - pnorm(qup, Exp.approx, 1) )
    rej.prob
    
    return(rej.prob)
  })
  
  rej.prob
}
  
