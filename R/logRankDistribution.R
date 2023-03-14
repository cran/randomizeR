#' @include getParameters.R

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
#' @title logRankDecSim
#' @importFrom coin logrank_test
#' @importFrom stats rweibull
#' @importFrom stats dweibull
#' @importFrom stats pweibull
#' @importFrom PwrGSD wtdlogrank
#'
logRankDecSim <- function(randSeq, bias, endp){
  stopifnot(is(randSeq, "randSeq"), randSeq@K == 2,
            is(bias, "issue"), is(endp, "endpoint"))
  # calculates the bias matrix
  biasM <- 1 / getExpectation(randSeq, bias, endp)
  followUp <- endp@cenTime - endp@accrualTime

  if (is(endp, "expEndp")){
    # calculates the bias matrix
    biasM <- 1 / getExpectation(randSeq, bias, endp)

    decision <- sapply(1:dim(randSeq@M)[1], function(i) {
      # time matrix
      timeVar <- rexp(randSeq@N, rate = biasM[i,])
      # random censoring, common exponential censoring rate
      randCenVar <- rexp(randSeq@N, rate = endp@cenRate)
      # hard censoring
      endCenVar <- runif(randSeq@N, min = followUp, max = endp@cenTime )
      # observed survival time
      randVar <- pmin(timeVar, randCenVar, endCenVar)
      # censoring status
      status <- (randVar == timeVar)*1

      if (sum(randSeq@M[i,]) == 0 || sum(randSeq@M[i,]) == length(biasM[i, ]) ||
          all(status == 0)) {
        return(FALSE)
      } else {
        # Compute LR test statistic
        sdf <- survdiff(Surv(randVar, status) ~ randSeq@M[i, ])
        # Compute p-values and test decision
        p.value <- 1 - pchisq(sdf$chisq, 1)
        return(as.numeric(p.value <= bias@alpha))
      }
    })
  }

  else if (is(endp, "survEndp")){

    # calculate the distribution parameters
    biasShape <- getDistributionPars(randSeq, bias, endp)$shape
    biasScale <- getDistributionPars(randSeq, bias, endp)$scale

    decision <- sapply(1:dim(randSeq@M)[1], function(i) {
      # time matrix
      timeVar <- rweibull(randSeq@N, shape = biasShape[i,], scale = biasScale[i,])
      # random censoring, common exponential censoring rate
      randCenVar <- rexp(randSeq@N, rate = endp@cenRate)
      # hard censoring
      endCenVar <- runif(randSeq@N, min = followUp, max = endp@cenTime )
      # observed survival time
      randVar <- pmin(timeVar, randCenVar, endCenVar)
      # censoring status
      status <- (randVar == timeVar)*1

      if(endp$maxcombo == TRUE){
        p.value <- maxcombo.pvalue(treated = randSeq@M[i,], time = randVar, status = status)
        return(as.numeric(p.value) <= bias$alpha)
      }

      else {
        if (sum(randSeq@M[i,]) == 0 || sum(randSeq@M[i,]) == randSeq@N ||
            all(status == 0)) {
          return(FALSE)
        }

        else {
          df <- data.frame(randVar, status, treat = randSeq@M[i,])

          if(endp@weights[2] == 0){
            # survival package
            sdf <- survdiff(Surv(randVar, status) ~ treat, data = df, rho = endp@weights[1])
            p.value <- 1 - pchisq(sdf$chisq, 1)
            return(as.numeric(p.value <= bias@alpha))
          }

          else {

            # Coin package
            if(all(status == 1)){
              sdf <- logrank_test(Surv(randVar, status) ~ factor(treat), data = df,
                                  type = "Fleming-Harrington", rho = endp@weights[1], gamma = endp@weights[2])
              p.value <- 1 - pchisq(sdf@statistic@teststatistic^2, 1)
              return(as.numeric(p.value <= bias@alpha))
            }

            # PwrGSD package
            sdf <- wtdlogrank(Surv(randVar, status) ~ treat, data = df, WtFun = "FH", param = endp@weights)
            p.value <- 1 - pchisq(sdf$Z^2, 1)
            return(as.numeric(p.value <= bias@alpha))
          }

        }
      }
    })
  }

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
  if (is(endp, "expEndp")){
    # calculate the bias matrix
    biasM <- 1 / getExpectation(randSeq, bias, endp)

    # function for calculating the rejection probability of each randomization sequence
    rej.prob <- sapply(1:dim(randSeq@M)[1], function(i) {

      # return zero if there were no observations in one treatment group
      if( sum(randSeq@M[i,] == 0) == 0 || sum(randSeq@M[i,] == 1) == 0){
        return(0)
      }

      else {
        # Define approximation functions
        weight <- function(t){ 1 }
        phi <- function(t){ sum( (1-randSeq@M[i,]) * dexp(t, rate = biasM[i,]) ) /
            sum( dexp(t, rate = biasM[i,]) )  }
        pi <- function(t){ sum( (1-randSeq@M[i,]) * (1-pexp(t, rate = biasM[i,])) ) /
            sum( (1-pexp(t, rate = biasM[i,])) )  }
        V <- function(t){ sum( dexp(t, rate = biasM[i,]) ) / randSeq@N *
            ( 1-pexp(t, rate = endp@cenRate) ) *
            ( 1-punif(t, min = followUp, max = endp@cenTime) )  }
        # Define combined functions
        f1 <- function(t){weight(t)*(phi(t)-pi(t))*V(t)}
        f2 <- function(t){weight(t)^2*pi(t)*(1-pi(t))*V(t)}
        # Define integrals
        up <- endp@cenTime
        int1 <- integrate(Vectorize(f1),0,up)$value
        int2 <- integrate(Vectorize(f2),0,up)$value

        # Compute expected value of the normal distribution
        Exp.approx  <- int1/sqrt(1/randSeq@N * int2)
        # Compute rejection probability
        qlow <- qnorm( alpha/2 )
        qup <- qnorm( 1 - alpha/2 )
        rej.prob <- pnorm( qlow, Exp.approx, 1 ) + (1 - pnorm(qup, Exp.approx, 1) )
        rej.prob

        return(rej.prob)
      }
    })
  }

  else if (is(endp, "survEndp")){
    # calculate the distribution parameters
    biasShape <- getDistributionPars(randSeq, bias, endp)$shape
    biasScale <- getDistributionPars(randSeq, bias, endp)$scale

    # function for calculating the rejection probability of each randomization sequence
    rej.prob <- sapply(1:dim(randSeq@M)[1], function(i) {

      # return zero if there were no observations in one treatment group
      if( sum(randSeq@M[i,] == 0) == 0 || sum(randSeq@M[i,] == 1) == 0){
        return(0)
      }

      else {
        # Define approximation functions
        KM <- function(t){ 1  - mean(pweibull(t, shape = biasShape[i,], scale = biasScale[i,])) }
        weight <- function(t){ KM(t)^endp@weights[1] * (1-KM(t))^endp@weights[2] }
        phi <- function(t){ sum( (1-randSeq@M[i,]) *
                                   dweibull(t, shape = biasShape[i,], scale = biasScale[i,]) ) /
            sum(dweibull(t, shape = biasShape[i,], scale = biasScale[i,])) }
        pi <- function(t){ sum((1-randSeq@M[i,]) *
                                 (1-pweibull(t, shape = biasShape[i,], scale = biasScale[i,]))) /
            sum((1-pweibull(t, shape = biasShape[i,], scale = biasScale[i,]))) }
        V <- function(t){ sum( dweibull(t, shape = biasShape[i,], scale = biasScale[i,])) / randSeq@N *
            (1-pexp(t, rate = endp@cenRate)) * (1-punif(t, min = followUp, max = endp@cenTime)) }

        # Define combined functions
        f1 <- function(t){weight(t)*(phi(t)-pi(t))*V(t)}
        f2 <- function(t){weight(t)^2*pi(t)*(1-pi(t))*V(t)}
        # Define integrals
        up <- endp@cenTime
        int1 <- integrate(Vectorize(f1),0,up)$value
        int2 <- integrate(Vectorize(f2),0,up)$value

        # Compute expected value of the normal distribution
        Exp.approx  <- int1/sqrt(1/randSeq@N * int2)
        # Compute rejection probability
        qlow <- qnorm( alpha/2 )
        qup <- qnorm( 1 - alpha/2 )
        rej.prob <- pnorm( qlow, Exp.approx, 1 ) + (1 - pnorm(qup, Exp.approx, 1) )
        rej.prob

        return(rej.prob)
      }
    })
  }
}
