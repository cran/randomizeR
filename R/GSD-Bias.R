#' @import mvtnorm
#' @import gsDesign
#' @import insight
#' @import reshape2
#' @import stats
#' @include rarSeq.R
#' @include rarPar.R
NULL
###########################################################


library(mvtnorm)      # Integral calculation
library(gsDesign)     # Adjusted significance levels
library(insight)      # nice look for dataframes
library(reshape2)     # Aggregate data

# Calculate the responses based on randomization sequence and biasing policy
#
# @param seq Vector indicating the randomization sequence, where 0 denotes Group A and 1 denotes Group B
# @param n_A Allocated patients to group A in previous stages
# @param n_B Allocated patients to group B in previous stages
#
# @return A list containing counts for positive responses in Group A and B, negative responses in Group A and B, neutral responses, and final counts of n_A and n_B.

response <- function(seq, n_A, n_B) {
  neutral = pos_A = pos_B = neg_A = neg_B = 0
  for (i in 1: length(seq)) {
    if (n_A == n_B) {
      neutral = neutral +1
    }
    if (seq[i]==0) 			  # Patient was allocated to Group A
    {
      if (n_A < n_B) { 	  	# Patient was expected to be allocated to Group A
        pos_A =pos_A+1
      } else if (n_A > n_B) {	# Patient was expected to be allocated to Group B
        neg_A =neg_A+1
      }
      n_A = n_A+1
    } else if (seq[i]==1) { 	# Patient was allocated to Group B
      if (n_A < n_B) { 		  # Patient was expected to be allocated to Group A
        pos_B =pos_B+1
      } else if (n_A > n_B) {	# Patient was expected to be allocated to Group B
        neg_B =neg_B+1
      }
      n_B = n_B+1
    }
  }
  return (list(pos_A,pos_B,neg_A,neg_B, neutral, n_A, n_B))
}

# Calculate normal parameters based on input conditions.
#
# @param seq A vector indicating the randomization sequence.
# @param nu The allocation bias.
# @param K The number of stages.
# @param n The sample size.
#
# @return A list containing the expected value (EV), covariance matrix (Cov), number of consecutive stages
# that allocated only to one group (zero_stages), and information for each stage (I).

normal_parameters <- function(seq, nu, K, n) {
  # Calculate zero stages
  k = n/K 					        	# Patients in each stage
  zero_stages = 0             # Consecutive stages that allocated only to one group
  subseq_sum = rep(0, K)      # counts the stages where this happens

  for (j in 1:K) {            # in that case we skip the test and move to the next stage
    subseq_sum[j] = sum(seq[0:(j*k)])
    if ((subseq_sum[j] ==0 ) | (subseq_sum[j] ==k*j )) {
      zero_stages = zero_stages + 1
    }
  }

  # Calculate Expected Value
  EV = rep(0, K)
  n_A = n_B = I = Summe = resp = vector("list", K)
  n_A[[1]] = n_B[[1]] = 0
  # Handle zero stages
  if ((zero_stages > 0) & (zero_stages < K)) {        # skip zero stages for EV and Cov calculation
    subseq = seq[1:(k*zero_stages)]                   # instead only calculate good / bad responses
    i= zero_stages
    resp[[i]] = response(subseq, n_A[[1]], n_B[[1]])  # list of good / bad responses in group A and B for stage i
    n_A[i+1] = resp[[i]][[6]]                         # total sample size in group A until stage i
    n_B[i+1] = resp[[i]][[7]]                         # total sample size in group B until stage i
  }

  for (i in (1+zero_stages):K) {                      # for each (nonzero) stage
    subseq = seq[((i-1)*k+1):(i*k)]                   # Allocations in Stage i
    resp[[i]] = response(subseq, n_A[[i]], n_B[[i]])  # list of good / bad responses in group A and B for stage i
    n_A[i+1] = resp[[i]][[6]]                         # total sample size in group A until stage i
    n_B[i+1] = resp[[i]][[7]]                         # total sample size in group B until stage i
    sigma = 0.73
    I[i] =  1 / ( sigma/ n_A[[i+1]]+ sigma /n_B[[i+1]] )   # Information for each stage

    # Calculate sum for expected value
    for (j in (1+zero_stages):i) {                   # Summands have to be recalculated for each stage because of because of division by n_A[[i+1]]
      Summe[[j]] = + nu[j]* (1/n_A[[i+1]] * (resp[[j]][[1]] - resp[[j]][[3]]) - 1/n_B[[i+1]] * (resp[[j]][[2]]- resp[[j]][[4]]))
      if (!(j==(1+zero_stages))) {                # sum up only if we have more than one summand
        Summe[[j]] = Summe[[j]] + Summe[[j-1]]
      }
    }
    EV[[i]] = sqrt(I[[i]]) * (Summe[[i]])
  }
  # Calculate Covariance Matrix
  Cov = matrix(0,K,K)

  for (i in (1+zero_stages):K) {
    for (j in i:K) {
      Cov[i,j] = Cov[j,i ] = sqrt(I[[i]]) / sqrt(I[[j]])
    }
  }
  return(list(EV, Cov, zero_stages, I))
}
# Calculate T1E for a given normal variable and testdesign
#
# @param mu Expected value
# @param cov Covariance
# @param sfu alpha spending function (Pocock, OF, sfLDOF)
# @param k stage number
# @param K total amount of stages
# @param ui update information in each stage (yes/no)
#
# @return A list of T1E probabilities for each stage
typeIerror <- function(mu, cov, sfu, k, K, I = 1:K / K, ui, zero_stages) {

  # If 0 patients have been allocated to one group the test will never reject H0
  if (zero_stages >= k) {
    return(0)
  }
  # For the case of K=1
  if (K==1) {
    set.seed(1)
    return(round(1 - pmvnorm(lower = -1.96, upper = 1.96, mean = mu, sigma = cov, abseps = 1e-16), digits = 5))
  }

  # Information Update
  if (ui == "yes") {
    # Replace NULL entries with 0 in I
    for (i in 1:length(I)) {
      if (is.null(I[[i]])) {
        warning("If there are only allocations to one group the information cannot be updated.")
        I[[i]] = 0.000001 # will be skipped, only set so that the numerical computation can be executed
      }
    }
    info_frac = as.numeric(I[1:K]) / as.numeric(I[K])  # Update information according to allocation
    testdesign = gsDesign(k = K, test.type = 2, sfu = sfu, n.I = info_frac)
  } else {
    testdesign = gsDesign(k = K, test.type = 2, sfu = sfu)  # Do not update information according to allocation
  }

  # Calculate integral
  set.seed(1)
  lower_bound = testdesign$lower$bound[(1 + zero_stages):k]
  upper_bound = testdesign$upper$bound[(1 + zero_stages):k]
  mean_values = mu[(1 + zero_stages):k]
  cov_matrix = cov[(1 + zero_stages):k, (1 + zero_stages):k]
  integral = pmvnorm(algorithm = Miwa(), lower = lower_bound, upper = upper_bound,
                     mean = mean_values, sigma = cov_matrix, abseps = 1e-16)
  alpha = 1-integral
  alpha = round(alpha, digits=5)
  return(alpha)
}


# Input: n = Sample size, K = amount of stages, sfu = design (Pocock, OF, sfLDOF, ...), nu = Selection Bias, ui = update information (yes/no)
# Output: type I error probability
GSD_selection_bias <- function (n, reps, K, randproc, sfu, nu =0, seed = 42, ui="No", rb = 4, mti = 3, p = 2/3)
{
  # error control
  if (!((n/K) %% 1 == 0)) {
    stop("The amount of stages is not divisible by the sample size.")
  }
  # Check if nu is a constant
  if (length(nu)==1) {
    nu <- rep(nu, K)
  } else if (length(nu) != K) {
    stop("Length of nu is not a constant or a vector of length K")
  }
  randobj <- switch(randproc,
                    "CR" = crPar(n, K = 2),
                    "RAR" = rarPar(n, K = 2, groups = c("0", "1")),
                    "BSD" = bsdPar(n, mti = mti, groups = c("0", "1")),
                    "EBC" = ebcPar(n, p, groups = c("0", "1")),
                    "CHEN" = chenPar(n, mti = mti, p = p, groups = c("0", "1")),
                    "PBR" = rpbrPar(n, rb = rb, groups = c("0", "1")),
                    "MP" = mpPar(n, mti = mti, ratio = c(1, 1)),
                    stop("Invalid randproc parameter.")
  )

  if ((identical(sfu, "POC") || identical(sfu, "OF")) && identical(ui, "yes")) {
    # Check if Pocock or OF designs are used with an information update
    warning("Pocock and OF designs can only be used without an information update. ui has been set to no.")
    ui <- "no"
  }

  seq = genSeq(randobj, reps, seed = seed)         # Generates a randomization sequence from randobj. Second parameter for amount of sequences to be created.
  #      seq@M = matrix(c(1,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0), nrow=1, byrow=TRUE)  #   Hardcoded sequence for testing
  typeIerror = matrix(, nrow=K, ncol=reps)        # rows: stage, column: repetition
  #  plotSeq(seq, plotAllSeq = T)                 # plotting of randomization sequences for debuggging
  gonogo = matrix(, nrow=K, ncol=reps)

  if (!(K==1)) {
    testdesign = gsDesign(k=K, test.type = 2 , sfu = sfu, alpha= 0.05, sfupar=0.25)
  }


  for (j in 1:reps) {
    sum = sum(seq@M[j,])

    if (((sum !=0 ) & (sum !=n ))) {
      EV_Cov = normal_parameters(seq@M[j, ], nu, K, n)      # Returns expected value, covariance
      zero_stages = EV_Cov[[3]]
    } else {
      zero_stages = K
    }

    for (k in 1:K) {
      if ((sum ==0 ) | (sum ==n )) {
        typeIerror[k,j] = 0
        stop("Warning: Every allocation has been to the same treatment group and T1E cannot be computed.")
      } else {
        typeIerror[k, j] = typeIerror(EV_Cov[[1]][1:k], EV_Cov[[2]][1:K, 1:K], sfu, k, K,EV_Cov[[4]],  ui, zero_stages)
      }
    }
  }
  alpha_mean = rowMeans(typeIerror, na.rm=TRUE)
  alpha_sums = rowSums(typeIerror, na.rm=TRUE) / (reps)

  return(list(typeIerror[K, ], seq))
}

#' Calculates the Type I error for different randomization sequences from a randomization procedure for a group sequential design
#'
#' @param n total sample size
#' @param reps number of simulations to be conducted
#' @param sfu Group sequential design used (currently available: \code{"Pocock"} - Pocock, \code{"OF"} - O'Brien & Fleming,
#' \code{sfLDPocock} - Lan & DeMets with Pocock like alpha spending function,
#' \code{sfLDOF} - Lan & DeMets with O'Brien & Fleming like alpha spending function)
#' @param K number of stages
#' @param rp the randomization procedure used (currently available: \code{'"CR"'},
#'  \code{'"RAR"'}, \code{'"BSD"'}, \code{'"CHEN"'}, \code{'"PBR"'}, \code{'"MP"'})
#' @param seed Randomization seed
#' @param ui for Lan & DeMets design. Update critical values after each stage according to allocation ratio observed if set to \code{"yes"}.
#' @param rb Block size for randomization procedure PBR.
#' @param mti Maximum tolerated imbalance for randomization procedure BSD and MP.
#' @param p Probability p in favor of the treatment with fewer allocations for EBC and CHEN.
#' @examples
#' #Simulate a group sequential design according to O'Brien and Fleming's design with 24 patients,
#' #10 simulation runs,3 Stages using Random Allocation Rule as a randomization procedure.
#' GSD_allocation(n=24, reps=10, sfu="OF", K=3, rp="RAR")
#' #Simulate a group sequential design according to Lan and deMets design with a Pocock
#' #like alpha spending function with 18 patients, 10 simulation runs,
#' #3 Stages using Permuted Block Randomization with block size 4
#' #as a randomization procedure without updating the critical values after each stage.
#' library(gsDesign)
#' GSD_allocation(n=18, reps=10, sfu=sfLDPocock, K=3, rp="PBR", ui="no", rb=4)
#' @returns A list consisting of a vector of Type I errors for each randomization sequence generated from the randomization procedure and a S4 object of the class of the randomization procedure.
#' @export
#'
GSD_allocation<- function (n, reps, sfu, K, rp, seed = 42, ui="No", rb = 4, mti = 3, p = 2/3)
{
  return(GSD_selection_bias(n, reps, K, randproc = rp, sfu, nu =0, seed, ui = ui, rb, mti, p))
}

#' Calculates the Type I error for a randomization sequence in a group sequential design
#'
#' @param sfu Group sequential design used (currently available: \code{"Pocock"} - Pocock, \code{"OF"} - O'Brien & Fleming,
#' \code{sfLDPocock} - Lan & DeMets with Pocock like alpha spending function,
#' \code{sfLDOF} - Lan & DeMets with O'Brien & Fleming like alpha spending function)
#' @param K number of stages
#' @param seq List of consecutive treatment allocations. 1 for first treatment A, 2 for second treatment.
#' @param ui Only for Lan & DeMets design. Update critical values after each stage according to allocation ratio observed if set to \code{"yes"}.

#' @examples
#' #Simulate a group sequential design according to Pocock's design with 24 patients
#' #and the following consecutive treatment allocation:
#' #A, A, B, A, A, B, A, B, A, B, A, B, A, B, A, B, B, B, A, B, B, A, B, B
#' GSD_allocation_seq(sfu ="Pocock", K=3, seq = c(1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0))
#' #Simulate a group sequential design according to Lan and DeMets with O'Brien & Fleming
#' #like alpha spending with 24 patients and the following consecutive treatment allocation:
#' #A, A, B, A, A, B, A, B, A, B, A, B, A, B, A, B, B, B, A, B, B, A, B, B
#' library(gsDesign)
#' GSD_allocation_seq(sfu =sfLDOF, K=3, seq = c(1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0))
#' @returns A list of type I error probabilities for each stage.
#' @export
#'
GSD_allocation_seq<- function (sfu, K, seq, ui="No")
{
  typeIerror = rep(0, K)
  sum = sum(seq)
  n= length(seq)
  nu = 0
  if (((sum !=0 ) & (sum !=n ))) {
    EV_Cov = normal_parameters(seq, K=K, nu=rep(0, K), n=n)      # Returns EV, COV
    zero_stages = EV_Cov[[3]]
  } else {
    zero_stages = K
  }
  for (k in 1:K) {
    if ((sum ==0 ) | (sum ==n )) {
      typeIerror[k] = 0
      stop("Warning: Every allocation has been to the same treatment group and T1E cannot be computed.")
    } else {
      typeIerror[k] = typeIerror(EV_Cov[[1]][1:k], EV_Cov[[2]][1:K, 1:K], sfu, k, K,EV_Cov[[4]],  ui, zero_stages)
    }
  }
  return(typeIerror)
}

# Calculate the go-no-go probability
#
# @param data_vector vector of T1E probabilites used to calculate the go-nog-go probability
#
# @returns vector of go-no-go probabilities
#
go_nogo = function(data_vector) {
  if (!is.numeric(data_vector)) {
    stop("Input data_vector must be numeric.")
  }

  data_vector = na.omit(data_vector)

  if (length(data_vector) == 0) {
    stop("Input data_vector does not contain any valid values.")
  }
  gonogo = mean(data_vector <= 0.05)
  return(gonogo)
}

#    GSD_allocation(n=24, reps=10, sfu="OF", K=3, rp="RAR")
#    GSD_allocation(n=18, reps=10, sfu=sfLDPocock, K=3, rp="PBR", ui="no", rb=4)
#    GSD_allocation_seq(sfu ="Pocock", K=3, seq = c(1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0))
#    GSD_allocation_seq(sfu =sfLDOF, K=3, seq = c(1,1,0,1,1,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0))
