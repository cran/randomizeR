#' Design Matrix
#' 
#' Calculate Design Matrix from randomization sequence
#' 
#' @param R randomization sequence, object of type randSeq
#' 
#' @return \code{makeDesignMatrix} converts the randomization sequence \code{R} 
#' to its Matrix form. The resulting matrix has \code{K} columns, one for 
#' each treatment group, and \code{N} rows, one for each subject. If a subject \code{i}
#' is randomized to a certain treatment {j}, the entry of \code{(i,j)} of the
#' matrix will be one, and all other entries in this row will be zero.
makeDesignMatrix <- function(R){
  M <- matrix(0, nrow = N(R), ncol = K(R))
  index <- lapply(0:(K(R)-1), function(j) which(as.vector(R@M) == j) )
  for(j in 1:K(R)) M[index[[j]], j]<- 1
  M
}


#' Biasing Policy for a Group of Favoured Treatments 
#' 
#' Calculate vector with the selection bias for each patient
#' 
#' @param R randomization sequence, a integer vector with entries 1, \dots, K, of length N
#' @param K number of treatment groups, a single integer value
#' @param pref preferred groups for the guessing
#' 
#' @return vector with the selection bias for each patient
getbiasCS1 <- function(R, K, pref){
  Gsize <- sapply(seq_len(K)-1, function(x) c(0,cumsum(R[-length(R)] == x)))
  as.vector(diff(apply(Gsize, 1, function(x) c(all(x[pref] > x[-pref]), all(x[pref] < x[-pref])))))
}

#' Biasing Policy for Avoiding the Placebo Treatment
#' 
#' Calculate vector with the selection bias for each patient
#' 
#' @param R randomization sequence, a integer vector with entries 1, \dots, K, of length N
#' @param K number of treatment groups, a single integer value
#' @param avoid avoided groups for the guessing
#' 
#' @return vector with the selection bias for each patient
getbiasCS2 <- function(R, K, avoid){
  Gsize <- sapply(seq_len(K)-1, function(x) c(0,cumsum(R[-length(R)] == x)))
  as.vector(diff(apply(Gsize, 1, function(x) c(any(x[-avoid] < x[avoid]), all(x[-avoid] > x[avoid])))))
}

#' Calculate Expectation vector 
#' 
#' @param R randomization sequence, an object of type randSeq
#' @param mu vector of length K containing the expectation values of groups 1, \dots, K
#' @param bias selection or chronological bias object, containing eta(or theta) and alpha
#' 
#' @return vector of length N with the biased expectation for each patient.
#' 
#' @export
makeBiasedExpectation <-  function(R, mu, bias){
  eta <- bias@eta
  means <- as.vector(makeDesignMatrix(R) %*% mu)
  if (is(bias, "selBias") && bias@type == "CS") {
    return(means + eta*getbiasCS1(R@M, K(R), 1))
  } else if (is(bias, "selBias") && bias@type == "CS2") {
    return(means + eta*getbiasCS2(R@M, K(R), 1)) 
  } else {
    warning("Please specify valid Biasing Policy")
  }
}


#' Calculate hat matrix
#' 
#' @param M Design Matrix representing the randomization sequence
#' 
#' @return Hat matrix, i.e. \code{M * (M^T M)^{-1} M^T} 
hatMatrix <- function(M){
  A <- t(M)%*%M
  A_inv <- solve(A)
  M%*% A_inv %*% t(M)
}

#' Calculate first non centrality parameter of the doubly non central F-distribution
#' 
#' @param H Hat Matrix
#' @param EY (Biased) expectation of the responses
#' 
#' @return First non centrality parameter, a single numeric value
lambda1 <- function(H, EY){
  ovMean <- matrix(rep(1/length(EY), length(EY)^2), ncol=length(EY))
  t(EY) %*% (H - ovMean) %*% EY
}

#' Calculate second non centrality parameter of the doubly non central F-distribution
#' 
#' @param H Hat Matrix
#' @param EY (Biased) expectation of the responses
#' 
#' @return Second non centrality parameter, a single numeric value
lambda2 <- function(H, EY){
  t(EY) %*% (diag(length(EY)) - H) %*% EY
}

#' Distribution function of the non central F-distribution
#' 
#' @param x quantile, a single numeric value
#' @param df1 first degree of freedom, a single integer value
#' @param df2 second degree of freedom, a single integer value
#' @param lambda1 first non centrality parameter, a single numeric value
#' @param lambda2 second non centrality parameter, a single numeric value
#' @param acc accuracy; last index of the approximation of the infinite sum
#' @param ex exactness; break early, if the summands are smaller than \code{10^(-ex)}
#' 
#' @return Probability of observing a value larger than x
doublyF_opt <- function(x, df1, df2, lambda1, lambda2, acc=100, ex = 5){
  c <- exp(-(lambda1 +lambda2)/2)
  f <- x * df1/df2
  s <- 0
  cj <- ck <- 1
  for(j in 0:acc){
    cj <- (lambda1/2)^j/factorial(j)
    for(k in 0:acc){
      ck <- (lambda2/2)^k/factorial(k)
      ibet <- pbeta(f/(1+f), df1/2 + j, df2/2 +k)
      s <- s + cj * ck * ibet
      if (ck < 10^(-ex)) break
    }
    if(cj < 10^(-ex)) break
  }
  s*c
}

#' Rejection probability for one sequence in the presence of selection bias
#' 
#' @param R object of type randSeq representing the randomization procedure
#' @param bias selection or chronological bias object, containing eta(or theta) and alpha
#' @param endp endpoint object, containing mu and sigma
#' 
#' @return data frame containing the non centrality parameters, the degrees of
#' freedom, the \code{1-alpha} quantile corresponding to the central f distribution and
#' a numeric value for probability of false rejection of the null 
#' hypothesis of no difference in the presence of selection bias.
doublyF_value <- function(R, bias, endp){
  alpha <- bias@alpha
  
  if(is(bias, "selBias")){
    EY <- makeBiasedExpectation(R, endp@mu, bias) 
  } else{
    EY <- as.vector(getExpectation(R, bias, endp))
  }
  
  M <- makeDesignMatrix(R)
  H <- hatMatrix(M)
  ovMean <- matrix(rep(1/N(R), N(R)^2), ncol=N(R)) # overall mean
  l1 <- lambda1(H, EY)
  l2 <- lambda2(H, EY)
  df1 <- K(R)-1
  df2 <- N(R)-K(R)
  
  x <- qf(1-alpha, df1, df2) # quantile of the central f distribution
  
  data.frame(
    x = x,
    df1 = df1,
    df2 = df2,
    lambda1 = l1,
    lambda2 = l2,
    p = 1 - doublyF_opt(x = x, df1 = df1, df2 = df2, lambda1 = l1, lambda2 = l2)
  )
}


#' Check function for occurance of all treatment groups in the sequence
#'
#' checks wheather each group has its value comming up at least once in the sequence
#' 
#' @param seq randomization sequence as inverted matrix
#' @param K number of treatment arms
#' 
#' @return TRUE if all groups represented, FAlSE otherwise
hasAllGroups <- function(seq, K){
  for(i in 0:(K-1)){
    if(!(i %in% seq)){
      return(FALSE)
    }
  }
  return(TRUE)
}



#' Rejection probability in case of selection bias in multi-arm trials
#'
#' calculates the non-centrality parameters of the F-distribution under third
#' order selection bias. 
#' @param randSeq the object containing the randomization sequences
#' @param bias selection bias object, containing eta and alpha
#' @param endp endpoint object, containing mu and sigma
#' 
#' @return data frame with the sequences and their corresponding ncps and rejection
#' probabilities.
#'
#' @export
doublyF_values <- function(randSeq, bias, endp){
  seq <- randSeq@M

  R <- genSeq(crPar(N=N(randSeq), K=K(randSeq)))
  prob <- apply(seq, 1, function(x){
    R@M <- matrix(x, ncol = N(R))
    if(!hasAllGroups(R@M, R@K)){
      data.frame(x = 0, df1 = 0, df2 = 0, lambda1 = 0, lambda2 = 0, p = 0)
    } else {
      doublyF_value(R, bias, endp)
    }
  })
  l1 <- sapply(prob, function(item) item$lambda1)
  l2 <- sapply(prob, function(item) item$lambda2)
  p <- sapply(prob, function(item) item$p)
  Sequences <- apply(getRandList(randSeq), 1, paste, collapse = " ")
  data.frame(Sequences, l1, l2, p)
}