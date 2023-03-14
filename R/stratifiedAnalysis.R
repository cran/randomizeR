#' @include doublyT.R
NULL




# ----------------------------------------------------
# Function for stratification and stratified analysis
# ----------------------------------------------------


#' Creates stratified sequences
#' Compares stratified sequences to their respective non-stratified version under the influence of bias.
#' @param endp object of class \code{endpoint}.
#' @param pr  at least one object of class \code{randPar} or just a list of objects of class \code{randPar}
#'
#' @details
#' Stratified and Non-stratified versions of a randomization sequence behave differently with respect to issues like selection bias, chronological bias or combined bias.
#' The \code{analyse} function creates both versions of a sequence for each of the specified randomization procedures and analyses them in relation to the bias created according to the \code{theta} and \code{eta} values.
#' The first argument should specify the total sample size of patients. The second argument should be one of class \code{normEndp} describing a normally distributed endpoint.
#' The third argument should be the allocation ratio for the different strata. The fourth argument should be the number of strata in the clinical trial. The fifth and sixth arguments should be the selection bias effect eta and the time trend theta.
#' The seventh argument should be a vector of strings representing different randomization procedures. The strings should be given as described by the \code{getDesign} function. Any additional parameters should be given after the design name of the procedure encapsulated in parenthesis.
#'
#' @return
#' The function returns a matrix that summarizes the performance of the randomization procedures. The values for each randomization procedure represent the percentage of sequences that kept the 5% level under the influence of bias.
#' @name analyse
#'
#' @export
NULL


analyse <- function(N, endp, allocRatio, strata, theta, eta, pr){


  if(!is(endp,"normEndp")){
    stop("Stratified Analysis only possible for normal endpoints")
  }


  if(length(allocRatio)!=strata){
    stop("Number of strata should be equal to length of allocation Ratio.")
  }

  if(!is.numeric(strata) || strata < 0 || strata > N){
    stop("Number of strata might not exceed total sample size.")
  }

  if(!is.vector(pr)){
    stop("Input is not a valid vector. Randomization procedures should be entered as a vector.")
  }
  #############################################################################
  ##### Initialize the amount of repetitions, eta, theta, as well as the bias##
  #############################################################################
  rep <- 1000
  eta <- eta*0.64
  theta <- theta/strata*endp@sigma[1]/strata
  biasVector <- list()
  seed <- genSeq(crPar(20))@seed
  for(i in 1:strata){
    if(eta == 0){
      biasVector[[i]] <- chronBias(type = "linT",theta = theta,method = "exact")
    }else if(theta == 0){
      biasVector[[i]] <- selBias(type = "CS",eta=eta,method="exact")
    }else{
      biasVector[[i]] <- combineBias(selBias(type = "CS", eta = eta, method = "exact"), chronBias(type = "linT", theta = theta, method = "exact"))
    }
  }
  bias <- biasVector[[1]]




  ###########################################################################################
  #### Put all of the randomization Procedures in a list and leave only the unique ones######
  ###########################################################################################
  L <- as.list(pr)
  L <- unique(L)
  ############################################
  ####Initialize all vectors for the output###
  ############################################


  R_strat <- list()
  sw<-list()
  sw_star<-list()
  us <- list()
  uus <- list()
  uw_star <-list()
  uw <-list()
  sus <- list()

  ####
  #### Main Loop
  ####
  for(i in 1:length(L)){
    r <- sum(allocRatio)
    oneunit <- as.integer(N/r)
    n<-vector(mode = "list", length = strata)
    sum<-0
    for(v in 1:(strata-1)){
      n[v] <- allocRatio[v]*oneunit
      sum <- sum + n[[v]]
    }
    n[[strata]] <- N - sum

    ############################################
    ######### Initialize process vector#########
    ############################################
    if(L[[i]] == "CR"){
      procs <- list()
      for(i in 1:strata){
        procs[[i]] <- genSeq(crPar(n[[i]]),seed = seed, r = rep)

      }
      whole_proc <- genSeq(crPar(N),r = rep)
    }else if(L[[i]] == "RAR"){
      procs <- list()
      for(i in 1:strata){
        procs[[i]] <- genSeq(rarPar(n[[i]]),seed = seed, r = rep)

      }
      whole_proc <- genSeq(rarPar(N),r = rep)

    }else if(L[[i]] == "HADA"){
      procs <- list()
      for(i in 1:strata){
        procs[[i]] <- genSeq(hadaPar(n[[i]]),seed = seed, r = rep)

      }
      whole_proc <- genSeq(hadaPar(N),r = rep)

    }else if(L[[i]] == paste(c("PBR(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      blockLength <- as.double(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))

      for(i in 1:strata){
        if(n[[i]] %% blockLength != 0){
          stop("Unable to partition strata of size ",paste0(n[[i]]), " into blocks of size ", paste0(blockLength))

        }
        procs[[i]] <- genSeq(pbrPar(bc = rep(blockLength, n[[i]]/blockLength)),seed = seed, r = rep)
      }

      whole_proc <- genSeq(pbrPar(bc = rep(blockLength,N/blockLength)),r = rep)

    }else if(L[[i]] == paste(c("RPBR(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      vals <-  gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]])
      rb <- as.double(gsub("^(.*?),.*", "\\1", vals))
      if(grepl(",", vals, fixed = TRUE)){
        filledBlock <- as.logical(sub('.*,\\s*', '', vals))
        if(filledBlock != TRUE || filledBlock != FALSE){
          stop("Not a valid logical operator.")
        }
      }else{
        filledBlock <- FALSE
      }



      for(i in 1:strata){
        procs[[i]] <- genSeq(rpbrPar(N=n[[i]],rb=rb,filledBlock = filledBlock),seed = seed, r = rep)

      }
      whole_proc <- genSeq(rpbrPar(N=N,rb=rb,filledBlock = filledBlock),r = rep)


    }else if(L[[i]] == "TBD"){
      procs <- list()
      for(i in 1:strata){
        procs[[i]] <- genSeq(tbdPar(bc = n[[i]]),seed = seed, r = rep)

      }
      whole_proc <- genSeq(tbdPar(bc = N),r = rep)

    }else if(L[[i]] == "RTBD" || L[[i]] == paste(c("RTBD(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "") ){

      if(L[[i]] == paste(c("RTBD(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
        filledBlock <- as.logical(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))
        if(filledBlock != TRUE || filledBlock != FALSE){
          stop("Not a valid logical operator.")
        }
      }
      procs <- list()

      for(i in 1:strata){
        procs[[i]] <- genSeq(rtbdPar(N = n[[i]],rb=n[[i]],filledBlock = FALSE),seed = seed, r = rep)

      }
      whole_proc <- genSeq(rtbdPar(N = N,rb=N,filledBlock = filledBlock),r = rep)


    }else if(L[[i]] == paste(c("EBC(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      p <- as.double(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))
      for(i in 1:strata){
        procs[[i]] <- genSeq(ebcPar(N=n[[i]],p=p),seed = seed, r = rep)

      }
      whole_proc <- genSeq(ebcPar(N=N,p=p),r = rep)


    }else if(L[[i]] == paste(c("BSD(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      mti <- as.double(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))
      procs <- list()
      for(i in 1:strata){
        procs[[i]] <- genSeq(bsdPar(N = n[[i]], mti = mti),seed = seed, r = rep)

      }
      whole_proc <- genSeq(bsdPar(N = N, mti = mti),r = rep)

    }else if(L[[i]] == paste(c("MP(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      mti <- as.double(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))
      for(i in 1:strata){
        procs[[i]] <- genSeq(mpPar(N=n[[i]],mti=mti),seed = seed, r = rep)

      }
      whole_proc <- genSeq(mpPar(N=N,mti=mti),r = rep)

    }else if(L[[i]] == paste(c("UD(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      vals <-  gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]])
      ini <- as.double(gsub("^(.*?),.*", "\\1", vals))
      add <- as.double(sub('.*,\\s*', '', vals))
      for(i in 1:strata){
        procs[[i]] <- genSeq(udPar(N=n[[i]],ini=ini,add=add),seed = seed, r = rep)

      }
      whole_proc <- genSeq(udPar(N=N,ini=ini,add=add),r = rep)

    }else if(L[[i]]=="CHEN"){
      procs <- list()
      for(i in 1:strata){
        procs[[i]] <- genSeq(chenPar(N=n[[i]]),seed = seed, r = rep)

      }
      whole_proc <- genSeq(chenPar(N=N),r = rep)

    }else if(L[[i]] == paste(c("GBCD(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      rho <- as.double(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))
      for(i in 1:strata){
        procs[[i]] <- genSeq(gbcdPar(N=n[[i]],rho=rho),seed = seed, r = rep)

      }
      whole_proc <- genSeq(gbcdPar(N=N,rho=rho),r = rep)

    }else if(L[[i]] == paste(c("ABCD(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      a <- as.double(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))
      for(i in 1:strata){
        procs[[i]] <- genSeq(abcdPar(N=n[[i]],a=a),seed = seed, r = rep)

      }
      whole_proc <- genSeq(abcdPar(N=N,a=a),r = rep)


    }else if(L[[i]] == paste(c("BBCD(",gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]),")"), collapse = "")){
      procs <- list()
      a <- as.double(gsub("[\\(\\)]", "", regmatches(L[[i]], gregexpr("\\(.*?\\)", L[[i]]))[[1]]))
      for(i in 1:strata){
        procs[[i]] <- genSeq(bbcdPar(N=n[[i]],a=a),seed = seed, r = rep)

      }
      whole_proc <- genSeq(bbcdPar(N=N,a=a),r = rep)

    }else{
      stop("Not a valid randomization procedure.")
    }
    #####################################
    #### Stratified randomization########
    #####################################
    #w*
    tmpsw <- doublyTValues_new(procs,biasVector,endp = endp,weight=FALSE)
    tmpsw <- sum(tmpsw<=0.05)/length(tmpsw)
    tmpsw <- format(round(tmpsw, 2), nsmall = 2)
    sw_star <- c(sw_star,tmpsw)
    #w = 1
    tmpw <- doublyTValues_new(procs,biasVector,endp = endp,weight=TRUE)
    tmpw <- sum(tmpw<=0.05)/length(tmpw)
    tmpw <- format(round(tmpw, 2), nsmall = 2)
    sw<-c(sw,tmpw)
    #normal t-test
    final <- 0
    values <-c()
    for(o in 1:strata){
      temp <- doublyTValues(procs[[o]],bias,endp = endp)
      values <- c(values, temp)
    }
    tmpsuus <- sum(values<=0.05)/length(values)
    tmpsuus <- format(round(tmpsuus, 2), nsmall = 2)
    sus <- c(sus,tmpsuus)


    #####################################
    #### Unstratified randomization######
    #####################################

    unstr_seq <- list()
    copy_of_N <- c(1,n)
    for(i in 2:length(copy_of_N)){
      unstr_seq[[i-1]] <- whole_proc
      unstr_seq[[i-1]]@N <- n[[i-1]]
      unstr_seq[[i-1]]@M <- whole_proc@M[,((i-2)*copy_of_N[[i-1]]+1):((i-2)*copy_of_N[[i-1]]+copy_of_N[[i]])]
    }




    tmpuw_star <- doublyTValues_new(unstr_seq,biasVector,endp = endp,weight=FALSE)
    tmpuw_star <- sum(tmpuw_star<=0.05)/length(tmpuw_star)
    tmpuw_star <- format(round(tmpuw_star, 2), nsmall = 2)
    uw_star<-c(uw_star,tmpuw_star)

    #w = 1
    tmpuw <- doublyTValues_new(unstr_seq,biasVector,endp = endp,weight=TRUE)
    tmpuw<-sum(tmpuw<=0.05)/length(tmpuw)
    tmpuw <- format(round(tmpuw, 2), nsmall = 2)
    uw<-c(uw,tmpuw)
    #normal t-test
    tmpuus <- doublyTValues(whole_proc,bias,endp=endp)
    tmpuus<-sum(tmpuus<=0.05)/length(tmpuus)
    tmpuus <- format(round(tmpuus, 2), nsmall = 2)
    uus <- c(uus,tmpuus)




  }
  #####################################
  ############# Combine Output ########
  #####################################
  tmp <- c()
  for(d in 1:length(allocRatio)-1){
    tmp <- c(tmp, allocRatio[i]*oneunit,":")

  }
  ar<-paste(allocRatio,sep = "",collapse = ":")
  th <- theta
  et <- eta
  M <- cbind(ar,th,et,L,sw_star,sw,sus,uw_star,uw,uus)
  colnames(M) <- c("[Allocation Ratio]", "[Theta]","[Eta]","[Randomization Procedure]","[w*]","[w=1]","[US]","[w*]","[w]","[US]")
  rownames(M) <- 1:length(L)
  M

}
