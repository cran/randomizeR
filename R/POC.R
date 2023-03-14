#' @import pracma
#' @import cubature
#' @include rarPar.R
NULL

######## Group Sequential Design Pocock using RAR ##############
POC_2_RAR_A1 <- function(n,norep) {
  #############            Pocock, 2 Stages, Random Allocation Rule, applying the randomization procedure to all patients







  #parameters
  K = 2
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))



  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    grpoc = 2.178
    for (i in 1:norep)
    {
      randseq = rarPar(n , K=2 , groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M[1:(n/K)]

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Expected Value Z_1  ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



      ################################### Type 1 Error 1st stage  ##################

      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -2.178, 2.178)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1



      ######################### Expected Value Z_2 ##########################
      no_T_2 = sum(gensequence@M)
      no_C_2 = n - no_T_2

      nu.stage2 = calc.nu(gensequence@M)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

      ######################## Covariance #######################
      ### Covariance Z_1 and Z_2

      cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



      ######################### Type 1 Error 2nd stage ##################################

      f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
      wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2
    }
    Ergebnismatrix_1[m,1] = mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)


  }






  cat("Pocock 2 Stages using RAR - applying the randomization procedure all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = round(Ergebnismatrix_1, digits = 3)
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "p<=0.05", "p<=0.05"))
  print(Ergebnis)







}

POC_3_RAR_A1 <- function(n,norep) {
  #############            Pocock, 3 Stages, Random Allocation Rule, applying the randomization procedure to all patients







  #parameters
  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)




  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))



  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    for (i in 1:norep)
    {
      randseq = rarPar(n , K=2 , groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M[1:(n/K)]

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Expected Value Z_1  ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




      ######################### Expected Value Z_2 ##########################
      Teilsequenz2 = gensequence$M[1:((n/K)*2)]
      no_T_2 = sum(Teilsequenz2)
      no_C_2 = 2*n/K - no_T_2

      nu.stage2 = calc.nu(Teilsequenz2)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





      ############################Expected Value Z_3  ##############################
      no_T_3 = sum(gensequence@M)
      no_C_3 = n - no_T_3

      nu.stage3 = calc.nu(gensequence@M)
      vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
      factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
      mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


      ########################### Covariance of Z1 and Z2 ###########################################
      cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Covariance of Z1 and Z3 ###########################################
      cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Covariance of Z2 and Z3 ###########################################
      cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

      #### Covariance Matrix  #############
      Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
      Sigma_inv = solve(Sigma)
      det_Sigma = det(Sigma)






      #### Type 1 Errors #####

      ######################### Type 1 Error 1st stage ##################################

      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -2.289, 2.289)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1


      ######################### Type 1 Error 2nd stage ##################################

      f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
      wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2

      ######################### Type 1 Error 3rd stage ##################################

      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
      wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
      wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
      alpha_3 = alpha_2 + wert_4 + wert_5

      ergebnis.vector.alpha.3[i] = alpha_3




    }
    Ergebnismatrix_1[m,1] = mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)

  }





  cat("Pocock 3 Stages using RAR - applying the randomization procedure all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = round(Ergebnismatrix_1, digits = 3)

  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05"))
  print(Ergebnis)






}

POC_4_RAR_A1 <- function(n,norep) {
  #############            Pocock, 4 Stages, Random Allocation Rule, applying the randomization procedure to all patients







  #parameters
  K = 4
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)



  Ergebnismatrix_1= matrix(rep(NA,length(eta_vector)*(K+4)), nrow = length(eta_vector))


  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    for (i in 1:norep)
    {
      randseq = rarPar(n , K=2 , groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M[1:(n/K)]

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Expected Value Z_1  ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




      ######################### Expected Value Z_2 ##########################
      Teilsequenz2 = gensequence$M[1:((n/K)*2)]
      no_T_2 = sum(Teilsequenz2)
      no_C_2 = 2*n/K - no_T_2

      nu.stage2 = calc.nu(Teilsequenz2)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

      ######################### Expected Value Z_3 ##########################
      Teilsequenz3 = gensequence$M[1:((n/K)*3)]
      no_T_3 = sum(Teilsequenz3)
      no_C_3 = 3*n/K - no_T_3

      nu.stage3 = calc.nu(Teilsequenz3)
      vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
      factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
      mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





      ############################ Expected Value Z_4 ##############################
      no_T_4 = sum(gensequence@M)
      no_C_4 = n - no_T_4

      nu.stage4 = calc.nu(gensequence@M)
      vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
      factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
      mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta


      ########################### Covariance of Z1 and Z2  ###########################################
      cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

      ###########################  Covariance of Z1 and Z3 ###########################################
      cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

      ###########################  Covariance of Z1 and Z4 ###########################################
      cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



      ########################### Covariance of Z2 and Z3 ###########################################
      cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
      ###########################  Covariance of Z2 and Z4 ###########################################
      cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

      ###########################  Covariance of Z3 and Z4 ###########################################
      cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






      #### Covariancematrix #############
      Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
      Sigma_inv = solve(Sigma)
      det_Sigma = det(Sigma)
      Sigma_3 = Sigma[1:3, 1:3]
      Sigma_inv_3 = solve(Sigma_3)
      det_Sigma_3 = det(Sigma_3)




      gpoc = 2.361

      #### Type 1 Errors  #####

      ######################### Type 1 Error 1st stage  ##################################

      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -gpoc, gpoc)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1


      #########################  Type 1 Error 2nd stage  ##################################

      f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
      wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2

      #########################  Type 1 Error 3rd stage  ##################################

      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
      wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
      wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
      alpha_3 = alpha_2 + wert_4 + wert_5

      ergebnis.vector.alpha.3[i] = alpha_3

      #########################  Type 1 Error 4th stage  ##################################

      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
      wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
      wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
      alpha_4 = alpha_3 + wert_6 + wert_7

      ergebnis.vector.alpha.4[i] = alpha_4




    }

    Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
    Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)

  }





  cat("Pocock 4 Stages using RAR - applying the randomization procedure all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3))
  Ergebnis = round(Ergebnismatrix_1, digits = 3)

  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "3","4", "p<=0.05", "p<=0.05", "p<=0.05","p<=0.05"))
  print(Ergebnis)






}

POC_2_RAR_A2 <- function(n,norep) {
  #############            Pocock, 2 Stages, Random Allocation Rule, applying the randomization procedure to each stage







  #parameters
  K = 2
  alpha = 0.05
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))



  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    grpoc = 2.178
    for (i in 1:norep)
    {
      randseq = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Erwartungswert für Teststatistik ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



      ################################### Fehlerwahrscheinlichkeit in der 1.Stufe ##################
      #
      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -2.178, 2.178)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1



      ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
      randseq_2 = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence_2 = genSeq(randseq_2,1)
      Teilsequenz2 = c(gensequence$M, gensequence_2$M)
      no_T_2 = sum(Teilsequenz2)
      no_C_2 = n - no_T_2


      nu.stage2 = calc.nu(Teilsequenz2)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

      ######################## Kovarianz #######################
      ### Kovarianz von Z1 und Z2

      cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



      ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################
      #
      f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
      wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2
    }
    Ergebnismatrix_1[m,1] = mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)

  }






  cat("Pocock 2 Stages using RAR - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = round(Ergebnismatrix_1, digits = 3)
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "p<=0.05", "p<=0.05"))
  print(Ergebnis)







}

POC_3_RAR_A2 <- function(n,norep) {
  #############            Pocock, 3 Stages, Random Allocation Rule, applying the randomization procedure to each stage






  #parameters
  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)




  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))



  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    for (i in 1:norep)
    {
      randseq = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1



      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




      ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
      randseq_2 = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence_2 = genSeq(randseq_2,1)
      Teilsequenz2 = c(gensequence$M, gensequence_2$M)
      no_T_2 = sum(Teilsequenz2)
      no_C_2 = 2*n/K - no_T_2

      nu.stage2 = calc.nu(Teilsequenz2)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





      ############################ Erwartungswert für Teststatistik in der 3. Stufe  ##############################
      randseq_3 = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence_3 = genSeq(randseq_3,1)
      Teilsequenz3 = c(Teilsequenz2, gensequence_3$M)
      no_T_3 = sum(Teilsequenz3)
      no_C_3 = n - no_T_3


      nu.stage3 = calc.nu(Teilsequenz3)
      vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
      factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
      mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


      ########################### Kovarianzen von Z1 und Z2 ###########################################
      cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Kovarianzen von Z1 und Z3 ###########################################
      cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Kovarianzen von Z2 und Z3 ###########################################
      cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

      #### Kovarianzmatrix #############
      Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
      Sigma_inv = solve(Sigma)
      det_Sigma = det(Sigma)






      #### Fehlerwahrscheinlichkeiten #####

      ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################
      #
      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -2.289, 2.289)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1


      ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################
      #
      f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
      wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2

      ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################
      #
      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
      wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
      wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
      alpha_3 = alpha_2 + wert_4 + wert_5

      ergebnis.vector.alpha.3[i] = alpha_3




    }
    Ergebnismatrix_1[m,1] = mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)


  }




  cat("Pocock 3 Stages using RAR - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = round(Ergebnismatrix_1, digits = 3)

  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05"))
  print(Ergebnis)






}

POC_4_RAR_A2 <- function(n,norep) {
  #############            Pocock, 4 Stages, Random Allocation Rule, applying the randomization procedure to each stage






  #parameters
  K = 4
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))



  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    for (i in 1:norep)
    {
      randseq = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {
          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




      ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
      randseq_2 = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence_2 = genSeq(randseq_2,1)
      Teilsequenz2 = c(gensequence$M, gensequence_2$M)
      no_T_2 = sum(Teilsequenz2)
      no_C_2 = 2*n/K - no_T_2

      nu.stage2 = calc.nu(Teilsequenz2)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

      ######################### Erwartungswert für Teststatistik in der 3. Stufe ##########################
      randseq_3 = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence_3 = genSeq(randseq_2,1)
      Teilsequenz3 = c(Teilsequenz2, gensequence_3$M)
      no_T_3 = sum(Teilsequenz3)
      no_C_3 = 3*n/K - no_T_3

      nu.stage3 = calc.nu(Teilsequenz3)
      vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
      factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
      mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





      ############################ Erwartungswert für Teststatistik in der 4. Stufe  ##############################
      randseq_4 = rarPar(n/K , K=2, groups = c("0","1"))
      gensequence_4 = genSeq(randseq_4,1)
      Teilsequenz4 = c(Teilsequenz3, gensequence_4$M)
      no_T_4 = sum(Teilsequenz4)
      no_C_4 = n - no_T_4


      nu.stage4 = calc.nu(Teilsequenz4)
      vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
      factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
      mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta


      ########################### Kovarianzen von Z1 und Z2 ###########################################
      cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Kovarianzen von Z1 und Z3 ###########################################
      cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Kovarianzen von Z1 und Z4 ###########################################
      cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



      ########################### Kovarianzen von Z2 und Z3 ###########################################
      cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
      ########################### Kovarianzen von Z2 und Z4 ###########################################
      cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Kovarianzen von Z3 und Z4 ###########################################
      cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






      #### Kovarianzmatrix #############
      Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
      Sigma_inv = solve(Sigma)
      det_Sigma = det(Sigma)
      Sigma_3 = Sigma[1:3, 1:3]
      Sigma_inv_3 = solve(Sigma_3)
      det_Sigma_3 = det(Sigma_3)




      gpoc = 2.361

      #### Fehlerwahrscheinlichkeiten #####

      ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################
      #
      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -gpoc, gpoc)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1


      ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################
      #
      f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
      wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2

      ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################
      #
      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
      wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
      wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
      alpha_3 = alpha_2 + wert_4 + wert_5

      ergebnis.vector.alpha.3[i] = alpha_3

      ######################### Fehlerwahrscheinlichkeit der 4. Stufe ##################################
      #
      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
      wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
      wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
      alpha_4 = alpha_3 + wert_6 + wert_7

      ergebnis.vector.alpha.4[i] = alpha_4

    }

    Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
    Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)

  }




  cat("Pocock 4 Stages using RAR - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3))
  Ergebnis = round(Ergebnismatrix_1, digits = 3)

  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "3","4", "p<=0.05", "p<=0.05", "p<=0.05","p<=0.05"))
  print(Ergebnis)






}



######## Group Sequential Design Pocock using PBD(2) ##############

POC_2_PBD2 <- function(n,norep) {
  ################### Pocock, 2 Stages, PBD(2), applying randomization procedure to all patients







  #parameters
  K = 2
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))



  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    grpoc = 2.178
    for (i in 1:norep)
    {
      randseq = ebcPar(n , p = 1 , groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M[1:(n/K)]

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Expected Value of Z_1 ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



      ################################### Type 1 Error 1st Stage ##################

      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -2.178, 2.178)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1



      ######################### Expected Value of Z_2 ##########################
      no_T_2 = sum(gensequence@M)
      no_C_2 = n - no_T_2

      nu.stage2 = calc.nu(gensequence@M)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

      ########################Covariance  #######################


      cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



      ######################### Type 1 Error 2nd Stage ##################################

      f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
      wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2
    }
    Ergebnismatrix_1[m,1] = mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)


  }





  cat("Pocock 2 Stages using PBD(2) - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = round(Ergebnismatrix_1, digits = 3)
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "p<=0.05", "p<=0.05"))
  print(Ergebnis)







}


POC_3_PBD2 <- function(n,norep) {

  ################### Pocock, 3 Stages, PBD(2), applying randomization procedure to all patients







  #parameters
  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)




  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))



  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    for (i in 1:norep)
    {
      randseq = ebcPar(n , p=1 , groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M[1:(n/K)]

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Expected Value of Z_1  ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




      #########################Expected Value of Z_2 ##########################
      Teilsequenz2 = gensequence$M[1:((n/K)*2)]
      no_T_2 = sum(Teilsequenz2)
      no_C_2 = 2*n/K - no_T_2

      nu.stage2 = calc.nu(Teilsequenz2)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





      ############################ Expected Value of Z_3  ##############################
      no_T_3 = sum(gensequence@M)
      no_C_3 = n - no_T_3

      nu.stage3 = calc.nu(gensequence@M)
      vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
      factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
      mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


      ########################### Covariance of Z_1 and Z_2  ###########################################
      cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Covariance of Z_1 and Z_3 ###########################################
      cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Covariance of Z_2 and Z_3 ###########################################
      cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

      #### Covariance matrix #############
      Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
      Sigma_inv = solve(Sigma)
      det_Sigma = det(Sigma)






      #### Type 1 Errors #####

      ######################### Type 1 Error 1st Stage  ##################################

      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -2.289, 2.289)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1


      ######################### Type 1 Error 2nd Stage ##################################

      f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
      wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2

      ######################### Type 1 Error 3rd Stage ##################################

      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
      wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
      wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
      alpha_3 = alpha_2 + wert_4 + wert_5

      ergebnis.vector.alpha.3[i] = alpha_3




    }
    Ergebnismatrix_1[m,1] = mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)


  }




  cat("Pocock 3 Stages using PBD(2) - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = round(Ergebnismatrix_1, digits = 3)

  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05"))

  print(Ergebnis)








}


POC_4_PBD2 <- function(n,norep) {

  ################### Pocock, 4 Stages, PBD(2), applying randomization procedure to all patients







  #parameters
  K = 4
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))


  for (m in 1:length(eta_vector))
  {


    eta = eta_vector[m]
    set.seed(NULL)

    for (i in 1:norep)
    {
      randseq = ebcPar(n , p = 1, groups = c("0","1"))
      gensequence = genSeq(randseq,1)
      subsequence_1 = gensequence$M[1:(n/K)]

      no_T_1 = sum(subsequence_1)
      no_C_1 = n/K-no_T_1




      calc.nu = function(seq)
      {
        nu = rep(0,6)
        noT = 0
        noC = 0
        for (i in 1: length(seq))
        {

          if ( (noT == noC) & (seq[i] == 0))
          {noC = noC +1
          nu[2] = nu[2]+1
          } else if ( (noT > noC) & (seq[i] == 0) )
          {noC = noC +1
          nu[1] = nu[1]+1
          } else if ( (noT < noC) & (seq[i] == 0))
          {noC = noC +1
          nu[3] = nu[3]+1
          }else if ( (noT == noC) & (seq[i] == 1))
          {noT = noT +1
          nu[5] = nu[5]+1
          }else if ( (noT > noC) & (seq[i] == 1))
          {noT = noT +1
          nu[4] = nu[4]+1
          }else if ( (noT < noC) & (seq[i] == 1))
          {noT = noT +1
          nu[6] = nu[6]+1
          }


        }
        return(nu)

      }


      nu.stage1 = calc.nu(subsequence_1)
      nu.stage1


      ############################ Expected Value of Z_1  ##############################
      vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
      factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
      mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




      ######################### Expected Value of Z_2 ##########################
      Teilsequenz2 = gensequence$M[1:((n/K)*2)]
      no_T_2 = sum(Teilsequenz2)
      no_C_2 = 2*n/K - no_T_2

      nu.stage2 = calc.nu(Teilsequenz2)
      vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
      factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
      mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

      ######################### Expected Value of Z_3 ##########################
      Teilsequenz3 = gensequence$M[1:((n/K)*3)]
      no_T_3 = sum(Teilsequenz3)
      no_C_3 = 3*n/K - no_T_3

      nu.stage3 = calc.nu(Teilsequenz3)
      vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
      factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
      mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





      ############################ Expected Value of Z_4  ##############################
      no_T_4 = sum(gensequence@M)
      no_C_4 = n - no_T_4

      nu.stage4 = calc.nu(gensequence@M)
      vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
      factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
      mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta


      ########################### Covariance of Z_1 and Z_2  ###########################################
      cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
      korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Covariance of Z_1 and Z_3  ###########################################
      cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Covariance of Z_1 and Z_4  ###########################################
      cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



      ########################### Covariance of Z_2 and Z_3  ###########################################
      cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
      korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
      ########################### Covariance of Z_2 and Z_4  ###########################################
      cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

      ########################### Covariance of Z_3 and Z_4  ###########################################
      cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
      korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






      #### Covariance matrix #############
      Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
      Sigma_inv = solve(Sigma)
      det_Sigma = det(Sigma)
      Sigma_3 = Sigma[1:3, 1:3]
      Sigma_inv_3 = solve(Sigma_3)
      det_Sigma_3 = det(Sigma_3)




      gpoc = 2.361

      ####  Type 1 Errors #####

      ######################### Type 1 Error 1st Stage ##################################

      ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
      wert_1 = integrate(ftest_1, -gpoc, gpoc)
      alpha_1=  1-wert_1$value
      ergebnis.vector.alpha.1[i] = alpha_1


      ######################### Type 1 Error 2nd Stage ##################################

      f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
      wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
      wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
      alpha_2 = alpha_1 + wert_2 + wert_3

      ergebnis.vector.alpha.2[i] = alpha_2

      ######################### Type 1 Error 3rd Stage ##################################

      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
      wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
      wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
      alpha_3 = alpha_2 + wert_4 + wert_5

      ergebnis.vector.alpha.3[i] = alpha_3

      ######################### Type 1 Error 4th Stage ##################################

      mu = c(mu_1, mu_2, mu_3)
      f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
      wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
      wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
      alpha_4 = alpha_3 + wert_6 + wert_7

      ergebnis.vector.alpha.4[i] = alpha_4




    }

    Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
    Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
    Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
    Ergebnismatrix_1[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
    Ergebnismatrix_1[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)

  }




  cat("Pocock 4 Stages using PBD(2) - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = round(Ergebnismatrix_1, digits = 3)
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05" ))
  print(Ergebnis)


}


######### Group Sequential Design Pocock using BSD(b) ################

POC_2_BSD_A1 <- function(n,norep) {
  #############         Pocock, 2 Stages, Big Stick Design, applying the randomization procedure to all patients






  #parameters
  K = 2
  alpha = 0.5
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3, 5, 10)

  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))
  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      grpoc = 2.178
      for (i in 1:norep)
      {
        randseq = bsdPar(n , mti = p , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



        ################################### Fehlerwahrscheinlichkeit in der 1.Stufe ##################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.178, 2.178)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1



        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        no_T_2 = sum(gensequence@M)
        no_C_2 = n - no_T_2

        nu.stage2 = calc.nu(gensequence@M)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

        ######################## Kovarianz #######################
        ### Kovarianz von Z1 und Z2

        cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
        wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2
      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }
    }
  }


  cat("Pocock 2 Stages using BSD(b) b=3/5/10 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("b=3 1", "2", "p<=0.05", "p<=0.05", "b=5 1", "2", "p<=0.05", "p<=0.05", "b=10 1", "2", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}


POC_3_BSD_A1 <- function(n,norep) {
  #############         Pocock, 3 Stages, Big Stick Design, applying the randomization procedure to all patients






  #parameters

  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3,5,10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))


  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = bsdPar(n , mti = p , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        Teilsequenz2 = gensequence$M[1:((n/K)*2)]
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = 2*n/K - no_T_2

        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





        ############################ Erwartungswert für Teststatistik in der 3. Stufe  ##############################
        no_T_3 = sum(gensequence@M)
        no_C_3 = n - no_T_3

        nu.stage3 = calc.nu(gensequence@M)
        vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
        factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
        mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


        ########################### Kovarianzen von Z1 und Z2 ###########################################
        cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z1 und Z3 ###########################################
        cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z2 und Z3 ###########################################
        cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

        #### Kovarianzmatrix #############
        Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
        Sigma_inv = solve(Sigma)
        det_Sigma = det(Sigma)






        #### Fehlerwahrscheinlichkeiten #####

        ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.289, 2.289)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1


        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
        wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2

        ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
        wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
        wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
        alpha_3 = alpha_2 + wert_4 + wert_5

        ergebnis.vector.alpha.3[i] = alpha_3




      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }
    }
  }

  cat("Pocock 3 Stages using BSD(b) b=3/5/10 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  Ergebnis
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("b=3 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05",
                               "b= 51", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05", "b=10 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}

POC_4_BSD_A1 <- function(n,norep) {
  #############         Pocock, 4 Stages, Big Stick Design, applying the randomization procedure to all patients






  #parameters
  K = 4
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3, 5, 10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))

  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = bsdPar(n , mti = p  , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        Teilsequenz2 = gensequence$M[1:((n/K)*2)]
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = 2*n/K - no_T_2

        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

        ######################### Erwartungswert für Teststatistik in der 3. Stufe ##########################
        Teilsequenz3 = gensequence$M[1:((n/K)*3)]
        no_T_3 = sum(Teilsequenz3)
        no_C_3 = 3*n/K - no_T_3

        nu.stage3 = calc.nu(Teilsequenz3)
        vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
        factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
        mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





        ############################ Erwartungswert für Teststatistik in der 4. Stufe  ##############################
        no_T_4 = sum(gensequence@M)
        no_C_4 = n - no_T_4

        nu.stage4 = calc.nu(gensequence@M)
        vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
        factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
        mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta


        ########################### Kovarianzen von Z1 und Z2 ###########################################
        cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z1 und Z3 ###########################################
        cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z1 und Z4 ###########################################
        cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



        ########################### Kovarianzen von Z2 und Z3 ###########################################
        cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
        ########################### Kovarianzen von Z2 und Z4 ###########################################
        cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z3 und Z4 ###########################################
        cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






        #### Kovarianzmatrix #############
        Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
        Sigma_inv = solve(Sigma)
        det_Sigma = det(Sigma)
        Sigma_3 = Sigma[1:3, 1:3]
        Sigma_inv_3 = solve(Sigma_3)
        det_Sigma_3 = det(Sigma_3)




        gpoc = 2.361

        #### Fehlerwahrscheinlichkeiten #####

        ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -gpoc, gpoc)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1


        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
        wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2

        ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
        wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
        wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
        alpha_3 = alpha_2 + wert_4 + wert_5

        ergebnis.vector.alpha.3[i] = alpha_3

        ######################### Fehlerwahrscheinlichkeit der 4. Stufe ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
        wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
        wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
        alpha_4 = alpha_3 + wert_6 + wert_7

        ergebnis.vector.alpha.4[i] = alpha_4




      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      }
    }
  }



  cat("Pocock 4 Stages using BSD(b) b=3/5/10 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3),
                   round(Ergebnismatrix_3, digits = 3))
  Ergebnis
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("b=3 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05",
                               "b=5 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05",
                               "b=10 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}

POC_2_BSD_A2 <- function(n,norep) {
  #############         Pocock, 2 Stages  Big Stick Design, applying the randomization procedure to each stage






  #parameters
  K = 2
  alpha = 0.05
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3,5,10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))


  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      grpoc = 2.178
      for (i in 1:norep)
      {
        randseq = bsdPar(n/K , mti = p, groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



        ################################### Fehlerwahrscheinlichkeit in der 1.Stufe ##################
        #
        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.178, 2.178)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1



        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        randseq_2 = bsdPar(n/K , mti = p, groups = c("0","1"))
        gensequence_2 = genSeq(randseq_2,1)
        Teilsequenz2 = c(gensequence$M, gensequence_2$M)
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = n - no_T_2


        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

        ######################## Kovarianz #######################
        ### Kovarianz von Z1 und Z2

        cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
        wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2
        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
        wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2
      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }

    }

  }


  cat("Pocock 2 Stages using BSD(b) b=3/5/10 - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2", "p<=0.05", "p<=0.05", "p=5 1", "2", "p<=0.05", "p<=0.05", "p=10 1", "2", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}

POC_3_BSD_A2 <- function(n,norep) {
  #############         Pocock, 3 Stages  Big Stick Design, applying the randomization procedure to each stage






  #parameters
  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3,5,10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))


  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = bsdPar(n/K , mti = p , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        randseq_2 = bsdPar(n/K , mti = p, groups = c("0","1"))
        gensequence_2 = genSeq(randseq_2,1)
        Teilsequenz2 = c(gensequence$M, gensequence_2$M)
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = 2*n/K - no_T_2

        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





        ############################ Erwartungswert für Teststatistik in der 3. Stufe  ##############################
        randseq_3 = bsdPar(n/K , mti = p, groups = c("0","1"))
        gensequence_3 = genSeq(randseq_3,1)
        Teilsequenz3 = c(Teilsequenz2, gensequence_3$M)
        no_T_3 = sum(Teilsequenz3)
        no_C_3 = n - no_T_3

        nu.stage3 = calc.nu(Teilsequenz3)
        vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
        factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
        mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


        ########################### Kovarianzen von Z1 und Z2 ###########################################
        cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z1 und Z3 ###########################################
        cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z2 und Z3 ###########################################
        cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

        #### Kovarianzmatrix #############
        Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
        Sigma_inv = solve(Sigma)
        det_Sigma = det(Sigma)






        #### Fehlerwahrscheinlichkeiten #####

        ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.289, 2.289)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1


        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
        wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2

        ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
        wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
        wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
        alpha_3 = alpha_2 + wert_4 + wert_5

        ergebnis.vector.alpha.3[i] = alpha_3




      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }
    }
  }

  cat("Pocock 3 Stages using BSD(b) b=3/5/10 - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05",
                               "p=5 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05", "p=10 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}

POC_4_BSD_A2 <- function(n,norep) {
  #############         Pocock, 4 Stages  Big Stick Design, applying the randomization procedure to each stage






  #parameters
  K = 4
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3, 5, 10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))

  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = bsdPar(n/K , mti = p  , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        randseq_2 = bsdPar(n/K , mti = p, groups = c("0","1"))
        gensequence_2 = genSeq(randseq_2,1)
        Teilsequenz2 = c(gensequence$M, gensequence_2$M)
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = 2*n/K - no_T_2


        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

        ######################### Erwartungswert für Teststatistik in der 3. Stufe ##########################
        randseq_3 = bsdPar(n/K , mti = p, groups = c("0","1"))
        gensequence_3 = genSeq(randseq_2,1)
        Teilsequenz3 = c(Teilsequenz2, gensequence_3$M)
        no_T_3 = sum(Teilsequenz3)
        no_C_3 = 3*n/K - no_T_3

        nu.stage3 = calc.nu(Teilsequenz3)
        vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
        factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
        mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





        ############################ Erwartungswert für Teststatistik in der 4. Stufe  ##############################
        randseq_4 = bsdPar(n/K , mti = p, groups = c("0","1"))
        gensequence_4 = genSeq(randseq_4,1)
        Teilsequenz4 = c(Teilsequenz3, gensequence_4$M)
        no_T_4 = sum(Teilsequenz4)
        no_C_4 = n - no_T_4

        nu.stage4 = calc.nu(Teilsequenz4)
        vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
        factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
        mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta


        ########################### Kovarianzen von Z1 und Z2 ###########################################
        cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z1 und Z3 ###########################################
        cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z1 und Z4 ###########################################
        cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



        ########################### Kovarianzen von Z2 und Z3 ###########################################
        cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
        ########################### Kovarianzen von Z2 und Z4 ###########################################
        cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z3 und Z4 ###########################################
        cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






        #### Kovarianzmatrix #############
        Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
        Sigma_inv = solve(Sigma)
        det_Sigma = det(Sigma)
        Sigma_3 = Sigma[1:3, 1:3]
        Sigma_inv_3 = solve(Sigma_3)
        det_Sigma_3 = det(Sigma_3)




        gpoc = 2.361

        #### Fehlerwahrscheinlichkeiten #####

        ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -gpoc, gpoc)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1


        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
        wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2

        ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
        wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
        wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
        alpha_3 = alpha_2 + wert_4 + wert_5

        ergebnis.vector.alpha.3[i] = alpha_3

        ######################### Fehlerwahrscheinlichkeit der 4. Stufe ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
        wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
        wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
        alpha_4 = alpha_3 + wert_6 + wert_7

        ergebnis.vector.alpha.4[i] = alpha_4




      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      }
    }
  }



  cat("Pocock 4 Stages using BSD(b) b=3/5/10 - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3),
                   round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("b=3 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05",
                               "b=5 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05",
                               "b=10 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}


###### Group Sequential Design Pocock using EBC(p) ###########

POC_2_EBC <- function(n,norep) {
  ################### Pocock, 2 Stages, Efron's Biased Coin, applying randomization procedure to all patients










  #Parameters
  K = 2
  alpha = 0.5
  sigma_sqrt = 1
  effect = 0
  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)
  p_vector = c(0.5, 2/3, 0.8)



  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+2)), nrow = length(eta_vector))


  #Calculating Type 1 Errors:

  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      grpoc = 2.178
      for (i in 1:norep)
      {
        randseq = ebcPar(n , p , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1


        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {

            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Expected Value of Z_1 ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



        ################################### Type 1 Error 1st Stage ##################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.178, 2.178)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1



        #########################    Expected Value of Z_2  ##########################
        no_T_2 = sum(gensequence@M)
        no_C_2 = n - no_T_2

        nu.stage2 = calc.nu(gensequence@M)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

        ########################     Covariance    #######################
        ###Covariance of Z1 and  Z2

        cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



        ######################### Type 1 Error of 2nd Stage  ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
        wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2
      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)

      }else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)

      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)

      }
    }
  }



  cat("POC 2 Stages using EBC(p) for p=0.5/(2/3)/0.8 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=0.5 1", "2", "p<=0.05", "p<=0.05", "p=(2/3) 1", "2", "p<=0.05", "p<=0.05", "p=0.8 1", "2", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}


POC_3_EBC <- function(n, norep) {
  ################### Pocock, 3 Stages, Efron's Biased Coin, applying randomization procedure to all patients









  #parameters
  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(0.5, 2/3, 0.8)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+3)), nrow = length(eta_vector))


  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = ebcPar(n , p , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {

            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Expected Value of Z_1   ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




        ######################### Expected Value of Z_2 ##########################
        Teilsequenz2 = gensequence$M[1:((n/K)*2)]
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = 2*n/K - no_T_2

        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





        ############################ Expected Value of Z_3  ##############################
        no_T_3 = sum(gensequence@M)
        no_C_3 = n - no_T_3

        nu.stage3 = calc.nu(gensequence@M)
        vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
        factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
        mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


        ########################### Covariance of Z_1 and Z_2 ###########################################
        cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

        ###########################  Covariance of Z_1 and Z_3 ###########################################
        cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

        ###########################  Covariance of Z_2 and Z_3 ###########################################
        cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

        ###Covariance Matrix  #############
        Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
        Sigma_inv = solve(Sigma)
        det_Sigma = det(Sigma)






        #### Type 1 Errors #####

        ######################### Type 1 Error 1st Stage##################################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.289, 2.289)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1


        ######################### Type 1 Error 2nd Stage ##################################
        #
        f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
        wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2

        ######################### Type 1 Error 3rd Stage ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
        wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
        wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
        alpha_3 = alpha_2 + wert_4 + wert_5

        ergebnis.vector.alpha.3[i] = alpha_3




      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }
    }
  }


  cat("POC 3 Stages using EBC(p) for p=0.5/(2/3)/0.8 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=0.5 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05",
                               "p=(2/3) 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05", "p=0.8 1", "2", "3", "p<=0.05", "p<=0.05", "p<=0.05"))

  print(Ergebnis)
}


POC_4_EBC <- function(n, norep) {

  ################### Pocock, 3 Stages, Efron's Biased Coin, applying randomization procedure to all patients








  #parameters
  K = 4
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3, 5, 10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K+4)), nrow = length(eta_vector))

  for (r in 1:length(p_vector))
  {
    p = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = bsdPar(n , mti = p  , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {

            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Expected Value of Z_1  ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




        #########################Expected Value of Z_2 ##########################
        Teilsequenz2 = gensequence$M[1:((n/K)*2)]
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = 2*n/K - no_T_2

        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

        ######################### Expected Value of Z_3 ##########################
        Teilsequenz3 = gensequence$M[1:((n/K)*3)]
        no_T_3 = sum(Teilsequenz3)
        no_C_3 = 3*n/K - no_T_3

        nu.stage3 = calc.nu(Teilsequenz3)
        vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
        factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
        mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





        ###########################Expected Value of Z_4  ##############################
        no_T_4 = sum(gensequence@M)
        no_C_4 = n - no_T_4

        nu.stage4 = calc.nu(gensequence@M)
        vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
        factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
        mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta


        ########################### Covariance of Z_1 and Z_2 ###########################################
        cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Covariance of Z_1 and Z_3 ###########################################
        cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Covariance of Z_1 and Z_4 ###########################################
        cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



        ########################### Covariance of Z_2 and Z_3 ###########################################
        cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
        ########################### Covariance of Z_2 and Z_4 ###########################################
        cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Covariance of Z_3 and Z_4 ###########################################
        cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
        korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






        #### Covariance Matrix #############
        Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
        Sigma_inv = solve(Sigma)
        det_Sigma = det(Sigma)
        Sigma_3 = Sigma[1:3, 1:3]
        Sigma_inv_3 = solve(Sigma_3)
        det_Sigma_3 = det(Sigma_3)




        gpoc = 2.361

        #### Type 1 Errors  #####

        #########################  Type 1 Error 1st Stage ##################################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -gpoc, gpoc)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1


        ######################### Type 1 Error 2nd Stage ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
        wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2

        #########################  Type 1 Error 3rd Stage ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
        wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
        wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
        alpha_3 = alpha_2 + wert_4 + wert_5

        ergebnis.vector.alpha.3[i] = alpha_3

        #########################  Type 1 Error 4th Stage ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
        wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
        wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
        alpha_4 = alpha_3 + wert_6 + wert_7

        ergebnis.vector.alpha.4[i] = alpha_4




      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = mean(ergebnis.vector.alpha.4)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
      }
    }
  }



  cat("POC 4 Stages using EBC(p) for p=0.5/(2/3)/0.8 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3),
                   round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=0.5 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05",
                               "p=(2/3) 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05",
                               "p=0.8 1", "2", "3", "4", "p<=0.05", "p<=0.05", "p<=0.05", "p<=0.05"))

  print(Ergebnis)

}


################ Group Sequential Design Pocock using BCDWIT(2/3,b) #################


POC_2_CHEN_A1 <- function(n,norep) {
  #############            Pocock, 2 Stages, Chen's Design, applying the randomization procedure to all patients






  #parameters
  K = 2
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3, 5, 10)

  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  for (r in 1:length(p_vector))
  {
    b = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      grpoc = 2.178
      for (i in 1:norep)
      {
        randseq = chenPar(n , mti = b , p = 2/3,  groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



        ################################### Fehlerwahrscheinlichkeit in der 1.Stufe ##################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.178, 2.178)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1



        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        no_T_2 = sum(gensequence@M)
        no_C_2 = n - no_T_2

        nu.stage2 = calc.nu(gensequence@M)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

        ######################## Kovarianz #######################


        cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
        wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2
      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }
    }
  }






  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2", "p<=0.05", "p<=0.05", "p=5 1", "2", "p<=0.05", "p<=0.05", "p=10 1", "2", "p<=0.05", "p<=0.05"))

  cat("POC 2 Stages using BCDWIT(2/3,b) for b=3/5/10 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  print(Ergebnis)
}

POC_3_CHEN_A1 <- function(n,norep) {
  #############            Pocock, 3 Stages, Chen's Design, applying the randomization procedure to all patients






  #parameters
  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3,5,10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))


  for (r in 1:length(p_vector))
  {
    b = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = chenPar(n , mti = b, p = 2/3 , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            #Case 1: n_T = n_C und seq = 0 (=Control)
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
        Teilsequenz2 = gensequence$M[1:((n/K)*2)]
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = 2*n/K - no_T_2

        nu.stage2 = calc.nu(Teilsequenz2)
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





        ############################ Erwartungswert für Teststatistik in der 3. Stufe  ##############################
        no_T_3 = sum(gensequence@M)
        no_C_3 = n - no_T_3

        nu.stage3 = calc.nu(gensequence@M)
        vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
        factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
        mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


        ########################### Kovarianzen von Z1 und Z2 ###########################################
        cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z1 und Z3 ###########################################
        cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

        ########################### Kovarianzen von Z2 und Z3 ###########################################
        cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
        korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

        #### Kovarianzmatrix #############
        Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
        Sigma_inv = solve(Sigma)
        det_Sigma = det(Sigma)






        #### Fehlerwahrscheinlichkeiten #####

        ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.289, 2.289)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1


        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
        wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2

        ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

        mu = c(mu_1, mu_2, mu_3)
        f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
        wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
        wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
        alpha_3 = alpha_2 + wert_4 + wert_5

        ergebnis.vector.alpha.3[i] = alpha_3




      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }
    }
  }



  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2","3", "p<=0.05", "p<=0.05","p<=0.05", "p=5 1", "2","3", "p<=0.05","p<=0.05", "p<=0.05", "p=10 1", "2","3", "p<=0.05","p<=0.05", "p<=0.05"))
  cat("POC 3 Stages using BCDWIT(2/3,b) for b=3/5/10 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  print(Ergebnis)
}

POC_4_CHEN_A1 <- function(n,norep) {
  #############            Pocock, 4 Stages, Chen's Design, applying the randomization procedure to all patients






  #parameters
  no_failure_1 = 0
  K = 4
  sigma_sqrt = 1
  effect = 0

  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3, 5, 10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))

  for (r in 1:length(p_vector))
  {
    b = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = chenPar(n , mti = b, p = 2/3  , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M[1:(n/K)]

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }

        if(no_C_1 == 0 || no_T_1 == 0)
        {
          no_failure_1 = no_failure_1 + 1
        }else{
          nu.stage1 = calc.nu(subsequence_1)
          nu.stage1


          ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
          vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
          factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
          mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




          ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
          Teilsequenz2 = gensequence$M[1:((n/K)*2)]
          no_T_2 = sum(Teilsequenz2)
          no_C_2 = 2*n/K - no_T_2

          nu.stage2 = calc.nu(Teilsequenz2)
          vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
          factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
          mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta

          ######################### Erwartungswert für Teststatistik in der 3. Stufe ##########################
          Teilsequenz3 = gensequence$M[1:((n/K)*3)]
          no_T_3 = sum(Teilsequenz3)
          no_C_3 = 3*n/K - no_T_3

          nu.stage3 = calc.nu(Teilsequenz3)
          vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
          factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
          mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





          ############################ Erwartungswert für Teststatistik in der 4. Stufe  ##############################
          no_T_4 = sum(gensequence@M)
          no_C_4 = n - no_T_4

          nu.stage4 = calc.nu(gensequence@M)
          vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
          factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
          mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta


          ########################### Kovarianzen von Z1 und Z2 ###########################################
          cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
          korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z1 und Z3 ###########################################
          cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
          korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z1 und Z4 ###########################################
          cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
          korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



          ########################### Kovarianzen von Z2 und Z3 ###########################################
          cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
          korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
          ########################### Kovarianzen von Z2 und Z4 ###########################################
          cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
          korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z3 und Z4 ###########################################
          cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
          korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






          #### Kovarianzmatrix #############
          Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
          Sigma_inv = solve(Sigma)
          det_Sigma = det(Sigma)
          Sigma_3 = Sigma[1:3, 1:3]
          Sigma_inv_3 = solve(Sigma_3)
          det_Sigma_3 = det(Sigma_3)




          gpoc = 2.361

          #### Fehlerwahrscheinlichkeiten #####

          ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################

          ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
          wert_1 = integrate(ftest_1, -gpoc, gpoc)
          alpha_1=  1-wert_1$value
          ergebnis.vector.alpha.1[i] = alpha_1


          ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

          f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
          wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
          wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
          alpha_2 = alpha_1 + wert_2 + wert_3

          ergebnis.vector.alpha.2[i] = alpha_2

          ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

          mu = c(mu_1, mu_2, mu_3)
          f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
          wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
          wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
          alpha_3 = alpha_2 + wert_4 + wert_5

          ergebnis.vector.alpha.3[i] = alpha_3

          ######################### Fehlerwahrscheinlichkeit der 4. Stufe ##################################

          mu = c(mu_1, mu_2, mu_3)
          f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
          wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
          wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
          alpha_4 = alpha_3 + wert_6 + wert_7

          ergebnis.vector.alpha.4[i] = alpha_4




        }
        if (r == 1)
        {
          Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
          Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
          Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
          Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
          Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
          Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
          Ergebnismatrix_1[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
          Ergebnismatrix_1[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
        } else if (r == 2)
        {
          Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
          Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
          Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
          Ergebnismatrix_2[m,4] = mean(ergebnis.vector.alpha.4)
          Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
          Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
          Ergebnismatrix_2[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
          Ergebnismatrix_2[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
        }else
        {
          Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
          Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
          Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
          Ergebnismatrix_3[m,4] = mean(ergebnis.vector.alpha.4)
          Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
          Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
          Ergebnismatrix_3[m,7] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
          Ergebnismatrix_3[m,8] = sum(ergebnis.vector.alpha.4<=0.05)/length(ergebnis.vector.alpha.4)
        }
      }
    }
  }



  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3),
                   round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2","3","4", "p<=0.05","p<=0.05", "p<=0.05","p<=0.05", "p=5 1", "2","3","4", "p<=0.05","p<=0.05","p<=0.05", "p<=0.05", "p=10 1", "2","3","4", "p<=0.05","p<=0.05","p<=0.05", "p<=0.05"))
  cat("POC 4 Stages using BCDWIT(2/3,b) for b=3/5/10 - applying the randomization procedure to all patients \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  print(Ergebnis)
}

POC_2_CHEN_A2 <- function(n,norep) {
  #############            Pocock, 2 Stages, Chen's Design, applying the randomization procedure to each stage






  #parameters
  K = 2
  alpha = 0.05
  sigma_sqrt = 1
  effect = 0

  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  alpha_1 = vector()
  alpha_2 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3,5,10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))


  for (r in 1:length(p_vector))
  {
    b = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      grpoc = 2.178
      for (i in 1:norep)
      {
        randseq = chenPar(n/K , mti = b, p = 2/3, groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1




        calc.nu = function(seq)
        {
          nu = rep(0,6)
          noT = 0
          noC = 0
          for (i in 1: length(seq))
          {
            if ( (noT == noC) & (seq[i] == 0))
            {noC = noC +1
            nu[2] = nu[2]+1
            } else if ( (noT > noC) & (seq[i] == 0) )
            {noC = noC +1
            nu[1] = nu[1]+1
            } else if ( (noT < noC) & (seq[i] == 0))
            {noC = noC +1
            nu[3] = nu[3]+1
            }else if ( (noT == noC) & (seq[i] == 1))
            {noT = noT +1
            nu[5] = nu[5]+1
            }else if ( (noT > noC) & (seq[i] == 1))
            {noT = noT +1
            nu[4] = nu[4]+1
            }else if ( (noT < noC) & (seq[i] == 1))
            {noT = noT +1
            nu[6] = nu[6]+1
            }


          }
          return(nu)

        }


        nu.stage1 = calc.nu(subsequence_1)
        nu.stage1


        ############################ Erwartungswert für Teststatistik ##############################
        vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
        factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
        mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta



        ################################### Fehlerwahrscheinlichkeit in der 1.Stufe ##################

        ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
        wert_1 = integrate(ftest_1, -2.178, 2.178)
        alpha_1=  1-wert_1$value
        ergebnis.vector.alpha.1[i] = alpha_1



        ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################

        randseq_2 = chenPar(n/K , mti = b, p=2/3, groups = c("0","1"))
        gensequence_2 = genSeq(randseq_2,1)

        nu_stage_temp = calc.nu(gensequence_2$M)
        Teilsequenz2 = c(gensequence$M, gensequence_2$M)
        no_T_2 = sum(Teilsequenz2)
        no_C_2 = n - no_T_2


        nu.stage2 = nu.stage1 + nu_stage_temp
        vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
        factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
        mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta


        ### Kovarianz von Z1 und Z2

        cov = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
        korr = cov/sqrt((sigma_sqrt * sigma_sqrt))



        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
        wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2
        ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

        f <- function(x) (1/(2*pi* sqrt(1-korr^2) )) * exp( -1/(2*(1-korr^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr * (x[1]-mu_1)*(x[2] - mu_2)) )
        wert_2 = cubature::pcubature(f, c(-2.178, -Inf), c(2.178, -2.178))$integral
        wert_3 = cubature::pcubature(f, c(-2.178, 2.178), c(2.178, Inf))$integral
        alpha_2 = alpha_1 + wert_2 + wert_3

        ergebnis.vector.alpha.2[i] = alpha_2
      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
      }

    }

  }


  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3),
                   round(Ergebnismatrix_3, digits = 3))
  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2","p<=0.05","p<=0.05", "p=5 1", "2","p<=0.05", "p<=0.05", "p=10 1", "2","p<=0.05", "p<=0.05"))
  cat("POC 2 Stages using BCDWIT(2/3,b) for b=3/5/10 - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  print(Ergebnis)
}

POC_3_CHEN_A2 <- function(n,norep) {
  #############            Pocock, 3 Stages, Chen's Design, applying the randomization procedure to each stage






  #parameters
  no_failure_1 = 0
  K = 3
  sigma_sqrt = 1
  effect = 0
  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()


  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3,5,10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))


  for (r in 1:length(p_vector))
  {
    b = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = chenPar(n/K , mti = b, p = 2/3 , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1

        if (no_T_1 == 0 || no_C_1 == 0)
        {
          no_failure_1 = no_failure_1 + 1

        }else
        {




          calc.nu = function(seq)
          {
            nu = rep(0,6)
            noT = 0
            noC = 0
            for (i in 1: length(seq))
            {
              if ( (noT == noC) & (seq[i] == 0))
              {noC = noC +1
              nu[2] = nu[2]+1
              } else if ( (noT > noC) & (seq[i] == 0) )
              {noC = noC +1
              nu[1] = nu[1]+1
              } else if ( (noT < noC) & (seq[i] == 0))
              {noC = noC +1
              nu[3] = nu[3]+1
              }else if ( (noT == noC) & (seq[i] == 1))
              {noT = noT +1
              nu[5] = nu[5]+1
              }else if ( (noT > noC) & (seq[i] == 1))
              {noT = noT +1
              nu[4] = nu[4]+1
              }else if ( (noT < noC) & (seq[i] == 1))
              {noT = noT +1
              nu[6] = nu[6]+1
              }


            }
            return(nu)

          }


          nu.stage1 = calc.nu(subsequence_1)
          nu.stage1


          ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
          vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
          factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
          mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta




          ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
          randseq_2 = chenPar(n/K , mti = b, p=2/3, groups = c("0","1"))
          gensequence_2 = genSeq(randseq_2,1)

          nu_stage_temp = calc.nu(gensequence_2$M)
          Teilsequenz2 = c(gensequence$M, gensequence_2$M)
          no_T_2 = sum(Teilsequenz2)
          no_C_2 = 2*(n/K) - no_T_2


          nu.stage2 = nu.stage1 + nu_stage_temp
          vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
          factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
          mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





          ############################ Erwartungswert für Teststatistik in der 3. Stufe  ##############################
          randseq_3 = chenPar(n/K , mti = b, p = 2/3 , groups = c("0","1"))
          gensequence_3 = genSeq(randseq_3,1)
          nu_stage_temp_2 = calc.nu(gensequence_3$M)
          Teilsequenz3 = c(Teilsequenz2, gensequence_3$M)
          no_T_3 = sum(Teilsequenz3)
          no_C_3 = n - no_T_3

          nu.stage3 = nu.stage2 + nu_stage_temp_2
          vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
          factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
          mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta


          ########################### Kovarianzen von Z1 und Z2 ###########################################
          cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
          korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z1 und Z3 ###########################################
          cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
          korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z2 und Z3 ###########################################
          cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
          korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))

          #### Kovarianzmatrix #############
          Sigma = matrix( c(1, cov_12, cov_13, cov_12, 1, cov_23, cov_13, cov_23, 1), nrow = 3)
          Sigma_inv = solve(Sigma)
          det_Sigma = det(Sigma)






          #### Fehlerwahrscheinlichkeiten #####

          ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################

          ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
          wert_1 = integrate(ftest_1, -2.289, 2.289)
          alpha_1=  1-wert_1$value
          ergebnis.vector.alpha.1[i] = alpha_1


          ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################

          f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
          wert_2 = cubature::pcubature(f, c(-2.289, -Inf), c(2.289, -2.289))$integral
          wert_3 = cubature::pcubature(f, c(-2.289, 2.289), c(2.289, Inf))$integral
          alpha_2 = alpha_1 + wert_2 + wert_3

          ergebnis.vector.alpha.2[i] = alpha_2

          ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

          mu = c(mu_1, mu_2, mu_3)
          f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
          wert_4 = cubature::pcubature(f, c(-2.289, -2.289, -Inf), c(2.289, 2.289, -2.289))$integral
          wert_5 = cubature::pcubature(f, c(-2.289, -2.289, 2.289), c(2.289, 2.289, Inf))$integral
          alpha_3 = alpha_2 + wert_4 + wert_5

          ergebnis.vector.alpha.3[i] = alpha_3

        }


      }
      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = sum(ergebnis.vector.alpha.1<=0.05)/length(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,5] = sum(ergebnis.vector.alpha.2<=0.05)/length(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,6] = sum(ergebnis.vector.alpha.3<=0.05)/length(ergebnis.vector.alpha.3)
      }
    }
  }


  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))

  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2","3","p<=0.05","p<=0.05","p<=0.05", "p=5 1", "2","3","p<=0.05","p<=0.05", "p<=0.05", "p=10 1", "2","3","p<=0.05","p<=0.05", "p<=0.05"))
  cat("POC 3 Stages using BCDWIT(2/3,b) for b=3/5/10 - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  print(Ergebnis)
}

POC_4_CHEN_A2 <- function(n,norep) {
  #############            Pocock, 4 Stages, Chen's Design, applying the randomization procedure to each stage






  #parameters
  no_failure = 0
  K = 4
  sigma_sqrt = 1
  effect = 0

  ergebnis.vector.alpha.2 = vector()
  ergebnis.vector.alpha.1 = vector()
  ergebnis.vector.alpha.3 = vector()
  ergebnis.vector.alpha.4 = vector()
  alpha_1 = vector()
  alpha_2 = vector()
  alpha_3 = vector()
  alpha_4 = vector()

  eta_vector = c(0, 0.01,0.1,0.2,0.3,0.4,0.5)

  p_vector = c(3, 5, 10)


  Ergebnismatrix_1= matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_2 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))
  Ergebnismatrix_3 = matrix(rep(0,length(eta_vector)*(K*2)), nrow = length(eta_vector))

  for (r in 1:length(p_vector))
  {
    b = p_vector[r]

    for (m in 1:length(eta_vector))
    {


      eta = eta_vector[m]
      set.seed(NULL)

      for (i in 1:norep)
      {
        randseq = chenPar(n/K , mti = b, p = 2/3  , groups = c("0","1"))
        gensequence = genSeq(randseq,1)
        subsequence_1 = gensequence$M

        no_T_1 = sum(subsequence_1)
        no_C_1 = n/K-no_T_1


        if (no_T_1 == 0 || no_C_1 == 0)
        {
          no_failure = no_failure + 1

        }else
        {

          calc.nu = function(seq)
          {
            nu = rep(0,6)
            noT = 0
            noC = 0
            for (i in 1: length(seq))
            {
              if ( (noT == noC) & (seq[i] == 0))
              {noC = noC +1
              nu[2] = nu[2]+1
              } else if ( (noT > noC) & (seq[i] == 0) )
              {noC = noC +1
              nu[1] = nu[1]+1
              } else if ( (noT < noC) & (seq[i] == 0))
              {noC = noC +1
              nu[3] = nu[3]+1
              }else if ( (noT == noC) & (seq[i] == 1))
              {noT = noT +1
              nu[5] = nu[5]+1
              }else if ( (noT > noC) & (seq[i] == 1))
              {noT = noT +1
              nu[4] = nu[4]+1
              }else if ( (noT < noC) & (seq[i] == 1))
              {noT = noT +1
              nu[6] = nu[6]+1
              }


            }
            return(nu)

          }


          nu.stage1 = calc.nu(subsequence_1)
          nu.stage1


          ############################ Erwartungswert für Teststatistik in der 1. Stufe  ##############################
          vector.exp = c(1/no_C_1, 0, -1/no_C_1, -1/no_T_1, 0 , 1/no_T_1)
          factor = as.numeric((vector.exp %*% nu.stage1)[1,1])
          mu_1 = sqrt( 1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1)) ) * factor * eta





          ######################### Erwartungswert für Teststatistik in der 2. Stufe ##########################
          randseq_2 = chenPar(n/K , mti = b, p=2/3, groups = c("0","1"))
          gensequence_2 = genSeq(randseq_2,1)

          nu_stage_temp = calc.nu(gensequence_2$M)
          Teilsequenz2 = c(gensequence$M, gensequence_2$M)
          no_T_2 = sum(Teilsequenz2)
          no_C_2 = 2*n/K - no_T_2


          nu.stage2 = nu.stage1 + nu_stage_temp
          vector.exp.secstage = c(1/no_C_2, 0, -1/no_C_2, -1/no_T_2, 0 , 1/no_T_2)
          factor2 = as.numeric((vector.exp.secstage %*% nu.stage2)[1,1])
          mu_2 = sqrt( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2)) ) * factor2 * eta





          ############################ Erwartungswert für Teststatistik in der 3. Stufe  ##############################
          randseq_3 = chenPar(n/K , mti = b, p = 2/3 , groups = c("0","1"))
          gensequence_3 = genSeq(randseq_3,1)
          nu_stage_temp_2 = calc.nu(gensequence_3$M)
          Teilsequenz3 = c(Teilsequenz2, gensequence_3$M)
          no_T_3 = sum(Teilsequenz3)
          no_C_3 = 3*n/K - no_T_3

          nu.stage3 = nu.stage2 + nu_stage_temp_2
          vector.exp.thirdstage = c(1/no_C_3, 0, -1/no_C_3, -1/no_T_3, 0 , 1/no_T_3)
          factor3 = as.numeric((vector.exp.thirdstage %*% nu.stage3)[1,1])
          mu_3 = sqrt( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3)) ) * factor3 * eta





          ############################ Erwartungswert für Teststatistik in der 4. Stufe  ##############################
          randseq_4 = chenPar(n/K , mti = b, p = 2/3 , groups = c("0","1"))
          gensequence_4 = genSeq(randseq_4,1)
          nu_stage_temp_3 = calc.nu(gensequence_4$M)
          Teilsequenz4 = c(Teilsequenz3, gensequence_4$M)
          no_T_4 = sum(Teilsequenz4)
          no_C_4 = n - no_T_4

          nu.stage4 = nu.stage3 + nu_stage_temp_3
          vector.exp.fourthstage = c(1/no_C_4, 0, -1/no_C_4, -1/no_T_4, 0 , 1/no_T_4)
          factor4 = as.numeric((vector.exp.fourthstage %*% nu.stage4)[1,1])
          mu_4 = sqrt( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4)) ) * factor4 * eta



          ########################### Kovarianzen von Z1 und Z2 ###########################################
          cov_12 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) ) * ( sigma_sqrt/no_C_2 + sigma_sqrt/no_T_2)
          korr_12 = cov_12/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z1 und Z3 ###########################################
          cov_13 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
          korr_13 = cov_13/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z1 und Z4 ###########################################
          cov_14 = sqrt(  (1/ (sigma_sqrt* (1/no_C_1 + 1/no_T_1))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
          korr_14 = cov_14/sqrt((sigma_sqrt * sigma_sqrt))



          ########################### Kovarianzen von Z2 und Z3 ###########################################
          cov_23 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) ) * ( sigma_sqrt/no_C_3 + sigma_sqrt/no_T_3)
          korr_23 = cov_23/sqrt((sigma_sqrt * sigma_sqrt))
          ########################### Kovarianzen von Z2 und Z4 ###########################################
          cov_24 = sqrt(  (1/ (sigma_sqrt* (1/no_C_2 + 1/no_T_2))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
          korr_24 = cov_24/sqrt((sigma_sqrt * sigma_sqrt))

          ########################### Kovarianzen von Z3 und Z4 ###########################################
          cov_34 = sqrt(  (1/ (sigma_sqrt* (1/no_C_3 + 1/no_T_3))) * ( 1/ (sigma_sqrt* (1/no_C_4 + 1/no_T_4))) ) * ( sigma_sqrt/no_C_4 + sigma_sqrt/no_T_4)
          korr_34 = cov_34/sqrt((sigma_sqrt * sigma_sqrt))






          #### Kovarianzmatrix #############
          Sigma = matrix( c(1, cov_12, cov_13, cov_14, cov_12, 1, cov_23, cov_24, cov_13, cov_23, 1, cov_34, cov_14, cov_24, cov_34, 1), nrow = 4)
          Sigma_inv = solve(Sigma)
          det_Sigma = det(Sigma)
          Sigma_3 = Sigma[1:3, 1:3]
          Sigma_inv_3 = solve(Sigma_3)
          det_Sigma_3 = det(Sigma_3)




          gpoc = 2.361

          #### Fehlerwahrscheinlichkeiten #####

          ######################### Fehlerwahrscheinlichkeit der 1. Stufe ##################################
          #
          ftest_1 = function(x) 1/(sqrt(2*pi)) * exp( - (( x-mu_1)^2) / 2 )
          wert_1 = integrate(ftest_1, -gpoc, gpoc)
          alpha_1=  1-wert_1$value
          ergebnis.vector.alpha.1[i] = alpha_1


          ######################### Fehlerwahrscheinlichkeit der 2. Stufe ##################################
          #
          f <- function(x) (1/(2*pi* sqrt(1-korr_12^2) )) * exp( -1/(2*(1-korr_12^2))  * (  (x[1]-(mu_1))^2 +  (x[2]-(mu_2))^2 - 2* korr_12 * (x[1]-mu_1)*(x[2] - mu_2)) )
          wert_2 = cubature::pcubature(f, c(-gpoc, -Inf), c(gpoc, -gpoc))$integral
          wert_3 = cubature::pcubature(f, c(-gpoc, gpoc), c(gpoc, Inf))$integral
          alpha_2 = alpha_1 + wert_2 + wert_3

          ergebnis.vector.alpha.2[i] = alpha_2

          ######################### Fehlerwahrscheinlichkeit der 3. Stufe ##################################

          mu = c(mu_1, mu_2, mu_3)
          f <- function(x) (1/ sqrt((2*pi)^3 * det_Sigma_3)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3)) %*% Sigma_inv_3 %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3))
          wert_4 = hcubature(f, c(-gpoc, -gpoc, -Inf), c(gpoc, gpoc, -gpoc))$integral
          wert_5 = hcubature(f, c(-gpoc, -gpoc, gpoc), c(gpoc, gpoc, Inf))$integral
          alpha_3 = alpha_2 + wert_4 + wert_5

          ergebnis.vector.alpha.3[i] = alpha_3

          ######################### Fehlerwahrscheinlichkeit der 4. Stufe ##################################

          mu = c(mu_1, mu_2, mu_3)
          f <- function(x) (1/ sqrt((2*pi)^4 * det_Sigma)) * exp( -0.5 * t(c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4)) %*% Sigma_inv %*% c(x[1]-mu_1, x[2]-mu_2, x[3]-mu_3, x[4]-mu_4))
          wert_6 = hcubature(f, c(-gpoc, -gpoc, -gpoc, -Inf), c(gpoc, gpoc, gpoc, -gpoc))$integral
          wert_7 = hcubature(f, c(-gpoc, -gpoc, -gpoc, gpoc), c(gpoc, gpoc, gpoc, Inf))$integral
          alpha_4 = alpha_3 + wert_6 + wert_7

          ergebnis.vector.alpha.4[i] = alpha_4


        }
      }


      if (r == 1)
      {
        Ergebnismatrix_1[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_1[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_1[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_1[m,4] = mean(ergebnis.vector.alpha.4)
      } else if (r == 2)
      {
        Ergebnismatrix_2[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_2[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_2[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_2[m,4] = mean(ergebnis.vector.alpha.4)
      }else
      {
        Ergebnismatrix_3[m,1] =  mean(ergebnis.vector.alpha.1)
        Ergebnismatrix_3[m,2] = mean(ergebnis.vector.alpha.2)
        Ergebnismatrix_3[m,3] = mean(ergebnis.vector.alpha.3)
        Ergebnismatrix_3[m,4] = mean(ergebnis.vector.alpha.4)
      }


    }
  }



  Ergebnis = cbind(round(Ergebnismatrix_1, digits = 3), round(Ergebnismatrix_2, digits = 3), round(Ergebnismatrix_3, digits = 3))

  dimnames(Ergebnis) = list( c("0.00", "0.01", "0.1", "0.2", "0.3", "0.4", "0.5"),
                             c("p=3 1", "2","3","4","p<=0.05","p<=0.05","p<=0.05","p<=0.05", "p=5 1", "2","3","4","p<=0.05","p<=0.05","p<=0.05", "p<=0.05", "p=10 1", "2","3","4","p<=0.05","p<=0.05","p<=0.05", "p<=0.05"))
  cat("POC 4 Stages using BCDWIT(2/3,b) for b=3/5/10 - applying the randomization procedure to each stage \n")
  cat("N: ", n)
  cat(" ")
  cat("NoReps:", norep)
  cat("\n")
  print(Ergebnis)
}

