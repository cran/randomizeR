###################################################################
# ----------------------------------------------------------------#
# Tests for the analyse function                                  #
# ----------------------------------------------------------------#
###################################################################


context("Analysis of Stratified Randomization Sequences")

test_that("Function Accepts only valid inputs ",{
  N <- 10
  endp <- normEndp(c(0,0),c(1,1))
  allocRatio <- c(1,3)
  strata <- 2
  theta <- 0
  eta = 0.05
  pr = c("CR")

  expect_error(analyse(N))
  expect_error(analyse(N,endp ,strata, theta, eta, pr))
  expect_error(analyse(N, allocRatio,strata, theta, eta, pr))


})
test_that("Function Accepts only valid randomization procedures ",{
  N <- 10
  endp <- normEndp(c(0,0),c(1,1))
  allocRatio <- c(1,3)
  strata <- 2
  theta <- 0
  eta = 0.05
  pr = c("PBR(3)")
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))
  N <- 73
  pr = c("PBR(4)")
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))


})


test_that("Allocation Ratio and Strata should be compatible ",{
  N <- 10
  endp <- normEndp(c(0,0),c(1,1))
  allocRatio <- c(1,1,3)
  strata <- 2
  theta <- 0
  eta = 0.05
  pr = c("CR")
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))
  strata <- 5
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))


})



test_that("Function Accepts only valid randomization procedures ",{
  N <- 10
  endp <- normEndp(c(0,0),c(1,1))
  allocRatio <- c(1,3)
  strata <- 2
  theta <- 0
  eta = 0.05
  pr = c("PBR")
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))
  pr = c("PBR(not_a_numerical)")
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))
  pr = c("RTBD(not_a_logical)")
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))
  pr = c("cr","bsd(8)")
  expect_error(analyse(N, endp, allocRatio,strata, theta, eta, pr))


})
