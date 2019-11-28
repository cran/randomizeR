#######################################################################
# --------------------------------------------------------------------#
# Tests for objects from the endpoint class and associated functions  #
# --------------------------------------------------------------------#
#######################################################################
context("Endpoints")

# Normally distributed endpoints
test_that("normEnd returns valid object", {
	mu <-rnorm(2); sigma <- sample(1:10,2, replace=T)
    expect_is(normEndp(mu,sigma), "endpoint")
    
  # Check with negative standard deviation  
	sigma <- -sample(1:10,2, replace=T)
	expect_error(normEndp(mu,sigma))
	# Check with negative standard deviation and wrong parameter lengths
	mu <-rnorm(3)
	expect_error(normEndp(mu,sigma))
	# Check with wrong parameter lengths
	mu <-rnorm(3); sigma <- sample(1:10,2, replace=T)
	expect_error(normEndp(mu,sigma))
	# Check with wrong parameter lengths
	mu <-rnorm(2); sigma <- sample(1:10,4, replace=T)
	expect_error(normEndp(mu,sigma))
  }
)


test_that("expEnd returns valid object", {
  # Check with expected lambda, cenRate, accrualTime, cenTime, 
  lambda <- abs(rnorm(2)); cenRate <- abs(rnorm(1))
  accrualTime <- sample(1:200, 1); cenTime <- accrualTime + sample(1:200, 1)
  expect_is(expEndp(lambda = lambda, cenRate = cenRate, accrualTime = accrualTime, cenTime = cenTime), "endpoint")
  
  # Check with vectors, not numericals
  cenRate <- abs(rnorm(2)); accrualTime = abs(rnorm(2)); cenTime = abs(rnorm(2))
  expect_error(expEndp(lambda = lambda, cenRate = cenRate))
  expect_error(expEndp(lambda = lambda, accrualTime = accrualTime))
  expect_error(expEndp(lambda = lambda, cenTime = cenTime))
  # Check with negative all negative separately and expect error
  cenRate <- abs(rnorm(1)); accrualTime <- sample(1:200, 1); cenTime <- accrualTime + sample(1:200, 1)
  expect_error(expEndp(lambda = -lambda, cenTime = cenTime))
  expect_error(expEndp(lambda = lambda, cenRate = -cenRate))
  expect_error(expEndp(lambda = lambda, accrualTime = -accrualTime))
  expect_error(expEndp(lambda = lambda, cenTime = -cenTime))
  # Check with centime smaller than accrualtime and expect error
  expect_error(expEndp(lambda = lambda, accrualTime = cenTime, cenTime = cenTime))
  expect_error(expEndp(lambda = lambda, accrualTime = cenTime+1, cenTime = cenTime))

  }
)
