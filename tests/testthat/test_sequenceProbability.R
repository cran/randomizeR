###############################################
# --------------------------------------------#
# Tests for the probabilities of sequences    #
# --------------------------------------------#
###############################################

context("Sequence Probability")

test_that("total probability adds up to 1", {
  
  N <- sample(seq(2, 12, 2), 1)                 # Sample number of patients
  mti <- sample(N/2, 1)                         # Sample maximum tolerated imbalance
  p <- sample(seq(0.5001, 1, 0.05), 1)          # biased coin parameter
  nr <- sample(3,1)                             # sample number of blocks 
  blocks <- sample(seq(2, 10, 2), nr)           # sample blocks
  
  
  allCRseq <- getAllSeq(crPar(N)) 
  allPBRseq <- getAllSeq(pbrPar(bc = blocks)) 
  allMPseq <- getAllSeq(mpPar(N, mti)) 
  allBSDseq <- getAllSeq(bsdPar(N, mti)) 
  allEBCseq <- getAllSeq(ebcPar(N, p)) 
  allTBDseq <- getAllSeq(tbdPar(bc = blocks)) 
  allHADAseq <- getAllSeq(hadaPar(N = 10))
  allRARseq <- getAllSeq(rarPar(N))
  
  expect_equal(sum(getProb(allCRseq)), 1)
  expect_equal(sum(getProb(allPBRseq)), 1)
  expect_equal(sum(getProb(allMPseq)), 1)
  expect_equal(sum(getProb(allBSDseq)), 1)
  expect_equal(sum(getProb(allEBCseq)), 1)
  expect_equal(sum(getProb(allTBDseq)), 1)
  expect_equal(sum(getProb(allHADAseq)), 1) 
  expect_equal(sum(getProb(allRARseq)), 1)
})
