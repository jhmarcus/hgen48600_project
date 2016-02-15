library(testthat)
library(Rcpp)

source("../simulation.R")
sourceCpp("../mcmc.cpp")

test_that("foward alg. likelihood is maximized by known parameter", {
  x0 <- .4
  s <- .1
  h <- 1.0
  N <- 5000
  nChrs <- c(1000, 1000, 1000, 1000, 1000, 1000)
  gens <- c(1, 75, 100, 150, 200, 250)
  df <- sample_wf(x0, N, s, nChrs, gens)
  O <- as.matrix(df)
  nObs <- nrow(O)
  states <- seq(0.0, 1.0, .025)
  nStates <- length(states)
  sCoefs <- c(s, s+seq(.05, .2, .05), s-seq(.05, .2, .05))
  likelihoods <- c()
  for (sCoef in sCoefs){
    likelihood <- sum(foward(O, states, nStates, sCoef, h, N)[,nObs])
    likelihoods <- c(likelihoods, likelihood)
  }
  print(likelihoods)
  expect_that(max(likelihoods), equals(likelihoods[1]))
})