library(testthat)
library(Rcpp)

source("R/simulation.R")
sourceCpp("src/mcmc.cpp")

test_that("normal approximation transition probability cpp code is working", {
  x_0 <- .4
  x_t <- .55
  n_gen <- 200
  s <- .1
  h <- 1.0
  N <- 5000
  trans_prob_cpp <- norm_trans_prob(x_t, x_0, n_gen, s, h, N)

  var <- x_0 * (1 - x_0)
  mu <- x_0 + (2 * s) * var * (x_0 + (h * (1 - (2 * x_0))))
  sigma <- sqrt( var * (n_gen / (2 * N)) )
  eps <- .001
  trans_prob_R <- pnorm(x_t + eps, mu, sigma, 1, 0) - pnorm(x_t - eps, mu, sigma, 1, 0);
  print(paste0(trans_prob_cpp, " ", trans_prob_R))
  expect_equal(trans_prob_cpp, trans_prob_R, tolerance=1e-3)
})


test_that("foward alg. likelihood is maximized by known parameter", {
  x0 <- .4
  s <- .1
  h <- 1.0
  N <- 5000
  nChrs <- c(1000, 1000, 1000, 1000)
  gens <- c(1, 25, 75, 100)
  df <- get_wf_samples(x0, N, s, nChrs, gens)
  #O <- as.matrix(df)
  O <- cbind(c(387, 673, 853, 888), c(1000, 1000, 1000, 1000),  c(1, 25, 75, 100))
  nObs <- nrow(O)
  states <- seq(0.0, 1.0, .05)
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