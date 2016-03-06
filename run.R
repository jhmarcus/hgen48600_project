library(Rcpp)
library(ggplot2)
library(dplyr)
source("R/simulation.R")
sourceCpp("src/mcmc.cpp")

x0 <- .3
s <- .004
h <- .5
N <- 5000
#n_chrs <- rep(10000, 20)
#gens <- seq(1, 200, 10)
n_chrs <- c(1000, 1000, 1000, 1000, 1000, 1000)
gens <- c(1, 75, 100, 150, 200, 250)
df <- get_wf_samples(x0, N, s, n_chrs, gens)
O <- as.matrix(df)
#O <- cbind(c(387, 673, 853, 888), c(1000, 1000, 1000, 1000),  c(1, 25, 75, 100))
n_obs <- nrow(O)
states <- seq(0.0, 1.0, .025)

s_0 <- 0.001
prop_sd <- .01
n_iter <- 10000

posterior_samples <- mcmc(O, states, s_0, h, N, prop_sd, n_iter) 
plot(posterior_samples,type="l")
hist(posterior_samples[4000:length(posterior_samples)])
