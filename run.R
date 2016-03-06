library(Rcpp)
library(ggplot2)
library(dplyr)
source("R/simulation.R")
sourceCpp("src/mcmc.cpp")

x0 <- .5
s <- .1
h <- .5
N <- 5000
n_chrs <- rep(10000, 20)
gens <- seq(1, 200, 10)
#n_chrs <- c(1000, 1000, 1000, 1000, 1000, 1000)
#gens <- c(1, 75, 100, 150, 200, 250)
df <- get_wf_samples(x0, N, s, n_chrs, gens)
O <- as.matrix(df)
n_obs <- nrow(O)
states <- seq(0.0, 1.0, .025)

s_0 <- 0.01
prop_sd <- .005
n_iter <- 20000

posterior_samples <- mcmc(O, states, s_0, h, N, prop_sd, n_iter) 
plot(posterior_samples,type="l")
hist(posterior_samples)
#hist(posterior_samples[3000:length(posterior_samples)])
