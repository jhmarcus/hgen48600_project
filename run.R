library(Rcpp)
library(ggplot2)
library(dplyr)
source("R/simulation.R")
sourceCpp("src/mcmc.cpp")

if(F){
x0 <- .4
s <- .1
h <- 1.0
N <- 5000
#nChrs <- c(1000, 1000, 1000, 1000, 1000, 1000)
#gens <- c(1, 75, 100, 150, 200, 250)
nChrs <- rep(1000, 10)
gens <- seq(1, 250, 25)
df <- get_wf_samples(x0, N, s, nChrs, gens)
O <- as.matrix(df)
print(O)

states <- seq(0.0, 1.0, .025)
nStates <- length(states)
}


sourceCpp("src/mcmc.cpp")
x0 <- .5
s <- .1
h <- 1.0
N <- 5000
n_chrs <- rep(10000, 20)
gens <- seq(1, 200, 10)
df <- get_wf_samples(x0, N, s, n_chrs, gens)
O <- as.matrix(df)
#O <- cbind(c(387, 673, 853, 888), c(1000, 1000, 1000, 1000),  c(1, 25, 75, 100))
n_obs <- nrow(O)
states <- seq(0.0, 1.0, .05)

s_0 <- 0.0
prop_sd <- .005
n_iter <- 50000

posterior_samples <- mcmc(O, states, s_0, h, N, prop_sd, n_iter) 
plot(posterior_samples,type="l")
hist(posterior_samples[5000:n_iter]/4)
