library(Rcpp)
library(ggplot2)
library(dplyr)
source("R/simulation.R")
sourceCpp("~/git/hgen48600_project/mcmc.cpp")

x0 <- .4
s <- .01
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

s0 <- 0.0
propSd <- .001
nIter <- 10000

posteriorSamples <- mcmc(O, states, nStates, s0, h, N, propSd, nIter) 

plot(posteriorSamples,type="l")
#hist(posteriorSamples[20000:nIter])
