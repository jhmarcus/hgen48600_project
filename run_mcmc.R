library(Rcpp)
library(ggplot2)
library(dplyr)
sourceCpp("mcmc.cpp")

simulate_wf <- function(x0, nGen, s, h, N){
  x <- c(x0)
  for(i in 1:nGen){
    p <- x0 + (s * x0 * (1 - x0) * (h + (1 - 2 * h) * x0))
    x0 <- rbinom(1, size = 2*N, prob = p) / (2 * N)
    x <- c(x, x0)
  }
  return(x)
}

sample_wf <- function(x0, N, s, nChrs, gens){
  x <- simulate_wf(x0, tail(gens, n=1), s, 1, 2*N)
  df <- data.frame(x = x[gens], gen = gens, nChr = nChrs)
  df <- df %>% mutate(alleleCount = rbinom(n(), nChr, x)) %>%
        select(alleleCount, nChr, gen)
  return(df)
}

x0 <- .4
s <- .1
h <- 1.0
N <- 5000
nChrs <- c(1000, 1000, 1000, 1000, 1000, 1000)
gens <- c(1, 75, 100, 150, 200, 250)
df <- sample_wf(x0, N, s, nChrs, gens)
O <- as.matrix(df)
print(O)

states <- seq(0.0, 1.0, .025)
nStates <- length(states)

s0 <- 0.0
propSd <- .005
nIter <- 10000

posteriorSamples <- mcmc(O, states, nStates, s0, h, N, propSd, nIter) 

plot(posteriorSamples,type="l")
#hist(posteriorSamples[20000:nIter])
