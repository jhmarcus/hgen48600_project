library(dplyr)

simulate_wf_selection <- function(x0, nGen, s, h, N){
  x <- c(x0)
  for(i in 1:nGen){
    p <- x0 + (s * x0 * (1 - x0) * (h + (1 - 2 * h) * x0))
    x0 <- rbinom(1, size = 2*N, prob = p) / (2 * N)
    x <- c(x, x0)
  }
  return(x)
}

get_wf_samples <- function(x0, N, s, nChrs, gens){
  x <- simulate_wf_selection(x0, tail(gens, n=1), s, 1, 2*N)
  df <- data.frame(x = x[gens], gen = gens, nChr = nChrs)
  df <- df %>% mutate(alleleCount = rbinom(n(), nChr, x)) %>%
      	select(alleleCount, nChr, gen)
  return(df)
}
