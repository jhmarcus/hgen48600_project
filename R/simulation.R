library(dplyr)

simulate_wf_selection <- function(x_0, n_gen, s, h, N){
  x <- c(x_0)
  for(i in 1:n_gen){
    # p <- x0 + (s * x0 * (1 - x0) * (h + (1 - 2 * h) * x0))
    var <- x_0 * (1 - x_0)
    p <- x_0 + (2 * s) * var * (x_0 + (h * (1 - (2 * x_0))))
    x_0 <- rbinom(1, size = 2*N, prob = p) / (2 * N)
    x <- c(x, x_0)
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
