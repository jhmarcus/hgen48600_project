g <- function(x, s, h){
  a <- s * x * (1 - x) 
  b <- h + ((1 - 2 * h) * x)
  return(x + (a * b))
}

moments_delta_method <- function(p0, t, s, h, N){
  mu <- p0
  sigma2 <- (1 / (2*N)) * p0 * (1 - p0)
  for(i in 1:t){
    a <- (1 / ( 2 * N)) * g(mu, s, h) * (1 - g(mu, s, h))
    b <- ((1 + s * (2 - (3 * mu)) * mu)^2) * sigma2
    sigma2 <- a + b
    mu <- g(mu, s, h)
  }
  return(c(mu, sigma2))
}


moments_delta_method(.1, 100, 0, .5, 50)