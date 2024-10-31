#' @importFrom stats runif
select_parents <- function(pop, nparents, fitvals) {
  # ranks: highest fitness gets largest rank = popsize - 1
  ranks <- rank(fitvals) - 1
  probs <- 1 - exp(-ranks) / sum(1 - exp(-ranks))  # exponential ranking
  cdf <- cumsum(probs)  # cumulative probability distribution for SUS

  # stochastic universal sampling (SUS) algorithm
  cr <- i <- 1
  r <- runif(1, min = 0, max = 1 / nparents); s <- vector("integer", nparents)
  while (cr <= nparents) {
    while (r <= cdf[i]) {
      s[cr] <- i; r <- r + (1 / nparents); cr <- cr + 1
    }
    i <- i + 1
  }

  # return parents
  return(pop[s,])
}
