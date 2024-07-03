# script to test the multitaper estimation functions for accuracy and
# plotting capability
nrep <- 30; n <- 1000
X <- sapply(seq_len(nrep), function(x) {
  arima.sim(list(ar = runif(1, min = 0.6, max = 0.7)), n = n)
  })
mtspec <- clustFBE::sine_multitaper(X)
