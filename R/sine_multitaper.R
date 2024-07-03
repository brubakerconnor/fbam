#' Multitaper estimation
#'
#' Multitaper estimation using sine tapers for one or more mutually independent
#' and stationary time series.
#'
#' @param X Vector of length `n` or matrix of shape (`n`, `nrep`) where `nrep`
#' is the number of independent time series. If a vector is provided, it is
#' treated as one series of observations from the same process.
#' @param ntapers Integer. The number of tapers to use in estimation. Defaults
#' to the floor of the square root of the number of observations `n`.
#'
#' @return A matrix of shape (`floor(n / 2) - 1`, `nrep`) even if `nrep = 1`.
#' @export
#'
#' @examples
#' # single replicate
#' x <- arima.sim(model(ar = 0.5), n = 1000)
#' spec <- TSA::ARMAspec(model(ar = 0.5))
#' mtspec <- sine_multitaper(x)
sine_multitaper <- function(X, ntapers = floor(sqrt(dim(as.matrix(X))[1]))) {
  if (is.vector(X)) X <- as.matrix(X)
  n <- dim(X)[1]
  if (n < dim(X)[2]) {
    warning("Wide matrix was provided. Ensure that the matrix is formatted such
            that columns (rather than rows) correspond to individual
            realizations.")
  }

  mtfreq <- seq(from = 1 / n, by = 1 / n, length = floor(n / 2) - 1)
  mtspec <- apply(X, 2, function(y) sine_multitaper_single(y, ntapers))

  return(list(
    mtfreq = mtfreq,
    mtspec = mtspec
  ))
}

sine_multitaper_single <- function(x, ntapers = floor(sqrt(length(x)))) {
  n <- length(x) # number of observations
  m <- floor(n / 2) - 1 # number of Fourier frequencies excluding 0 and n / 2

  sine_tapers <- outer(seq_len(ntapers), seq_len(n), function(x, y) {
    sqrt(2 / (n + 1)) * sin(pi * y * x / (n + 1))
  })
  direct_spec <- apply(sine_tapers, 1, function(r) Mod(fft(r * x)[2:(m + 1)])^2)
  return(rowMeans(direct_spec))
}


