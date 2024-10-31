#' Sine multitaper estimation
#'
#' Compute the multitaper estimator for one or more replicate time series
#' using the orthogonal sine tapers.
#'
#' @importFrom stats fft
#' @param X Column-wise matrix of replicate time series. A vector is treated as a single replicate.
#' @param ntapers Number of tapers to use in constructing the multitaper estimator. Default value is
#' the floor of the square root of the length of the time series (number of rows in \code{X}).
#'
#' @return A list with the following components:
#' \itemize{
#'  \item \code{mtfreq}: Vector of frequencies
#'  \item \code{mtspec}: Column-wise matrix of spectral estimates in correspondence to the rows of \code{X}.
#' }
#' @details
#' The zero and one-half frequencies are not included in the output estimates.
#' The choice of \code{ntapers = floor(sqrt(nrow(X)))} ensures consistency of the resulting estimator.
#'
#' @references
#' Riedel, K.S., and A. Sidorenko. 1995. *Minimum Bias Multiple Taper Spectral Estimation.* IEEE Transactions on Signal Processing 43 (1): 188–95.
#' Available at \url{https://doi.org/10.1109/78.365298}.
#'
#' Walden, A. T. 2000. *A Unified View of Multitaper Multivariate Spectral Estimation.* Biometrika 87 (4): 767–88.
#' Available at \url{https://doi.org/10.1093/biomet/87.4.767}.
#'
#' Percival, Donald B, and Andrew T Walden. 2020. *Spectral Analysis for Univariate Time Series*. Vol. 51. Cambridge University Press.
#' @export
#'
#' @examples
#' # single replicate
#' x <- arima.sim(list(ar = 0.6), n = 500)
#' sine_mt(x)
#'
#' # many replicates
#' X <- matrix(nrow = 500, ncol = 20)
#' for (i in 1:20) X[,i] <- arima.sim(list(ar = runif(1, 0.2, 0.8)), n = 500)
#' sine_mt(X)
sine_mt <- function(X, ntapers = floor(sqrt(dim(as.matrix(X))[1]))) {
  if (!is.matrix(X)) X <- as.matrix(X) # vector treated as single replicate
  len <- dim(X)[1]
  if (len < dim(X)[2]) {
    warning("Wide matrix was provided. Ensure that the matrix is formatted such
            that columns (rather than rows) correspond to replicates.")
  }
  mtfreq <- seq(from = 1 / len, by = 1 / len, length = floor(len / 2) - 1)
  mtspec <- apply(X, 2, function(x) {
    m <- floor(len / 2) - 1 # number of frequencies excluding 0 and len / 2
    sine_tapers <- outer(seq_len(ntapers), seq_len(len), function(x, y) {
      sqrt(2 / (len + 1)) * sin(pi * y * x / (len + 1))
    })
    ds <- apply(sine_tapers, 1, function(r) Mod(fft(r * x)[2:(m + 1)])^2)
    return(rowMeans(ds))
  })
  return(list(mtfreq = mtfreq, mtspec = mtspec))
}
