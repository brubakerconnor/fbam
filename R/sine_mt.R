#' Sine multitaper estimation
#'
#' Sine multitaper estimation for one or more replicate time series.
#'
#' @importFrom stats fft
#' @param x Vector or matrix (columns) of replicate time series. A vector is treated as a single replicate.
#' @param ntapers Number of sine tapers. Default is the floor of the square root of the length of the time series (number of rows in \code{x}).
#'
#' @return A list with the following components:
#' \describe{
#'  \item{mtfreq}{Vector of Fourier frequencies}
#'  \item{mtspec}{Matrix of spectral estimates (columns)}
#' }
#' @details
#' The zero frequency and Nyquist frequency (folding frequency) are not included in the resulting estimates.
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
#' x <- arima.sim(list(ar = 0.6), n = 128)
#' sine_mt(x)
#'
#' # many replicates
#' x <- matrix(nrow = 128, ncol = 8)
#' for (i in 1:ncol(x)) x[,i] <- arima.sim(list(ar = runif(1, 0.2, 0.8)), n = 128)
#' sine_mt(x)
sine_mt <- function(x, ntapers = NULL) {

  # input checks ------------------------------------------------------------

  x <- as.matrix(x)
  if (is.null(ntapers)) ntapers <- floor(sqrt(nrow(x)))

  len <- dim(x)[1]
  if (dim(x)[1] < dim(x)[2]) {
    warning("Wide matrix was provided. Double check that the matrix is such
            that columns (rather than rows) correspond to replicates.")
  }


  # estimation --------------------------------------------------------------

  mtfreq <- (1 / len) * 1:(floor(len / 2) - 1)
  nfreq <- length(mtfreq)
  mtspec <- apply(x, 2, function(z) {
    sine_tapers <- outer(seq_len(ntapers), seq_len(len), function(z, y) {
      sqrt(2 / (len + 1)) * sin(pi * y * z / (len + 1))
    })
    ds <- apply(sine_tapers, 1, function(r) Mod(fft(r * z)[2:(nfreq + 1)])^2)
    return(rowMeans(ds))
  })

  return(list(mtfreq = mtfreq, mtspec = mtspec))
}
