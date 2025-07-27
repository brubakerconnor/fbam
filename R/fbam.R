#' Frequency band analysis for multiple stationary time series
#'
#' @param x Matrix of replicate time series (columns). A vector is treated as a single replicate.
#' @param nbands Integer or vector of numbers of frequency bands to consider.
#' @param nsubpop Integer or vector of numbers of subpopulations to consider.
#' @param ntapers Number of tapers to use in sine multitaper estimation.
#' @param ncores Number of cores to use for parallelization.
#' @param ... Additional arguments to [genetic_algorithm].
#'
#' @returns A list with the following components:
#' \describe{
#'  \item{grid}{All solutions across the grid of parameters 'nsubpop' > 1 and 'nbands'.}
#'  \item{solution}{Selected solution from 'grid' using validation criteria.}
#'  \item{grid1}{All solutions across the grid of parameters 'nsubpop' = 1 and 'nbands'.}
#'  \item{solution1}{Selection solution from 'grid1' using validation criteria.}
#' }
#' @export
#'
#' @examples
#' x <- matrix(nrow = 128, ncol = 12)
#' for (i in 1:ncol(x)) x[,i] <- arima.sim(list(ar = runif(1, 0.2, 0.8)), n = 128)
#' fbam(x, nbands = 2:5, nsubpop = 1)
#' fbam(x, nbands = 2:5, nsubpop = 3)
#' fbam(x, nbands = 2, nsubpop = 1:5)
#' fbam(x, nbands = 2:5, nsubpop = 1:5)
fbam <- function(x, nbands = 2, nsubpop = 1, ntapers = NULL, ncores = 1, ...) {


  # input checks ------------------------------------------------------------

  if (is.vector(x)) x <- matrix(x)

  nbands <- nbands[nbands >= 2]
  if (length(nbands) == 0) {
    warning("Invalid input for 'nbands'. Setting to 2.")
    nbands <- 2
  }

  nsubpop <- nsubpop[nsubpop >= 1]
  if (length(nsubpop) == 0) {
    warning("Invalid input for 'nsubpop'. Setting to 1.")
    nsubpop <- 1
  }


  # spectrum estimation -----------------------------------------------------

  sine_mt_out <- sine_mt(x, ntapers)
  freq <- sine_mt_out$mtfreq
  spec <- sine_mt_out$mtspec

  # grid search -------------------------------------------------------------

  # if nsubpop = 1 is one of the requested values:
  # cannot compare solutions with 1 subpop and those with more than 1
  # grid search over varying nbands with nsubpop = 1 fixed done first
  if (1 %in% nsubpop) {
    nsubpop1_grid <- parallel::mclapply(nbands, function(L) {
      genetic_algorithm(spec, nbands = L, nsubpop = 1, ...)
    }, mc.cores = ncores)
    nsubpop1_validation <- unlist(lapply(nsubpop1_grid, function(x) x$validation))
    nsubpop1_selected <- nsubpop1_grid[[which.min(nsubpop1_validation)]]
  } else {
    nsubpop1_grid <- NULL
    nsubpop1_selected <- NULL
  }


  # if nsubpop contains values greater than 1:
  # can compare solutions corresponding to any value of nsubpop > 1
  # grid search over varying nbands AND nsubpop as requested by user
  if (length(nsubpop[nsubpop > 1]) > 0) {
    param_grid <- expand.grid(nbands = nbands, nsubpop = nsubpop[nsubpop > 1])
    grid <- parallel::mclapply(1:nrow(param_grid), function(i) {
      nsubpop <- param_grid$nsubpop[i]
      nbands <- param_grid$nbands[i]
      genetic_algorithm(spec, nbands = nbands, nsubpop = nsubpop, ...)
    }, mc.cores = ncores)
    validation <- unlist(lapply(grid, function(x) x$validation))
    solution <- grid[[which.min(validation)]]
  } else {
    grid <- NULL
    solution <- NULL
  }

  return(list(grid = grid, solution = solution,
              grid1 = nsubpop1_grid, solution1 = nsubpop1_selected))
}
