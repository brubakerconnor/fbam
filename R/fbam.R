#' Frequency band analysis of multiple time series
#'
#' @param X Matrix of shape (`n`, `nrep`) where `n` is the length of the time
#' series observations and `nrep` is the number of independent replicates.
#' @param nclust_grid Vector of possible number of clusters. Duplicate values
#' are removed.
#' @param nbands_grid Vector of possible number of bands. Duplicate values are
#' removed.
#' @param parallel Boolean. Should parallelization be used? Defaults to `TRUE`.
#' If a cluster has already been created, that cluster will be used.
#' @param popsize Integer. Size of the population carried through each
#' generation of the GA.
#' @param nparents Integer. The number of parents to select that will be
#' subsequently subjected to crossover and/or mutation in each generation. Must
#' be at least 2. Defaults to twice the population size.
#' @param nelite Integer. The number of elites (most fit) solutions to carry
#' over from the previous generation into the current generation at each
#' iteration of the GA. Default is `2`. This can be set to `0` to prevent
#' the implementation of an elitism policy.
#' @param pcrossover Probability that any two selected parents will be crossed.
#' Setting this to `0` means no crossover will be performed. Default is `0.5`.
#' @param pmutate Probability that any endpoint or collapsed measure of a single
#' solution will be subjected to mutation (independent of the other values and
#' other solutions). Default is `0.25`.
#' @param endpoint_range Integer or vector of integers. The total maximum width
#' of the window in which any new endpoint is selected. See details below for
#' default values.
#' @param collapsed_range Double or vector of double. the total width of the
#' window in which any new collapsed measure is shifted. See details below for
#' default values.
#' @param maxgen Integer. The maximum number of generations the GA is allowed
#' to run. Defaults to `500`.
#' @param maxrun Integer. The maxumum number of generations the GA is allowed
#' to run without significant improvement in the maximum fitness of the
#' population. Defaults to `50`.
#' @param tol Double. The percent difference above which a significant
#' improvement in maximum fitness is declared. Larger values mean that larger
#' improvements are needed to increase the run count and therefore lead to
#' shorter convergence times. Defaults to `1e-3` (0.1%).
#' @param s the size of the window to shift each evenly spaced frequency within
#' when initializing endpoints for each solution. Defaults to 10% of the number
#' of Fourier frequencies.
#' @param ntapers Integer. The number of tapers to use in estimation of the
#' replicate spectra from `X`. Defaults to the floor of the square root of
#' length of the time series in `X`.
#' @param verbose Boolean. Whether to print progress of the algorithm to the
#' console. Defaults to `FALSE`.
#'
#' @return A list with the following components.
#' \item{\code{X}}{Input matrix of time series data}
#' \item{\code{spec}}{Matrix of shape `(nfreq, nrep)` of spectral estimates used
#' in optimization.}
#' \item{\code{labels}}{Clustering labels}
#' \item{\code{endpoints_index}}{A matrix of shape `(nclust, nbands - 1)`
#' containing the Fourier indices of the frequency bands associated to each
#' cluster.}
#' \item{\code{endpoints}}{A matrix of shape `(nclust, nbands - 1)`
#' containing the Fourier frequencies of the frequency bands associated to each
#' cluster.}
#' \item{\code{endpoints}}{A matrix of shape `(nclust, nbands)`
#' containing the estimated collapsed measures using the frequency bands
#' associated to each cluster.}
#' \item{\code{loss}}{The value of the optimized loss function.}
#' \item{\code{avgfit}}{A history of the average fitness value in each
#' generation over the course of evolution.}
#' \item{\code{maxfit}}{A history of the maximum fitness value in each
#' generation over the course of evolution.}
#' \item{\code{ninfeasible}}{A history of the number of infeasible solutions in
#' each generation over the course of evolution.}
#' \item{\code{params}}{A list of the input parameters used in the genetic
#' algorithm(s).}
#' @export
#'
#' @examples
#' data(model1)
#' # search over a grid of 2 to 4 clusters, and 2 to 6 bands
#' out <- clustFBE(model1$x, 2:4, 2:6, parallel = T, maxrun = 10)
#'
#' # fix clusters but vary bands
#' out <- clustFBE(model1$x, 3, 2:6, parallel = T, maxrun = 10)
#'
#' # fix both
#' out <- clustFBE(model1$x, 3, 3, parallel = T, maxrun = 10)
#'
#' # no clustering
#' out <- clustFBE(model1$x, 1, 3, parallel = T, maxrun = 10)
fbam <- function(X,
                 nclust_grid,
                 nbands_grid,
                 parallel = TRUE,
                 popsize = 50,
                 nparents = 2 * popsize,
                 nelite = 2,
                 pcrossover = 0.5,
                 pmutate = 0.25,
                 endpoint_range = NULL,
                 collapsed_range = NULL,
                 maxgen = 500,
                 maxrun = 50,
                 tol = 1e-3,
                 s = NULL,
                 ntapers = floor(sqrt(dim(X)[1])),
                 verbose = T) {
  # enable parallelization, if specified
  if (parallel) {
    if (!foreach::getDoParRegistered()) {
      numCores <- parallel::detectCores() - 2
      cat("Registering a cluster with ", numCores, " threads.\n")
      cl <- parallel::makeCluster(numCores, type = "FORK")
      doParallel::registerDoParallel(cl)
      parallel::clusterCall(cl, function() {
        library(foreach); library(doParallel)
      })
    } else {
      message("Registered cluster detected; using ", foreach::getDoParName(),
              " with ", foreach::getDoParWorkers(), " workers.")
    }
  } else {
    message("Parallelization not specified. Optimizing sequentially.")
  }

  nclust_grid <- unique(nclust_grid)
  nbands_grid <- unique(nbands_grid)
  nclust_grid_size <- length(nclust_grid)
  nbands_grid_size <- length(nbands_grid)

  if (is.vector(X)) X <- as.matrix(X)
  nfreq <- dim(X)[1]; nrep <- dim(X)[2]
  if (any(nclust_grid < 1) | any(nclust_grid > nrep)) {
    stop("Ensure valid nclust grid.")
  }
  if (any(nbands_grid < 1) | any(nbands_grid > nfreq - 1)) {
    stop("Ensure valid nbands grid.")
  }

  # optimize loss over grid
  `%:%` <- foreach::`%:%`
  `%dopar%` <- foreach::`%dopar%`
  grid <- foreach::foreach(j = 1:nclust_grid_size) %:%
    foreach::foreach(b = 1:nbands_grid_size) %dopar% {
      return(ga(
        X = X, nclust = nclust_grid[j], nbands = nbands_grid[b],
        popsize = popsize, pop = NULL, nparents = nparents, nelite = nelite,
        pcrossover = pcrossover, pmutate = pmutate,
        endpoint_range = endpoint_range, collapsed_range = collapsed_range,
        maxgen = maxgen, maxrun = maxrun, tol = tol, s = s, ntapers = ntapers,
        verbose = FALSE)
      )
    }

  # compute selection index for each solution to find "best" solution
  if (nclust_grid_size + nbands_grid_size > 2) {
    current_best <- grid[[1]][[1]]
    current_best_si <- selection_index(current_best$spec,
                                       current_best$labels,
                                       current_best$endpoints_index,
                                       current_best$collapsed)
    for (j in 1:nclust_grid_size) {
      for (b in 1:nbands_grid_size) {
        current_out <- grid[[j]][[b]]
        current_out_si <- selection_index(current_out$spec,
                                          current_out$labels,
                                          current_out$endpoints_index,
                                          current_out$collapsed)
        if (current_out_si < current_best_si) {
          current_best <- current_out; current_best_si <- current_out_si
        }
      }
    }
    return(list(selected = current_best, grid = grid))
  }
  return(grid[[1]][[1]])
}
