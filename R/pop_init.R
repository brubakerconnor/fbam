#' Population initialization
#'
#' Initialize a population for evolution in a genetic algorithm using
#' spectral estimates to guide initialization.
#'
#' @param popsize Integer. Number of solutions in the population.
#' @param nclust Integer. Number of clusters in each solution.
#' @param nbands Integer. Number of frequency bands associated to each cluster.
#' @param spec Vector of length `nfreq` in the case of one replicate or matrix
#' of shape (`nfreq`, `nrep`) in the case of multiple replicates where `nfreq`
#' is the number of Fourier frequencies from estimation (typically
#' `floor(T / 2) - 1`).
#' @param s the size of the window to shift each evenly spaced frequency within
#' when initializing endpoints for each solution. Defaults to 10% of the number
#' of Fourier frequencies.
#'
#' @details
#' Population initialization proceeds in the following way. First, an evenly
#' spaced grid of `nbands + 1` boundaries is set down (this includes the first
#' index `1` and the last index `nfreq + 1`). For each solution, these evenly
#' spaced endpoints are randomly shifted around within a window of width `s`,
#' sorted, and used as the initial endpoints for each cluster. To initialize
#' cluster-average collapsed measures for each cluster, a single replicate is
#' sampled and collapsed according to these endpoints.
#'
#' Using the same endpoints for each cluster initially was observed to result
#' in fewer infeasible solutions, that is, solutions that resulted in empty
#' clusters when replicates were assigned based on minimum L2 distance.
#'
#'
#' @return A list of length `popsize` where each entry in the list is a matrix
#' of shape (`nclust`, `2 * nbands - 1`). Each row of these matrices represents
#' the cluster center of a single cluster. These rows are organized such that
#' the first `nbands - 1` entries are the indices of the frequency band
#' boundaries and the remaining entries are the corresponding cluster-average
#' collapsed measures.
#'
#' @examples
#' # 10 replicates of white noise where T = 400
#' # spec is simulated using the asymptotic distribution of the periodogram
#' # for a white noise process with variance 1
#' length <- 400; nfreq <- floor(length / 2) - 1
#' nrep <- 10
#' spec <- matrix(rchisq(nfreq * nrep, df = 2) / 2, nrow = nfreq)
#' pop <- pop_init(popsize = 100, nclust = 2, nbands = 3, spec = spec)
#' print(pop[[1]])
pop_init <- function(popsize, nclust, nbands, spec,
                     s = floor(0.1 * dim(spec)[1])) {
  nfreq <- dim(spec)[1]
  nrep <- dim(spec)[2]
  nendpoints <- nbands - 1

  # checks on input parameters
  if (popsize < 1) {
    stop(paste(c("Requested population size not valid. Got", popsize, "."),
               collapse = " "))
  }
  if (nclust < 1) {
    stop(paste(c("Requested number of clusters is not positive. Got", nclust,
                 "."), collapse = " "))
  }
  if (nclust > nrep) {
    stop(paste(c("Requested number of clusters is too large. Got", nclust,
                 ". Must be no greater than the number of replicates in spec."),
               collapse = " "))
  }
  if (nclust > nrep / 2) {
    warning("Requested number of clusters might be too large (exceeds half the
            number of replicates).")
  }
  if (nbands < 2) {
    stop(paste(c("Requested number of bands is too small. Got", nbands,
                 ". Must be at least 2."),
               collapse = " "))
  }
  if (nbands > nfreq) {
    stop(paste(c("Requested number of bands is too large Got", nbands,
                 ". Must be no greater than number of Fourier frequencies."),
               collapse = " "))
  }
  if (nbands > nfreq / 2) {
    warning("Requested number of bands might be too large (exceeds half the
            number of frequencies).")
  }
  if (s < 1) {
    s <- 0
    warning("Shift window for population initialization was set to zero due
            a small provided value; all solutions will initialize to the
            evenly spaced grid.")
  }
  if (s > nfreq) {
    stop("Shift window for population initialization is set too high.")
  }


  # equally spaced grid to base initial endpoints around
  eqspaced <- floor(seq(1, nfreq + 1, length = nbands + 1))
  if (any(diff(eqspaced) == 0)) {
    stop(paste0("Population initialization failed. ",
         "Parameter nbands probably set too high."))
  }

  # initialize each member of the population
  pop <- list()
  for (i in seq_len(popsize)) {
    ch <- matrix(nrow = nclust, ncol = 2 * nbands - 1)

    # initialize endpoints; same for each cluster
    endpoints <- rep(0, nendpoints)
    for (b in 2:nbands) {
      left <- max(eqspaced[b - 1] + 1, eqspaced[b] - s)
      right <- min(eqspaced[b + 1] - 1, eqspaced[b] + s)
      endpoints[b - 1] <- sample(seq(left, right), 1)
    }
    endpoints <- sort(endpoints) # ensure they're sorted

    # initialize collapsed measures of power by randomly selecting a
    # replicate spectrum and collapsing it using the endpoints set above
    collapsed <- matrix(nrow = nclust, ncol = nbands)
    for (j in seq_len(nclust)) {
      rspec <- spec[, sample(ncol(spec), 1)]
      ch[j, 1:nendpoints] <- endpoints
      ch[j, (nendpoints + 1):ncol(ch)] <- avg_collapsed_by_index(rspec, endpoints)
    }

    pop[[i]] <- ch # add to population
  }
  return(pop)
}
