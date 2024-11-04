pop_init <- function(popsize, nsubpop, nbands, spec) {
  nfreq <- dim(spec)[1]; nrep <- dim(spec)[2]
  nendpoints <- nbands - 1  # number of frequency band boundaries

  # checks on input parameters
  if (popsize < 1) stop("Population size must be at least 1.")
  if (nsubpop < 1) stop("Number of subpopulations must be at least 1.")
  if (nsubpop > nrep) stop("Number of subpopulations cannot exceed the number of
                          subjects.")
  if (nsubpop > nrep / 2 & nsubpop > 1) warning("Number of subpopulations exceeds half the
                                 number of replicates.")
  if (nbands < 2) stop("Number of frequency bands must be at least 2.")
  if (nbands > nfreq) stop("Number of frequency bands cannot exceed the number
                           of frequencies.")
  if (nbands > nfreq / 2) warning("Number of frequency bands exceeds half the
                                  number of frequencies.")

  # equally spaced grid for the initial boundaries
  # should be unique as long as nbands < nfreq (checked above)
  # eqspaced <- floor(seq(1, nfreq + 1, length = nbands + 1))
  # eqspaced <- eqspaced[-c(1, length(eqspaced))] # keep only interior boundaries

  # initialize each member of the population
  pop <- matrix(nrow = popsize, ncol = nsubpop * (2 * nbands - 1))
  for (i in seq_len(popsize)) {
    # initialize current chromosome
    p <- matrix(nrow = nsubpop, ncol = 2 * nbands - 1)
    # initialize collapsed measures of power by randomly selecting a
    # replicate spectrum and collapsing it using the endpoints set above
    collapsed <- matrix(nrow = nsubpop, ncol = nbands)
    for (j in seq_len(nsubpop)) {
      rand_spec <- spec[, sample(nrep, 1), drop = FALSE]
      p[j, 1:nendpoints] <- sort(sample(2:nfreq, size = nbands - 1))
      p[j, (nendpoints + 1):ncol(p)] <- avg_summary(rand_spec, p[j, 1:nendpoints])
    }
    pop[i,] <- as.vector(t(p)) # add to population
  }
  return(pop)
}
