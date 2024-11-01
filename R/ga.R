ga <- function(X, nsubpop, nbands, popsize, pmutate, maxgen, maxrun, tol,
               ntapers, verbose) {

  # estimate spectra from matrix of time series data X
  if (is.vector(X)) X <- as.matrix(X) # treat vector as one replicate
  mtout <- sine_mt(X, ntapers)
  mtfreq <- mtout$mtfreq; spec <- mtout$mtspec
  nfreq <- length(mtfreq); nrep <- dim(spec)[2]

  # checks on maxgen, maxrun, and tolerance
  if (maxgen < 1) {
    warning("maxgen must be at least 1. Using recommended setting of 500.")
    maxgen <- 500
  }
  if (maxrun < 1) {
    stop("maxrun must be at least 1. Using recommended setting of 100")
    maxrun <- 100
  }
  if (maxrun > maxgen) {
    warning("maxrun greater than maxgen so convergence will be declared only
            when maxgen number of generations have been completed.")
    maxrun <- maxgen
  }
  if (tol > 1) {
    tol <- 1
    warning("Tolerance set to maximum of 1. Recommended setting is 1e-3.")
  }
  if (tol < 1e-3) {
    warning("Tolerance was set below recommended value of 1e-3 which may
            lead to longer convergence times.")
  }

  # step 1: population initialization
  if (popsize < 1) {
    stop("popsize must be at least 1. Using recommended setting of 50.")
    popsize <- 50
  }
  pop <- pop_init(popsize, nsubpop, nbands, spec)

  # step 2: initial evaluation of population and initialization of elite
  fitvals <- evaluate(pop, spec, nsubpop, nbands)
  elite <- pop[which.max(fitvals), ]
  elite_fitval <- fitvals[which.max(fitvals)]

  # initialize monitoring
  gen <- run <- 0
  avgfit <- maxfit <- ninfeasible <- rep(0, maxgen + 1)
  avgfit[gen + 1] <- mean(fitvals)
  maxfit[gen + 1] <- max(fitvals)
  ninfeasible[gen + 1] <- sum(fitvals == 1e-50)
  if (verbose) {
    cat(paste0("START \t\t AVG FIT: ", round(avgfit[gen + 1], 5),
      "\t\t MAX FIT: ", round(maxfit[gen + 1], 5)), "\n")
  }

  # iterate until convergence
  while (gen < maxgen & run <= maxrun) {
    gen <- gen + 1 # increase generation count
    pop <- select_parents(pop, 2 * popsize, fitvals)

    # mutation
    for (i in 1:nrow(pop)) {
      # mutate the current member
      p <- matrix(pop[i,], nrow = nsubpop, byrow = TRUE)
      pop[i,] <- mutate(p, pmutate, spec)
      # cluster the spectra on the basis of this mutated solution
      # if the mutated solution results in an infeasible solution,
      # the mutations are discarded and the current solution is left as is
      # labels_ <- apply(dist_mat(spec, out), 1, which.min)
      # if (length(unique(mutated_labels)) == nclust) parents[[i]] <- mutated
    }

    # evaluation and rank selection of offspring
    # evaluate the children that were created via crossover and/or mutation
    # of the selected parents.
    # children are evaluated and ranked. only the top popsize of them are
    # allowed to survive into the next generation.
    # the previous generation from which parents were selected is discarded.
    fitvals <- evaluate(pop, spec, nsubpop, nbands)
    ranks <- rank(-fitvals, ties.method = "random")
    pop <- pop[ranks <= popsize,]
    fitvals <- fitvals[ranks <= popsize]

    # elitism
    if (elite_fitval > max(fitvals)) {
      # the previous elite is better than the current elite
      # replace current least fit member by the previous elite
      pop[which.min(fitvals), ] <- elite
      fitvals[which.min(fitvals)] <- elite_fitval
    } else {
      # the current elite is better than the previous elite
      # replace previous elite with current elite
      elite <- pop[which.max(fitvals),]
      elite_fitval <- fitvals[which.max(fitvals)]
    }

    # monitoring
    avgfit[gen + 1] <- mean(fitvals)
    maxfit[gen + 1] <- max(fitvals)
    ninfeasible[gen + 1] <- sum(fitvals == 1e-50)
    if ((maxfit[gen + 1] - maxfit[gen]) / maxfit[gen] < tol) {
      run <- run + 1
    } else {
      run <- 0
    }
    if (verbose & gen %% 25 == 0) {
      cat(paste0("GEN ", gen, "\t\t AVG FIT: ", round(avgfit[gen + 1], 5),
        "\t\t MAX FIT: ", round(maxfit[gen + 1], 5)), "\n")
    }
  }
  solution <- matrix(pop[which.max(fitvals),], nrow = nsubpop, byrow = T)
  solution_fitness <- fitvals[which.max(fitvals)]
  if (nsubpop > 1) {
    labels <- l2_assign(spec, solution)
  } else {
    labels <- rep(1, nrep)
  }
  endpoints_index <- solution[, 1:(nbands - 1), drop = FALSE]
  endpoints <- endpoints_index
  avg_summary_mat <- matrix(nrow = nsubpop, ncol = nbands)
  rep_summary_mat <- matrix(nrow = dim(spec)[2], ncol = nbands)
  for (j in 1:nsubpop) {
    subpop <- spec[, labels == j, drop = FALSE]
    endpoints[j,] <- mtfreq[endpoints_index[j,]]
    avg_summary_mat[j,] <- avg_summary(subpop, endpoints_index[j,])
    rep_summary_mat[labels == j,] <- rep_summary(subpop, endpoints_index[j,])
  }
  validation <- validation_index(spec, labels, endpoints_index, avg_summary_mat)
  return(list(
    freq = mtfreq,
    spec = spec,
    labels = labels,
    endpoints = endpoints,
    avg_summary = avg_summary_mat,
    rep_summary = rep_summary_mat,
    objective = 1 / fitvals[which.max(fitvals)],
    validation = validation,
    avgfit = avgfit[1:(gen + 1)],
    maxfit = maxfit[1:(gen + 1)],
    ninfeasible = ninfeasible[1:(gen + 1)],
    params = list(
      nsubpop = nsubpop,
      nbands = nbands,
      popsize = popsize,
      pmutate = pmutate,
      maxgen = maxgen,
      maxrun = maxrun,
      tol = tol
    )
  ))
}
