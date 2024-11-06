ga <- function(X, nbands, nsubpop, popsize, pmutate, maxgen, maxrun, tol,
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
    cat(paste0("START \t\t AVG FIT: ", round(avgfit[gen + 1], 6),
      "\t\t MAX FIT: ", round(maxfit[gen + 1], 6)), "\n")
  }

  # iterate until convergence
  while (gen < maxgen & run <= maxrun) {
    gen <- gen + 1 # increase generation count
    # pop <- select_parents(pop, 2 * popsize, fitvals)
    R <- 4
    new_pop <- matrix(nrow = R * popsize, ncol = ncol(pop))
    for (i in 1:nrow(pop)) {
      for (r in 1:R) {
        p <- matrix(pop[i,], nrow = nsubpop, byrow = TRUE)
        new_pop[R * (i - 1) + r,] <- mutate(p, pmutate, spec)
      }
      # p1 <- matrix(pop[i,], nrow = nsubpop, byrow = TRUE)
      # p2 <- matrix(pop[i,], nrow = nsubpop, byrow = TRUE)
      # new_pop[2 * (i - 1) + 1,] <- mutate(p1, pmutate, spec)
      # new_pop[2 * (i - 1) + 2,] <- mutate(p2, pmutate, spec)
    }
    # evaluation and rank selection of offspring
    # evaluate the children that were created via crossover and/or mutation
    # of the selected parents.
    # children are evaluated and ranked. only the top popsize of them are
    # allowed to survive into the next generation.
    # the previous generation from which parents were selected is discarded.
    fitvals <- evaluate(new_pop, spec, nsubpop, nbands)
    ranks <- rank(-fitvals, ties.method = "random")
    pop <- new_pop[ranks <= popsize,]
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
      cat(paste0("GEN ", gen, "\t\t AVG FIT: ", round(avgfit[gen + 1], 6),
        "\t\t MAX FIT: ", round(maxfit[gen + 1], 6)), "\n")
    }
  }
  solution <- matrix(pop[which.max(fitvals),], nrow = nsubpop, byrow = T)
  solution_fitness <- fitvals[which.max(fitvals)]
  if (nsubpop > 1) {
    labels <- l2assign(spec, solution)
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
  ga <- list(freq = mtfreq, spec = spec, labels = labels, endpoints = endpoints,
    avg_summary = avg_summary_mat, rep_summary = rep_summary_mat,
    objective = 1 / fitvals[which.max(fitvals)], validation = validation,
    avgfit = avgfit[1:(gen + 1)], maxfit = maxfit[1:(gen + 1)],
    ninfeasible = ninfeasible[1:(gen + 1)],
    params = list(nsubpop = nsubpop, nbands = nbands, popsize = popsize,
                  pmutate = pmutate, maxgen = maxgen, maxrun = maxrun,
                  tol = tol))
  class(ga) <- "ga_output"
  return(ga)
}

#' @export
print.ga_output <- function(ga) {
  cat('Number of replicates provided:', ncol(ga$spec), '\n')
  cat('nsubpop =', ga$params$nsubpop, ', nbands =', ga$params$nbands)
  for (j in 1:ga$params$nsubpop) {
    cat('\n\nSubpopulation j =', j, '\n\tBoundaries:',
        round(ga$endpoints[j,], 2),
        '\n\tAverage Summary Measures:', round(ga$avg_summary[j,], 2))
  }
}

#' @export
plot.ga_output <- function(ga, type = 'solution') {
  if (type == 'solution') {
    nsubpop <- ga$params$nsubpop; nbands <- ga$params$nbands
    if (ga$params$nsubpop > 1) oldpar <- par(mfrow = c(ceiling(nsubpop / 2), 2))
    for (j in 1:nsubpop) {
      subpop_spec <- ga$spec[, ga$labels == j, drop = FALSE]
      matplot(ga$freq, subpop_spec,
              type = 'l', lty = 1, lwd = 0.35, col = 'grey80',
              ylim = c(0, max(ga$spec)),
              xlab = 'Frequency', ylab = 'Power',
              main = ifelse(nsubpop > 1, paste0('Subpopulation ', j), 'All Replicates'))
      abline(v = ga$endpoints[j,], lwd = 1.5, lty = 2)
      endpoints <- c(0, ga$endpoints[j,], 0.5)
      for (l in 1:nbands) {
        band_freq <- ga$freq[ga$freq >= endpoints[l] & ga$freq < endpoints[l + 1]]
        lines(band_freq, rep(ga$avg_summary[j, l], length(band_freq)),
              col = 'red', lwd = 1.5)
      }
    }
    if (ga$params$nsubpop > 1) par(oldpar)
  } else if (type == 'history') {
    num_gen <- length(ga$avgfit)
    plot(1:num_gen, 1 / ga$avgfit, type = 'o', col = 'grey50',
         ylim = c(min(c(1/ga$avgfit, 1/ga$maxfit)), max(c(1/ga$avgfit, 1/ga$maxfit))),
         xlab = 'Generation Number', ylab = 'Objective Function',
         main = 'Evolution History')
    lines(1:num_gen, 1 / ga$maxfit, type = 'o', col = 'red')
    grid()
    legend('topright', legend = c('average objective', 'minimum objective'),
           lty = 1, col = c(1, 2))
  } else {
    stop("Unknown plot type. Please specify 'solution' or 'history'.")
  }

}
