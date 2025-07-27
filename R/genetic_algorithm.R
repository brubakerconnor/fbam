#' Genetic algorithm
#'
#' Optimization of FBAM loss function
#'
#' @param spec Replicate-specific spectrum estimates (columns). A single vector is treated as one replicate.
#' @param nbands Number of frequency bands.
#' @param nsubpop Number of subpopulations.
#' @param popsize Number of chromosomes per population.
#' @param factor Multiplier of the population size that sets how many offspring to generate each generation for subsequent selection.
#' @param pmutate_init Initial mutation rate.
#' @param rate Rate at which to decrease mutation rate in terms of number of generations it takes to halve the rate.
#' @param nislands Number of independent islands.
#' @param nmigrants Number of migrants to move between islands.
#' @param epoch How often (in generations) to perform migration.
#' @param maxgen Terminate evolution after this number of generations.
#' @param maxrun Terminate evolution after non-improvement in fitness after this many generations.
#' @param tol Sets the threshold for 'improvement'.
#' @param verbose Print information to console?
#'
#' @returns A object of class "ga" with the following elements:
#' \describe{
#'  \item{nfreq}{Number of frequencies (number of rows in 'spec').}
#'  \item{labels}{Subpopulation assignment vector.}
#'  \item{endpoints}{Frequency band endpoints for each subpopulation as integer-valued indices of rows in spec (row-wise).}
#'  \item{avg_summary}{Average summary measures for each subpopulation (rows). Columns correspond to frequency band.}
#'  \item{rep_summary}{Replicate-speciifc summary measures (rows). Columns correspond to frequency bands.}
#'  \item{loss}{Loss function value at solution.}
#'  \item{validation}{Validation criterion evaluated on solution.}
#'  \item{avg_loss}{Vector of average loss function values for each generation.}
#'  \item{min_loss}{Vector of minimum loss function values for each generation.}
#'  \item{ninfeasible}{Number of infeasible solutions in each generation.}
#'  \item{params}{List of algorithm parameters used.}
#' }
#' @export
#'
#' @examples
genetic_algorithm <- function(
    spec, nbands = 2, nsubpop = 1, popsize = 50, factor = 2,
    pmutate_init = 0.5, rate = 50, nislands = 6, nmigrants = 5,
    epoch = 50, maxgen = 500, maxrun = 500, tol = 1e-3, verbose = FALSE
) {


  # input checks ---------------------------------------------------------------

  if (is.vector(spec)) spec <- as.matrix(spec)
  if (nbands < 2) {
    stop("'nbands' must be at least 2.")
  }
  if (nsubpop < 1 | nsubpop > ncol(spec)) {
    stop("'nsubpop' must be at least 1 and no greater than the number of
         provided replicates (columns in 'spec').")
  }
  if (popsize < 1) {
    stop("'popsize' must be at least 1. Using recommended setting of 50.")
    popsize <- 50
  }
  if (maxgen < 1) {
    warning("'maxgen' must be at least 1. Using recommended setting of 500.")
    maxgen <- 500
  }
  if (maxrun < 1) {
    stop("'maxrun' must be at least 1. Setting to value of 'maxgen'.")
    maxrun <- maxgen
  }
  if (maxrun > maxgen) {
    warning("maxrun greater than maxgen so convergence will be declared only
            when maxgen number of generations have been completed.")
    maxrun <- maxgen
  }
  if (tol > 1) {
    tol <- 1e-3
    warning("'tol' too high. Setting to 1e-3.")
  }
  if (tol < 1e-3) {
    warning("'tol' was set below recommended value of 1e-3 which may
            lead to longer convergence times.")
  }
  nrep <- ncol(spec)

  # initialization -------------------------------------------------------------

  if(verbose) cat('[INFO] Initialization population(s)\n')
  pops <- lapply(seq(nislands), function(x) {
    pop <- pop_init(popsize, nsubpop, nbands, spec)
    fitvals <- evaluate(pop, spec, nsubpop, nbands)
    elite <- pop[which.max(fitvals), ]
    elite_fitval <- fitvals[which.max(fitvals)]
    return(list(pop = pop, fitvals = fitvals, elite = elite,
                elite_fitval = elite_fitval))
  })

  gen <- run <- 0
  avgfit <- maxfit <- ninfeasible <- rep(0, maxgen + 1)
  avgfit[gen + 1] <- mean(unlist(lapply(pops, function(x) mean(x$fitvals))))
  maxfit[gen + 1] <- max(unlist(lapply(pops, function(x) max(x$fitvals))))
  ninfeasible[gen + 1] <- sum(unlist(lapply(pops, function(x) sum(x$fitvals == 1e-50))))
  if (verbose) {
    cat(paste0("[INFO] Generation: ", gen,
               "\tAverage Loss: ", round(1 / avgfit[gen + 1], 6),
               "\tMinimum Loss: ", round(1 / maxfit[gen + 1], 6)), "\n\n")
  }

  pmutate <- pmutate_init
  alpha <- log(2) / rate


  # evolution ------------------------------------------------------------------

  if(verbose) cat('[INFO] Starting evolution\n')
  while (gen < maxgen & run <= maxrun) {
    gen <- gen + 1 # increase generation count

    pops <- lapply(pops, function(x) {
      # mutation step
      new_pop <- matrix(nrow = factor * popsize, ncol = ncol(x$pop))
      for (i in 1:nrow(x$pop)) {
        for (r in seq(factor)) {
          p <- matrix(x$pop[i,], nrow = nsubpop, byrow = TRUE)
          tmp <- mutate(p, pmutate, spec)
          tmp_labels <- l2assign(spec, matrix(tmp, nrow = nsubpop, byrow = TRUE))
          if (length(unique(tmp_labels)) != nsubpop) {
            new_pop[factor * (i - 1) + r,] <- x$pop[i,]
          } else {
            new_pop[factor * (i - 1) + r,] <- tmp
          }
        }
      }
      fitvals <- evaluate(new_pop, spec, nsubpop, nbands)
      ranks <- rank(-fitvals, ties.method = "random")
      pop <- new_pop[ranks <= popsize,]
      fitvals <- fitvals[ranks <= popsize]

      # elitism step
      if (x$elite_fitval > max(fitvals)) {
        # the previous elite is better than the current elite
        # replace current least fit member by the previous elite
        pop[which.min(fitvals), ] <- x$elite
        fitvals[which.min(fitvals)] <- x$elite_fitval
      } else {
        # the current elite is better than the previous elite
        # replace previous elite with current elite
        x$elite <- pop[which.max(fitvals),]
        x$elite_fitval <- fitvals[which.max(fitvals)]
      }
      return(list(pop = pop, fitvals = fitvals, elite = x$elite,
                  elite_fitval = x$elite_fitval))
    })

    if (gen %% epoch == 0) {
      if (verbose) cat('[INFO] Performing migration\n')
      # perform migration
      ind <- sample(1:popsize, size = nmigrants, replace = F)
      samp_prev <- pops[[1]]$pop[ind,]
      for (i in 2:nislands) {
        samp <- pops[[i]]$pop[ind,] <- samp_prev
        samp_prev <- samp
      }
      pops[[1]]$pop[ind,] <- samp_prev
    }

    # update mutation rate
    pmutate <- pmutate_init * exp(-alpha * gen)

    # reporting
    avgfit[gen + 1] <- mean(unlist(lapply(pops, function(x) mean(x$fitvals))))
    maxfit[gen + 1] <- max(unlist(lapply(pops, function(x) max(x$fitvals))))
    ninfeasible[gen + 1] <- sum(unlist(lapply(pops, function(x) sum(x$fitvals == 1e-50))))
    if ((maxfit[gen + 1] - maxfit[gen]) / maxfit[gen] < tol) {
      run <- run + 1
    } else {
      run <- 0
    }
    if (verbose & gen %% 25 == 0) {
      cat(paste0("[INFO] Generation: ", gen,
                 "\tAverage Loss: ", round(1 / avgfit[gen + 1], 6),
                 "\tMinimum Loss: ", round(1 / maxfit[gen + 1], 6)), "\n")
    }
  }
  if(verbose) cat('[INFO] Evolution terminated\n')

  # termination of evolution, extract solution ---------------------------------

  ## extract raw solution ----
  pop <- do.call(rbind, lapply(pops, function(x) x$pop))
  fitvals <- unlist(lapply(pops, function(x) x$fitvals))
  solution <- matrix(pop[which.max(fitvals),], nrow = nsubpop, byrow = T)
  solution_fitness <- fitvals[which.max(fitvals)]

  ## assign replicates to subpopulations ----
  if (nsubpop > 1) {
    labels <- l2assign(spec, solution)
  } else {
    labels <- rep(1, nrep)
  }

  ## extract frequency bands and summary measures ----
  endpoints <- solution[, 1:(nbands - 1), drop = FALSE]
  avg_summary_mat <- matrix(nrow = nsubpop, ncol = nbands)
  rep_summary_mat <- matrix(nrow = dim(spec)[2], ncol = nbands)
  for (j in 1:nsubpop) {
    subpop <- spec[, labels == j, drop = FALSE]
    avg_summary_mat[j,] <- avg_summary(subpop, endpoints[j,])
    rep_summary_mat[labels == j,] <- rep_summary(subpop, endpoints[j,])
  }

  ## compute validation index ----
  validation <- validation_index(spec, labels, endpoints, avg_summary_mat)

  ## output ----
  out <- list(
    nfreq = nrow(spec),
    labels = labels,
    endpoints = endpoints,
    avg_summary = avg_summary_mat,
    rep_summary = rep_summary_mat,
    loss = 1 / fitvals[which.max(fitvals)],
    validation = validation,
    avg_loss = 1 / avgfit[1:(gen + 1)],
    min_loss = 1 / maxfit[1:(gen + 1)],
    ninfeasible = ninfeasible[1:(gen + 1)],
    params = list(
      nsubpop = nsubpop,
      nbands = nbands,
      popsize = popsize,
      factor = factor,
      pmutate_init = pmutate_init,
      pmutate_final = pmutate,
      maxgen = maxgen,
      maxrun = maxrun,
      tol = tol))
  class(out) <- "ga"
  return(out)
}

#' Print method for object of class "ga"
#'
#' @param out An object of class "ga"
#' @param freq Vector of observed frequencies (useful when sampling rate is different from 1)
#' @method print ga
#' @export
print.ga <- function(out, freq = NULL) {
  if (is.null(freq)) {
    # assumed that sampling rate was 1
    freq <- seq(0, 0.5, length = out$nfreq + 2)
    freq <- freq[freq != 0 & freq != 0.5]
  }
  cat('Total number of replicates:', length(out$labels), '\n')
  cat('Requested number of subpopulations:', out$params$nsubpop, '\n')
  cat('Requested number of frequency bands:', out$params$nbands, '\n')
  cat('Loss function at solution:', round(out$loss, 4), '\n\n')
  for (j in 1:out$params$nsubpop) {
    bands <- round(freq[out$endpoints[j,]], 4)
    summaries <- round(out$avg_summary[j,], 2)
    nrep_subpop <- sum(out$labels == j)
    cat('Subpopulation', j, '\n')
    cat('Number of assigned replicates: ', nrep_subpop, '\n')
    cat('Frequency bands: ', paste0(bands, ' Hz'), '\n')
    cat('Average summary measures of power: ', summaries, '\n\n')
  }
}

#' Print method for object of class "ga"
#'
#' Plot evolution history (average and minimum loss functions in each generation) of the algorithm
#'
#' @param out An object of class "ga"
#' @method plot ga
#' @export
plot.ga <- function(out) {
  ngen <- length(out$avg_loss)
  gen <- seq(ngen)
  plot(gen, out$avg_loss, type = 'o',
       ylim = c(min(out$min_loss), max(out$avg_loss)),
       xlab = "Generation", ylab = "Loss", main = "Evolution History")
  lines(gen, out$min_loss, type = 'o', col = 'red')
  grid()
  legend('topright', legend = c('Average loss', 'Minimum loss'),
         col = 1:2, lty = 1)
}
