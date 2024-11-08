ea <- function(X, nbands, nsubpop, popsize, maxgen, maxrun, tol,
               ntapers, sample_rate, verbose) {

  n_islands <- 6
  n_migrants <- 5
  epoch <- 50 # migration frequency

  # estimate spectra from matrix of time series data X
  if (is.vector(X)) X <- as.matrix(X) # treat vector as one replicate
  mtout <- sine_mt(X, ntapers, sample_rate)
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

  # step 1: population initialization for each of the islands
  # initialize fitness values and elites
  if (popsize < 1) {
    stop("popsize must be at least 1. Using recommended setting of 50.")
    popsize <- 50
  }
  pops <- lapply(1:n_islands, function(x) {
    pop <- pop_init(popsize, nsubpop, nbands, spec)
    fitvals <- evaluate(pop, spec, nsubpop, nbands)
    elite <- pop[which.max(fitvals), ]
    elite_fitval <- fitvals[which.max(fitvals)]
    return(list(pop = pop, fitvals = fitvals,
                elite = elite, elite_fitval = elite_fitval))
  })

  # initialize monitoring (aggregated across islands)
  gen <- run <- 0
  avgfit <- maxfit <- ninfeasible <- rep(0, maxgen + 1)
  avgfit[gen + 1] <- mean(unlist(lapply(pops, function(x) mean(x$fitvals))))
  maxfit[gen + 1] <- max(unlist(lapply(pops, function(x) max(x$fitvals))))
  ninfeasible[gen + 1] <- sum(unlist(lapply(pops, function(x) sum(x$fitvals == 1e-50))))
  if (verbose) {
    cat(paste0("GEN ", gen, "\t\t AVG FIT: ", round(avgfit[gen + 1], 6),
               "\t\t MAX FIT: ", round(maxfit[gen + 1], 6)), "\n")
  }

  # adaptive mutation rate
  pmutate_init <- pmutate <- 0.5
  alpha <- log(2) / 50  # halve mutation rate every 50 generations

  # iterate until convergence
  while (gen < maxgen & run <= maxrun) {
    gen <- gen + 1 # increase generation count

    pops <- lapply(pops, function(x) {
      R <- 2 # popsize
      new_pop <- matrix(nrow = R * popsize, ncol = ncol(x$pop))
      for (i in 1:nrow(x$pop)) {
        for (r in 1:R) {
          p <- matrix(x$pop[i,], nrow = nsubpop, byrow = TRUE)
          tmp <- mutate(p, pmutate, spec)
          tmp_labels <- l2assign(spec, matrix(tmp, nrow = nsubpop, byrow = TRUE))
          if (length(unique(tmp_labels)) != nsubpop) {
            new_pop[R * (i - 1) + r,] <- x$pop[i,]
          } else {
            new_pop[R * (i - 1) + r,] <- tmp
          }
        }
      }
      fitvals <- evaluate(new_pop, spec, nsubpop, nbands)
      ranks <- rank(-fitvals, ties.method = "random")
      pop <- new_pop[ranks <= popsize,]
      fitvals <- fitvals[ranks <= popsize]

      # elitism
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
      return(list(pop = pop, fitvals = fitvals,
                  elite = x$elite, elite_fitval = x$elite_fitval))
    })

    if (gen %% epoch == 0) {
      # perform migration
      ind <- sample(1:popsize, size = n_migrants, replace = F)
      samp_prev <- pops[[1]]$pop[ind,]
      for (i in 2:n_islands) {
        samp <- pops[[i]]$pop[ind,] <- samp_prev
        samp_prev <- samp
      }
      pops[[1]]$pop[ind,] <- samp_prev
    }

    pmutate <- pmutate_init * exp(-alpha * gen)

    # monitoring
    avgfit[gen + 1] <- mean(unlist(lapply(pops, function(x) mean(x$fitvals))))
    maxfit[gen + 1] <- max(unlist(lapply(pops, function(x) max(x$fitvals))))
    ninfeasible[gen + 1] <- sum(unlist(lapply(pops, function(x) sum(x$fitvals == 1e-50))))
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
  pop <- do.call(rbind, lapply(pops, function(x) x$pop))
  fitvals <- unlist(lapply(pops, function(x) x$fitvals))
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
  ea_output <- list(freq = mtfreq, spec = spec, labels = labels, endpoints = endpoints,
    avg_summary = avg_summary_mat, rep_summary = rep_summary_mat,
    objective = 1 / fitvals[which.max(fitvals)], validation = validation,
    avgfit = avgfit[1:(gen + 1)], maxfit = maxfit[1:(gen + 1)],
    ninfeasible = ninfeasible[1:(gen + 1)],
    params = list(nsubpop = nsubpop, nbands = nbands, popsize = popsize,
                  pmutate = pmutate, maxgen = maxgen, maxrun = maxrun,
                  tol = tol, sample_rate = sample_rate))
  class(ea_output) <- "ea_output"
  return(ea_output)
}

#' @export
print.ea_output <- function(ea) {
  cat('Number of replicates provided:', ncol(ea$spec), '\n')
  cat('nsubpop = ', ea$params$nsubpop, ', nbands = ', ea$params$nbands, sep = "")
  for (j in 1:ea$params$nsubpop) {
    nrep_subpop <- ncol(ea$spec[, ea$labels == j, drop = FALSE])
    cat('\n\nSubpopulation j =', j, '\nNumber of replicates:', nrep_subpop,
        '\n\tBoundaries:',
        paste0(round(ea$endpoints[j,], 4), ' Hz'),
        '\n\tAverage Summary Measures:', round(ea$avg_summary[j,], 2))
  }
}

#' @export
plot.ea_output <- function(ea_output, type = 'solution') {
  if (type == 'solution') {
    `%>%` <- dplyr::`%>%`
    nrep <- length(ea_output$labels)
    nsubpop <- ea_output$params$nsubpop; nbands <- ea_output$params$nbands
    # data frame for multitaper spectral estimates
    mt <- cbind(as.data.frame(t(ea_output$spec)), ifelse(
      rep(nsubpop > 1, nrep),
      paste0("Subpopulation ", ea_output$labels),
      rep("All Replicates", length(ea_output$labels))
    ))
    mt <- mt %>% dplyr::mutate(row_id = dplyr::row_number())
    names(mt)[1:(ncol(mt) - 2)] <- ea_output$freq
    names(mt)[ncol(mt) - 1] <- "label"
    mt_long <- mt %>%
      tidyr::pivot_longer(cols = -c(label, row_id),
                          names_to = "freq", values_to = "spec") %>%
      dplyr::mutate(freq = as.numeric(freq)) # for plotting

    # data frame for internal band boundaries
    boundaries <- data.frame(
      label = ifelse(
        rep(nsubpop > 1, nsubpop * (nbands - 1)),
        paste0("Subpopulation ", rep(1:nsubpop, each = nbands - 1)),
        rep("All Replicates", nsubpop * (nbands - 1))),
      boundary = as.vector(t(ea_output$endpoints))
    )

    # data frame for collapsed measures with start and end points
    endpoints <- matrix(min(ea_output$freq), nrow = nsubpop, ncol = nbands + 1)
    endpoints[1:nsubpop, 2:nbands] <- ea_output$endpoints
    endpoints[1:nsubpop, nbands + 1] <- max(ea_output$freq)
    collapsed <- data.frame(
      label = ifelse(
        rep(nsubpop > 1, nsubpop * nbands),
        paste0("Subpopulation ", rep(1:nsubpop, each = nbands)),
        rep("All Replicates", nsubpop * nbands)
      ),
      xstart = as.vector(t(endpoints[, 1:nbands])),
      xend = as.vector(t(endpoints[, 2:(nbands + 1)])),
      yvalue = as.vector(t(ea_output$avg_summary))
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = freq, y = spec, group = row_id), mt_long,
                linewidth = 0.25) +
      ggplot2::scale_x_continuous(
        breaks = seq(from = 0, to = 0.5 * ea_output$params$sample_rate, by = 0.25)
      ) +
      ggplot2::facet_wrap(~label) +
      ggplot2::geom_segment(ggplot2::aes(x = xstart, xend = xend, y = yvalue, yend = yvalue),
                   collapsed, color = 'red', size = 1) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = boundary), boundaries, linetype = 'dashed',
                 size = 1) +
      ggplot2::labs(x = "Hz", y = "Power") +
      ggplot2::theme_bw()

    return(p)


  } else if (type == 'history') {
    num_gen <- length(ea_output$avgfit)
    plot(1:num_gen, 1 / ea_output$avgfit, type = 'o', col = 'grey50',
         ylim = c(min(c(1/ea_output$avgfit, 1/ea_output$maxfit)), max(c(1/ea_output$avgfit, 1/ea_output$maxfit))),
         xlab = 'Generation Number', ylab = 'Objective Function',
         main = 'Evolution History')
    lines(1:num_gen, 1 / ea_output$maxfit, type = 'o', col = 'red')
    grid()
    legend('topright', legend = c('average objective', 'minimum objective'),
           lty = 1, col = c(1, 2))
  } else {
    stop("Unknown plot type. Please specify 'solution' or 'history'.")
  }

}
