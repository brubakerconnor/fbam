library(future.apply)

fbam <- function(x, nbands = 2, nsubpop = 1, ntapers = NULL, ncores = 1, ...) {
  if (is.vector(x)) {
    x <- matrix(x)
  }
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
  sine_mt_out <- sine_mt(x, ntapers)
  freq <- sine_mt_out$mtfreq
  spec <- sine_mt_out$mtspec
  
  ## set up parallelism
  plan(multisession, workers = ncores)

  # if nsubpop = 1 is one of the requested values:
  # cannot compare solutions with 1 subpop and those with more than 1
  # grid search over varying nbands with nsubpop = 1 fixed done first
  if (1 %in% nsubpop) {
    nsubpop1_grid <- future_lapply(
      nbands,
      \(x) ga(spec, nbands = x, nsubpop = 1, ...),
      future.seed = TRUE
    )
    nsubpop1_s <- unlist(lapply(nsubpop1_grid, function(x) x$s_band))
    nsubpop1_selected <- nsubpop1_grid[[which.min(nsubpop1_s)]]
  } else {
    nsubpop1_grid <- NULL
    nsubpop1_selected <- NULL
  }

  # if nsubpop contains values greater than 1:
  # can compare solutions corresponding to any value of nsubpop > 1
  # grid search over varying nbands AND nsubpop as requested by user
  if (length(nsubpop[nsubpop > 1]) > 0) {
    param_grid <- expand.grid(nbands = nbands, nsubpop = nsubpop[nsubpop > 1])
    grid <- future_lapply(
      1:nrow(param_grid),
      function(i) {
        nsubpop <- param_grid$nsubpop[i]
        nbands <- param_grid$nbands[i]
        ga(spec, nbands = nbands, nsubpop = nsubpop, ...)
      },
      future.seed = TRUE
    )

    if (length(nsubpop) == 1) {
      ## only one nsubpop value - use s_band for selection
      s <- unlist(lapply(grid, function(x) x$s_band))
      solution <- grid[[which.min(s)]]
    } else if (length(nbands) == 1) {
      ## only one nbands value - use s_pop for selection
      s <- unlist(lapply(grid, function(x) x$s_pop))
      solution <- grid[[which.min(s)]]
    } else {
      ## use s_joint after scaling
      s_band <- unlist(lapply(grid, function(x) x$s_band))
      s_band_scaled <- s_band / max(s_band)

      s_pop <- unlist(lapply(grid, function(x) x$s_pop))
      s_pop_scaled <- s_pop / max(s_pop)

      s_joint <- s_band + s_pop
      solution <- grid[[which.min(s_joint)]]
    }
  } else {
    grid <- NULL
    solution <- NULL
  }
  return(list(
    grid = grid,
    solution = solution,
    nsubpop1_grid = nsubpop1_grid,
    nsubpop1_solution = nsubpop1_selected
  ))
}

ga <- function(
  spec,
  nbands = 2,
  nsubpop = 1,
  popsize = 50,
  pmutate = 0.15,
  nislands = 10,
  nmigrants = 5,
  epoch = 50,
  maxgen = 500,
  maxrun = 100,
  tol = 0.01,
  verbose = TRUE
) {
  if (is.vector(spec)) {
    spec <- as.matrix(spec)
  }
  if (nbands < 2) {
    stop("'nbands' must be at least 2.")
  }
  if (nbands > nrow(spec)) {
    stop("'nbands' cannot exceed number of frequencies (rows in 'spec').")
  }
  if (nsubpop < 1 | nsubpop > ncol(spec)) {
    stop("'nsubpop' must be between 1 and number of columns in 'spec'.")
  }
  if (popsize < 1) {
    stop("'popsize' must be positive.")
  }
  if (maxgen < 1) {
    stop("'maxgen' must be positive.")
  }
  if (maxrun < 1) {
    stop("'maxrun' must be positive.")
  }
  if (maxrun > maxgen) {
    warning(
      "maxrun greater than maxgen so convergence will be declared only
            when maxgen number of generations have been completed."
    )
  }
  if (tol > 1) {
    tol <- 1e-3
    warning("'tol' too high. Setting to 1e-3.")
  }
  if (tol < 1e-3) {
    warning(
      "'tol' was set below recommended value of 1e-3 which may
            result in longer convergence time."
    )
  }
  nrep <- ncol(spec)

  if (verbose) {
    cat('[INFO] Initializing population(s)...\n')
  }
  pops <- lapply(seq(nislands), function(x) {
    pop <- pop_init(popsize, nsubpop, nbands, spec)
    loss <- loss_fn(pop, spec)
    elite <- pop[[which.min(loss)]]
    elite_loss <- loss[which.min(loss)]
    return(list(pop = pop, loss = loss, elite = elite, elite_loss = elite_loss))
  })

  gen <- run <- 0
  avgloss <- minloss <- ninfeasible <- rep(0, maxgen + 1)
  avgloss[gen + 1] <- mean(unlist(lapply(pops, function(x) {
    mean(x$loss, na.rm = TRUE)
  })))
  minloss[gen + 1] <- min(unlist(lapply(pops, function(x) {
    min(x$loss, na.rm = TRUE)
  })))
  ninfeasible[gen + 1] <- sum(unlist(lapply(pops, function(x) {
    sum(is.na(x$loss))
  })))
  if (verbose) {
    cat(
      paste0(
        "[INFO] Generation: ",
        gen,
        "\tAverage Loss: ",
        round(avgloss[gen + 1], 2),
        "\tMinimum Loss: ",
        round(minloss[gen + 1], 2)
      ),
      "\n\n"
    )
  }

  if (verbose) {
    cat('[INFO] Starting evolution...\n')
  }
  while (gen < maxgen & run <= maxrun) {
    gen <- gen + 1

    # selection, mutation, and replacement; elitism
    pops <- lapply(pops, function(x) {
      probs = 1 / x$loss
      probs[is.na(probs)] <- 0
      child <- x$pop[sample(popsize, popsize, TRUE, prob = probs)]
      pop <- mutate(child, spec, pmutate)
      loss <- loss_fn(pop, spec)
      if (x$elite_loss < min(loss, na.rm = TRUE)) {
        # the previous elite is better than the current el  ite
        pop[[which.max(loss)]] <- x$elite
        loss[which.max(loss)] <- x$elite_loss
      } else {
        # the current elite is better than the previous elite
        x$elite <- pop[[which.min(loss)]]
        x$elite_loss <- loss[which.min(loss)]
      }
      return(list(
        pop = pop,
        loss = loss,
        elite = x$elite,
        elite_loss = x$elite_loss
      ))
    })

    # migration policy
    if (gen %% epoch == 0 & nislands > 1) {
      if (verbose) {
        cat('[INFO] Performing migration...\n')
      }
      ind <- sample(popsize, nmigrants, FALSE)
      samp_prev <- pops[[1]]$pop[ind]
      for (i in 2:nislands) {
        samp <- pops[[i]]$pop[ind] <- samp_prev
        samp_prev <- samp
      }
      pops[[1]]$pop[ind] <- samp_prev
    }

    avgloss[gen + 1] <- mean(unlist(lapply(pops, function(x) {
      mean(x$loss, na.rm = TRUE)
    })))
    minloss[gen + 1] <- min(unlist(lapply(pops, function(x) {
      min(x$loss, na.rm = TRUE)
    })))
    ninfeasible[gen + 1] <- sum(unlist(lapply(pops, function(x) {
      sum(is.na(x$loss))
    })))
    if (verbose & gen %% 5 == 0) {
      cat(
        paste0(
          "[INFO] Generation: ",
          gen,
          "\tAverage Loss: ",
          round(avgloss[gen + 1], 2),
          "\tMinimum Loss: ",
          round(minloss[gen + 1], 2)
        ),
        "\n"
      )
    }

    if (abs((minloss[gen + 1] - minloss[gen]) / minloss[gen]) < tol) {
      run <- run + 1
    } else {
      run <- 0
    }
  }
  if (verbose) {
    cat('[INFO] Evolution terminated.\n')
  }
  pop <- unlist(lapply(pops, function(x) x$pop), recursive = FALSE)
  loss <- unlist(lapply(pops, function(x) x$loss))
  solution <- pop[[which.min(loss)]]
  loss <- loss[which.min(loss)]

  # final assignment
  labels <- rep(1, ncol(spec))
  if (nsubpop > 1) {
    labels <- l2_assign(spec, solution$cuts, solution$avg_power)
  }

  ## selection index
  s_band <- band_sim(spec, labels, solution$cuts, solution$avg_power)
  s_pop <- NULL
  if (nsubpop > 1) {
    s_pop <- subpop_sim(spec, labels, solution$cuts, solution$avg_power)
  }

  ## output ----
  out <- list(
    spec = spec,
    cuts = solution$cuts,
    avg_power = solution$avg_power,
    labels = labels,
    loss = loss,
    s_band = s_band,
    s_pop = s_pop,
    avgloss = avgloss[1:(gen + 1)],
    minloss = minloss[1:(gen + 1)],
    ninfeasible = ninfeasible[1:(gen + 1)],
    params = list(
      nsubpop = nsubpop,
      nbands = nbands,
      popsize = popsize,
      pmutate = pmutate,
      nislands = nislands,
      epoch = epoch,
      maxgen = maxgen,
      maxrun = maxrun,
      tol = tol
    )
  )
  class(out) <- "ga"
  return(out)
}

# 'which' can be any combination of 1 or 2
#    1 - plot of average and minimum loss over generations
#    2 - plot of number of infeasible solutions over generations
plot.ga <- function(gaout, which = c(1, 2)) {
  ngen <- length(gaout$avgloss)
  if (1 %in% which) {
    plot_max <- max(gaout$avgloss, gaout$minloss)
    plot_min <- min(gaout$avgloss, gaout$minloss)
    plot(
      seq(1, ngen),
      gaout$avgloss,
      type = 'o',
      xlab = "generation",
      ylab = "loss",
      main = "Loss history",
      ylim = c(plot_min, plot_max)
    )
    lines(seq(1, ngen), gaout$minloss, type = 'o', col = 2)
    abline(v = seq(gaout$params$epoch, ngen, by = gaout$params$epoch), lty = 2)
    legend(
      "topright",
      lty = c(1, 1, 2),
      col = c(1, 2, 1),
      legend = c("average loss", "minimum loss", "migration event")
    )
  }
  if (2 %in% which) {
    plot(
      seq(1, ngen),
      gaout$ninfeasible,
      type = 'o',
      xlab = "generation",
      ylab = "number infeasible",
      main = "Infeasible solutions per generation"
    )
    abline(v = seq(gaout$params$epoch, ngen, by = gaout$params$epoch), lty = 2)
    legend(
      "topright",
      lty = 1:2,
      col = 1,
      legend = c("# infeasible sols", "migration event")
    )
  }
}

summarize_power <- function(spec, cuts) {
  nfreq <- nrow(spec)
  nrep <- ncol(spec)
  if (!1 %in% cuts) {
    cuts <- c(1, cuts)
  }
  if (!(nfreq + 1) %in% cuts) {
    cuts <- c(nfreq + 1, cuts)
  }
  cuts <- sort(cuts)
  nbands <- length(cuts) - 1
  rep <- t(do.call(
    rbind,
    lapply(seq(nbands), function(l) {
      colMeans(spec[seq(cuts[l], cuts[l + 1] - 1), , drop = F])
    })
  ))
  avg <- colMeans(rep)
  list(rep = rep, avg = avg)
}


pop_init <- function(popsize, nsubpop, nbands, spec) {
  nfreq <- nrow(spec)
  nrep <- ncol(spec)
  ncuts <- nbands - 1
  eq <- floor(seq(1, nfreq + 1, length = nbands + 1)) # equally spaced cuts
  lapply(seq(popsize), function(i) {
    avg_power <- do.call(
      rbind,
      lapply(seq(nsubpop), function(j) {
        x <- spec[, sample(nrep, 1), drop = FALSE] # random replicate
        summarize_power(x, eq)$avg
      })
    )
    list(
      cuts = matrix(rep(eq, nsubpop), nrow = nsubpop, byrow = T),
      avg_power = avg_power
    )
  })
}

mutate <- function(pop, spec, pmutate) {
  nfreq <- nrow(spec)
  nrep <- ncol(spec)
  nsubpop <- nrow(pop[[1]]$cuts)
  nbands <- ncol(pop[[1]]$avg_power)
  lapply(pop, function(p) {
    cuts <- p$cuts
    avg_power <- p$avg_power

    # cuts
    for (j in seq(nsubpop)) {
      for (l in seq(2, nbands)) {
        if (!all(diff(cuts[j, seq(l - 1, l + 1)]) == 1) & runif(1) < pmutate) {
          cuts[j, l] <- round(truncnorm::rtruncnorm(
            n = 1,
            a = cuts[j, l - 1] + 1,
            b = cuts[j, l + 1] - 1,
            mean = cuts[j, l],
            sd = max(diff(cuts[j, seq(l - 1, l + 1)])) / 8
          ))
        }
      }
    }
    cuts <- apply(cuts, 2, sort)
    if (nsubpop == 1) {
      cuts <- t(matrix(cuts))
    }

    # avg_power
    lottery <- matrix(runif(nsubpop * nbands) < pmutate, nrow = nsubpop)
    if (sum(lottery) > 0) {
      avg_power[lottery] <- truncnorm::rtruncnorm(
        n = sum(lottery),
        a = 0,
        b = max(spec),
        mean = avg_power[lottery],
        sd = sqrt(log(max(spec)))
      )
    }

    list(cuts = cuts, avg_power = avg_power)
  })
}

loss_fn <- function(pop, spec) {
  nfreq <- nrow(spec)
  nrep <- ncol(spec)
  nsubpop <- nrow(pop[[1]]$cuts)
  nbands <- ncol(pop[[1]]$collapsed) - 1
  if (nsubpop == 1) {
    assignments <- rep(1, nrep)
  }
  unlist(lapply(pop, function(p) {
    if (nsubpop > 1) {
      assignments <- l2_assign(spec, p$cuts, p$avg_power)
    }
    if (length(unique(assignments)) != nsubpop) {
      return(NA)
    }
    sum(unlist(lapply(seq(nsubpop), function(j) {
      x <- rep(p$avg_power[j, ], diff(p$cuts[j, ]))
      x <- matrix(rep(x, sum(assignments == j)), ncol = sum(assignments == j))
      sum((spec[, assignments == j] - x)^2) / nfreq
    })))
  }))
}

l2_assign <- function(spec, cuts, avg_power) {
  nfreq <- nrow(spec)
  nrep <- ncol(spec)
  nsubpop <- nrow(cuts)
  nbands <- ncol(cuts) - 1
  dist_mat <- do.call(
    rbind,
    lapply(seq(nsubpop), function(j) {
      x <- rep(avg_power[j, ], diff(cuts[j, ])) # average collapsed spectrum
      x <- matrix(rep(x, nrep), ncol = nrep)
      colSums((spec - x)^2)
    })
  )
  apply(dist_mat, 2, which.min)
}

band_sim <- function(spec, labels, cuts, avg_power) {
  nsubpop <- nrow(cuts)
  nbands <- ncol(avg_power)
  mean(unlist(lapply(seq(nsubpop), function(j) {
    mean(unlist(lapply(seq(nbands - 1), function(l) {
      # variability in band l
      v1 <- sqrt(sum(
        (spec[seq(cuts[j, l], cuts[j, l + 1] - 1), labels == j] -
          avg_power[j, l])^2
      ))

      # variability in band l+1
      v2 <- sqrt(sum(
        (spec[seq(cuts[j, l + 1], cuts[j, l + 2] - 1), labels == j] -
          avg_power[j, l + 1])^2
      ))

      # variability across bands (l, l+1)
      avg_power_across <- mean(avg_power[j, seq(l, l + 1)])
      v3 <- sqrt(sum(
        (spec[seq(cuts[j, l], cuts[j, l + 2] - 1), labels == j] -
          avg_power_across)^2
      ))

      # ratio
      (v1 + v2) / v3
    })))
  })))
}

subpop_sim <- function(spec, labels, cuts, avg_power) {
  nsubpop <- nrow(cuts)
  nbands <- ncol(avg_power)

  rmat <- matrix(0, nsubpop, nsubpop)
  for (j in seq(1, nsubpop - 1)) {
    for (i in seq(j + 1, nsubpop)) {
      yj <- rep(avg_power[j, ], diff(cuts[j, ]))
      yi <- rep(avg_power[i, ], diff(cuts[i, ]))
      vj <- sqrt(sum((spec[, labels == j] - yj)^2 / sum(labels == j)))
      vi <- sqrt(sum((spec[, labels == i] - yi)^2 / sum(labels == i)))
      dij <- sqrt(sum((yi - yj)^2))
      rmat[j, i] <- (vj + vi) / dij
    }
  }
  rmat <- rmat + t(rmat)
  sum(apply(rmat, 1, max)) / nsubpop
}

sine_mt <- function(x, ntapers = NULL) {
  x <- as.matrix(x)
  if (is.null(ntapers)) {
    ntapers <- floor(sqrt(nrow(x)))
  }
  len <- dim(x)[1]
  if (dim(x)[1] < dim(x)[2]) {
    warning(
      "Wide matrix was provided. Double check that the matrix is such
            that columns (rather than rows) correspond to replicates."
    )
  }
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

# simulation models -------------------------------------------------------

model1 <- function(nrep, len) {
  # locations of marginal breakpoints
  breaks <- matrix(
    c(
      0.1,
      0.25, #lower
      0.2,
      0.3, #middle
      0.25,
      0.4 #upper
    ),
    nrow = 3,
    byrow = TRUE
  )
  delta1 <- 0.025 # spacing for cubic spline

  # marginal values of the constant segments
  values <- matrix(
    c(
      15,
      7.5,
      2, #lower
      25,
      12.5,
      4, #middle
      35,
      17.5,
      6 #upper
    ),
    nrow = 3,
    byrow = TRUE
  )
  delta2 <- 2 # spacing around average summary measures

  # underlying spectra
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = floor(len / 2))
  spec <- matrix(nrow = floor(len / 2), ncol = 3 * nrep)
  for (j in 1:3) {
    for (k in 1:nrep) {
      spec[, nrep * (j - 1) + k] <- piecewise_smooth_step(
        freq,
        values[j, ] + runif(ncol(values), min = -delta2, max = delta2),
        breaks[j, ],
        delta1
      )
    }
  }
  x <- apply(rbind(spec, spec[dim(spec)[1]:1, ]), 2, spec_sim)
  mtout <- sine_mt(x, ntapers = floor(sqrt(len / 2)))
  return(list(
    x = x,
    labels = labels,
    freq = freq,
    spec = spec,
    mtfreq = mtout$mtfreq,
    mtspec = mtout$mtspec
  ))
}

# peaks:        vector of 3 marginal peak locations
# bw:           vector of 3 marginal bandwidths
model2 <- function(nrep, len, peaks, bw, sd = 2.25) {
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = 250) # true frequencies
  spec <- matrix(nrow = 250, ncol = 3 * nrep) # true spectra
  x <- matrix(nrow = len, ncol = 3 * nrep) # time series data
  for (i in 1:nrep) {
    # realizations of parameters using peak/bw parameterization from
    # Granados-Garcia et al. (2022)
    peaks_rep <- peaks + runif(3, -0.02, 0.02)
    bw_rep <- bw + runif(3, -0.005, 0.005)
    phi1 <- 2 * cos(2 * pi * peaks_rep) * exp(-bw_rep)
    phi2 <- -exp(-2 * bw_rep)
    # underlying spectra
    spec[, i] <- ar2_spec(freq, phi1[1], phi2[1], sd)
    spec[, nrep + i] <- ar2_spec(freq, phi1[2], phi2[2], sd)
    spec[, 2 * nrep + i] <- ar2_spec(freq, phi1[3], phi2[3], sd)
    # time series data
    x[, i] <- arima.sim(list(ar = c(phi1[1], phi2[1])), n = len, sd = sd)
    x[, nrep + i] <- arima.sim(list(ar = c(phi1[2], phi2[2])), n = len, sd = sd)
    x[, 2 * nrep + i] <- arima.sim(
      list(ar = c(phi1[3], phi2[3])),
      n = len,
      sd = sd
    )
  }
  mtout <- sine_mt(x)
  return(list(
    x = x,
    labels = labels,
    freq = freq,
    spec = spec,
    mtfreq = mtout$mtfreq,
    mtspec = mtout$mtspec
  ))
}

# model2a - most overlap between peaks
model2a <- function(nrep, len) {
  peaks <- c(0.23, 0.25, 0.27)
  bw <- rep(0.15, 3)
  return(model2(nrep, len, peaks, bw))
}

# model2b - medium overlap between peaks
model2b <- function(nrep, len) {
  peaks <- c(0.21, 0.25, 0.29)
  # bw <- c(0.155, 0.15, 0.155)
  bw <- c(0.15, 0.16, 0.2)
  return(model2(nrep, len, peaks, bw))
}

# model2c - least overlap between peaks
model2c <- function(nrep, len) {
  peaks <- c(0.19, 0.25, 0.31)
  # bw <- c(0.165, 0.15, 0.165)
  bw <- c(0.15, 0.165, 0.205)
  return(model2(nrep, len, peaks, bw))
}

model3 <- function(nrep, len) {
  freq <- seq(0, 0.5, length = 250)

  # ar1 processes
  x <- matrix(nrow = len, ncol = 3 * nrep)
  xspec <- matrix(nrow = 250, ncol = 3 * nrep)
  for (i in 1:nrep) {
    peaks_rep <- rep(0, 3)
    bw_rep <- rep(0.5, 3) + runif(3, -0.02, 0.02)
    phi1 <- 2 * cos(2 * pi * peaks_rep) * exp(-bw_rep)
    phi2 <- -exp(-2 * bw_rep)
    xspec[, i] <- ar2_spec(freq, phi1[1], phi2[1], 1)
    xspec[, nrep + i] <- ar2_spec(freq, phi1[2], phi2[2], 1)
    xspec[, 2 * nrep + i] <- ar2_spec(freq, phi1[3], phi2[3], 1)
    x[, i] <- arima.sim(list(ar = c(phi1[1], phi2[1])), n = len, sd = 1)
    x[, nrep + i] <- arima.sim(list(ar = c(phi1[2], phi2[2])), n = len, sd = 1)
    x[, 2 * nrep + i] <- arima.sim(
      list(ar = c(phi1[3], phi2[3])),
      n = len,
      sd = 1
    )
  }

  # ar2 processes
  y <- matrix(nrow = len, ncol = 3 * nrep)
  yspec <- matrix(nrow = 250, ncol = 3 * nrep)
  for (i in 1:nrep) {
    peaks_rep <- c(0.25, 0.3, 0.35) + runif(3, -0.015, 0.015)
    bw_rep <- c(0.15, 0.15, 0.15) #+ runif(3, -0.005, 0.005)
    phi1 <- 2 * cos(2 * pi * peaks_rep) * exp(-bw_rep)
    phi2 <- -exp(-2 * bw_rep)
    yspec[, i] <- ar2_spec(freq, phi1[1], phi2[1], 1)
    yspec[, nrep + i] <- ar2_spec(freq, phi1[2], phi2[2], 1)
    yspec[, 2 * nrep + i] <- ar2_spec(freq, phi1[3], phi2[3], 1)
    y[, i] <- arima.sim(list(ar = c(phi1[1], phi2[1])), n = len, sd = 2)
    y[, nrep + i] <- arima.sim(list(ar = c(phi1[2], phi2[2])), n = len, sd = 2)
    y[, 2 * nrep + i] <- arima.sim(
      list(ar = c(phi1[3], phi2[3])),
      n = len,
      sd = 2
    )
  }

  # add ar1 and ar2 processes
  labels <- rep(1:3, each = nrep)
  z <- x + y
  zspec <- xspec + yspec
  mtout <- sine_mt(z)
  return(list(
    x = z,
    labels = labels,
    freq = freq,
    spec = zspec,
    mtfreq = mtout$mtfreq,
    mtspec = mtout$mtspec
  ))
}

# piecewise smooth step function
# x         inputs to at which to evaluate function
# values    vector of length L
# breaks    vector of length L-1
# delta     double, spacing around breaks for cubic spline
piecewise_smooth_step <- function(x, values, breaks, delta = 0.025) {
  coefs <- matrix(nrow = length(breaks), ncol = 4)
  for (i in 1:length(breaks)) {
    left <- (breaks[i] - delta / 2)
    right <- (breaks[i] + delta / 2)
    A <- matrix(
      c(
        left^3,
        left^2,
        left,
        1,
        right^3,
        right^2,
        right,
        1,
        3 * left^2,
        2 * left,
        1,
        0,
        3 * right^2,
        2 * right,
        1,
        0
      ),
      nrow = 4,
      ncol = 4,
      byrow = TRUE
    )
    coefs[i, ] <- solve(A, c(values[i], values[i + 1], 0, 0))
  }
  y <- rep(values[length(values)], length(x))
  marks <- sort(c(0, breaks - (delta / 2), breaks + (delta / 2), 0.5))
  for (i in 1:length(breaks)) {
    ind <- 2 * (i - 1) + 1
    y[x >= marks[ind] & x < marks[ind + 1]] <- values[i]
    x_gap <- x[x >= marks[ind + 1] & x < marks[ind + 2]]
    y[x %in% x_gap] <- matrix(
      c(
        x_gap^3,
        x_gap^2,
        x_gap,
        rep(1, length(x_gap))
      ),
      ncol = 4
    ) %*%
      coefs[i, ]
  }
  return(y)
}

# ar1 spectrum
ar1_spec <- function(x, phi, sd) {
  sd^2 / (1 + phi^2 - 2 * phi * cos(2 * pi * x))
}

# ar2 spectrum
ar2_spec <- function(x, phi1, phi2, sd) {
  sd^2 /
    (1 +
      phi1^2 +
      phi2^2 -
      2 * phi1 * (1 - phi2) * cos(2 * pi * x) -
      2 * phi2 * cos(4 * pi * x))
}

# simulate time series realization from a single theoretical spectrum
# See: Guo and Dai (2006) Multivariate time-dependent spectral analysis using
# Cholesky decomposition
# spec - vector of spectrum values
spec_sim <- function(spec) {
  nx <- length(spec)
  sd <- sqrt(1 / (2 * nx))
  z <- vector("complex", nx)
  y <- vector("complex")
  i <- complex(imaginary = 1)
  for (j in 1:nx) {
    if (j / nx == 0.5 || j / nx == 1) {
      z[j] <- rnorm(1, sd = sd)
    } else if (j < floor(nx / 2) + 1) {
      z[j] <- complex(real = rnorm(1, sd = sd), imaginary = rnorm(1, sd = sd))
    } else {
      z[j] <- Conj(z[nx - j])
    }
  }
  for (t in 1:nx) {
    y[t] <- sum(sqrt(spec) * exp(2 * pi * i * (1:nx) * t / nx) * z)
  }
  return(Re(y))
}
