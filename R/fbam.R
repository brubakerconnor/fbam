library(doParallel)

#### frequency band analysis for multiple time series ##########################
fbam <- function(X, nclust_grid, nbands_grid, parallel = TRUE, popsize = 50,
                 nparents = 2 * popsize, nelite = 2, pcrossover = 0.5,
                 pmutate = 0.15, endpoint_range = NULL, collapsed_range = NULL,
                 maxgen = 500, maxrun = 50, tol = 1e-3, s = NULL,
                 ntapers = floor(sqrt(dim(X)[1])), verbose = TRUE) {

  # enable parallelization, if specified
  if (parallel & !foreach::getDoParRegistered()) {
    num_cores <- parallel::detectCores() - 1
    cat("Registering a cluster with ", num_cores, " workers.\n")
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function() {
      library(doParallel); source('R/fbam.R')
    })
  } else if (foreach::getDoParRegistered()) {
    message("Using previously registered cluster.")
  } else {
    message("Parallelization not specified. Optimizing sequentially.")
  }

  nclust_grid <- unique(nclust_grid)
  nbands_grid <- unique(nbands_grid)
  nclust_grid_size <- length(nclust_grid)
  nbands_grid_size <- length(nbands_grid)

  nfreq <- dim(X)[1]; nrep <- dim(X)[2]
  if (any(nclust_grid < 1) | any(nclust_grid > nrep)) {
    stop("Ensure valid nclust grid.")
  }
  if (any(nbands_grid < 1) | any(nbands_grid > nfreq - 1)) {
    stop("Ensure valid nbands grid.")
  }

  # optimize objective over grid
  # `%:%` <- foreach::`%:%`
  # `%dopar%` <- foreach::`%dopar%`
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
  # parallel::stopCluster(cl)

  # compute selection index for each solution to find "best" solution
  # if more than one pair of values for clusters and bands were provided
  if (nclust_grid_size + nbands_grid_size > 2) {
    current_best <- grid[[1]][[1]]
    current_best_si <- grid[[1]][[1]]$si
    for (j in 1:nclust_grid_size) {
      for (b in 1:nbands_grid_size) {
        current_out <- grid[[j]][[b]]
        current_out_si <- grid[[j]][[b]]$si
        if (current_out_si < current_best_si) {
          current_best <- current_out; current_best_si <- current_out_si
        }
      }
    }
    return(list(selected = current_best, grid = grid))
  }
  return(grid[[1]][[1]])
}

#### genetic algorithm #########################################################
ga <- function(X, nclust, nbands, popsize = 50, pop = NULL,
               nparents = 2 * popsize, nelite = 2, pcrossover = 0.5,
               pmutate = 0.15, endpoint_range = NULL, collapsed_range = NULL,
               maxgen = 500, maxrun = 50, tol = 1e-3, s = NULL,
               ntapers = floor(sqrt(dim(X)[1])), verbose = F) {
  # estimate spectra from X
  if (is.vector(X)) X <- as.matrix(X)
  mtout <- sine_multitaper(X, ntapers = ntapers)
  mtfreq <- mtout$mtfreq
  spec <- mtout$mtspec
  nfreq <- dim(spec)[1]

  # ensure that the number of parents is at least 2
  if (nparents < 2) {
    warning(paste0("Number of parents specified (nparents = ", nparents,
                   ") was too small. nparents = 2 will be used."))
    nparents <- 2
  }

  # mutation parameters default to adaptive values if none are specified
  # endpoint range
  if (is.null(endpoint_range)) {
    # monotonically decreasing endpoint mutation range over generations
    endpoint_range <- floor(seq(from = floor(0.1 * nrow(spec)),
                                to = 1,
                                length = maxgen))
  } else if (is.integer(endpoint_range)) {
    # use same value for each generation
    endpoint_range <- rep(endpoint_range, maxgen)
  } else if (is.vector(endpoint_range)) {
    if (length(endpoint_range) < maxgen) {
      stop(paste0("If providing a vector of endpoint mutation ranges for
                  adaptive mutation, ensure the length of this vector is
                  at least the maximum number of allowed generations (",
                  maxgen, ")."))
    }
  }

  # collapsed measure range
  if (is.null(collapsed_range)) {
    # monotonically decreasing endpoint mutation range over generations
    collapsed_range <- seq(from = 0.2 * max(spec),
                           to = 0.01 * sd(spec),
                           length = maxgen)
  } else if (is.double(collapsed_range)) {
    # use same value for each generation
    collapsed_range <- rep(collapsed_range, maxgen)
  } else if (is.vector(collapsed_range)) {
    if (length(collapsed_range) < maxgen) {
      stop(paste0("If providing a vector of collapsed measure mutation ranges
                  for adaptive mutation, ensure the length of this vector is
                  at least the maximum number of allowed generations (",
                  maxgen, ")."))
    }
  }

  # checks on maxgen and maxrun
  if (maxgen < 1) {
    stop("maxgen must be at least 1. Recommended setting is 500.")
  }
  if (maxrun < 1) {
    stop("maxrun must be at least 1. Recommended setting is 100")
  }
  if (maxrun > maxgen) {
    maxrun <- maxgen
    warning("maxrun greater than maxgen so convergence will be declared only
            when mxagen number of generations have been completed.
            Recommended setting is maxrun = 100.")
  }
  if (tol > 1) {
    tol <- 1
    warning("Tolerance set to 1. Recommended setting is 1e-3.")
  }
  if (tol < 1e-3) {
    warning("Tolerance was set below recommended value of 1e-3 which may
            cause longer than expected convergence times.")
  }

  # set s
  if (is.null(s)) s <- floor(0.1 * dim(spec)[1])

  # STEP 1 - population initialization
  if (popsize < 1) {
    stop("popsize must be at least 1. Recommended setting is 50.")
  }
  if (is.null(pop)) {
    pop <- pop_init(popsize, nclust, nbands, spec, s)
  } else {
    popsize <- length(pop)
  }

  # STEP 2 - initial evaluation of population
  fitness_values <- evaluate(pop, spec)
  if (nelite > 0) {
    # rank the population by fitness value
    u <- rank(-fitness_values, ties.method = "random")

    # get the top nelite ranked members and their fitness values
    elites <- pop[which(u %in% seq_len(nelite))]
    elite_fitness_values <- fitness_values[which(u %in% seq_len(nelite))]
  }

  # initialize monitoring
  gennum <- runlength <- 0
  avgfit <- maxfit <- ninfeasible <- rep(0, maxgen + 1)
  avgfit[gennum + 1] <- mean(fitness_values)
  maxfit[gennum + 1] <- max(fitness_values)
  ninfeasible[gennum + 1] <- sum(fitness_values == 1e-50)
  if (verbose) {
    cat(paste0(
      "START \t\t AVG FIT: ", round(avgfit[gennum + 1], 5),
      "\t\t MAX FIT: ", round(maxfit[gennum + 1], 5)
    ), "\n")
  }

  # iterate until convergence
  while (gennum < maxgen & runlength <= maxrun) {
    gennum <- gennum + 1 # increase generation count

    # parent selection by FPS via SUS
    parents <- parent_selection(pop, nparents, fitness_values)

    # crossover
    if (pcrossover > 0) {
      # number of pairs
      nmating <- floor(nparents / 2)
      # randomly pair the members of the current population
      mating <- matrix(sample(1:(2 * nmating), size = 2 * nmating), ncol = 2)
      for (i in seq_len(nmating)) {
        if (runif(1) < pcrossover) {
          # perform mutation (otherwise they are left alone)
          cout <- crossover(parents[[mating[i, 1]]], parents[[mating[i, 2]]])
          parents[[mating[i, 1]]] <- cout$c1
          parents[[mating[i, 2]]] <- cout$c2
        }
      }
    }

    # mutation
    for (i in seq_len(nparents)) {
      # mutate the current member
      mutated <- mutate(parents[[i]],
                        pmutate,
                        endpoint_range[gennum],
                        collapsed_range[gennum],
                        nfreq)
      # cluster the spectra on the basis of this mutated solution
      # if the mutated solution results in an infeasible solution,
      # the mutations are discarded and the current solution is left as is
      mutated_labels <- apply(dist_mat(spec, mutated), 1, which.min)
      if (length(unique(mutated_labels)) == nclust) parents[[i]] <- mutated
    }

    # evaluation and rank selection of offspring
    # evaluate the children that were created via crossover and/or mutation
    # of the selected parents.
    # children are evaluated and ranked. only the top popsize of them are
    # allowed to survive into the next generation.
    # the previous generation from which parents were selected is discarded.
    fitness_values <- evaluate(parents, spec)
    u <- rank(-fitness_values, ties.method = "random")
    pop <- parents[u <= popsize]
    fitness_values <- fitness_values[u <= popsize]

    # elitism (if specified)
    if (nelite > 0) {
      # rank the current population and the elites of the previous generation
      u <- rank(-c(fitness_values, elite_fitness_values),
                ties.method = "random")

      # use rank selection to permit elites that are better than at least one
      # of the children selected above to be carried over
      pop <- c(pop, elites)[u <= popsize]
      fitness_values <- c(fitness_values, elite_fitness_values)[u <= popsize]

      # re-rank the new population which may contain the elites from the
      # previous generation to get the new elites of the current generation
      u <- rank(-c(fitness_values, elite_fitness_values),
                ties.method = "random")
      elites <- c(pop, elites)[which(u %in% seq_len(nelite))]
      elite_fitness_values <- c(fitness_values, elite_fitness_values)[which(u %in% seq_len(nelite))]
    }

    # monitoring
    avgfit[gennum + 1] <- mean(fitness_values)
    maxfit[gennum + 1] <- max(fitness_values)
    ninfeasible[gennum + 1] <- sum(fitness_values == 1e-50)
    prevmfit <- maxfit[gennum]
    currmfit <- maxfit[gennum + 1]
    if ((currmfit - prevmfit) / prevmfit < tol) {
      runlength <- runlength + 1
    } else {
      runlength <- 0
    }
    if (verbose & gennum %% 25 == 0) {
      cat(paste0(
        "GEN ", gennum, "\t\t AVG FIT: ", round(avgfit[gennum + 1], 5),
        "\t\t MAX FIT: ", round(maxfit[gennum + 1], 5)
      ), "\n")
    }
  }

  # determine solution
  solution <- pop[[which.max(fitness_values)]]
  labels <- apply(dist_mat(spec, solution), 1, which.min)
  endpoints_index <- solution[, 1:(nbands - 1), drop = FALSE]
  endpoints <- endpoints_index
  for (j in 1:nclust) endpoints[j, ] <- mtfreq[endpoints_index[j, ]]
  collapsed <- avg_collapsed_by_index(spec, endpoints_index, labels)

  # compute replicate specific collapsed measures
  rep_collapsed <- rep_collapsed_by_index(spec, endpoints_index, labels)

  # compute internal validation index
  si <- selection_index(spec, labels, endpoints_index, collapsed)

  return(list(
    X = X,
    spec = spec,
    labels = labels,
    endpoints_index = endpoints_index,
    endpoints = endpoints,
    collapsed = collapsed,
    rep_collapsed = rep_collapsed,
    objective = 1 / fitness_values[which.max(fitness_values)],
    avgfit = avgfit[1:(gennum + 1)],
    maxfit = maxfit[1:(gennum + 1)],
    ninfeasible = ninfeasible[1:(gennum + 1)],
    si = si,
    params = list(
      nclust = nclust,
      nbands = nbands,
      popsize = popsize,
      nparents = nparents,
      nelite = nelite,
      pcrossover = pcrossover,
      pmutate = pmutate,
      endpoint_range = endpoint_range,
      collapsed_range = collapsed_range,
      maxgen = maxgen,
      maxrun = maxrun,
      tol = tol
    )
  ))
}

#### collapsed measures of power ###############################################
rep_collapsed_by_index <- function(spec, endpoints_index, labels = NULL) {
  spec <- as.matrix(spec)
  nfreq <- dim(spec)[1]
  nrep <- dim(spec)[2]
  if (is.vector(endpoints_index)) {
    endpoints_index <- matrix(endpoints_index, nrow = 1, byrow = T)
  }
  nclust <- nrow(endpoints_index)

  # if labels were provided, ensure that an appropriate matrix of endpoints
  # was provided, too
  if (!is.null(labels)) {
    if (length(labels) < nrep) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
                  "not enough labels provided."))
    }
    if (length(labels) > nrep) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
                  "too many labels provided."))
    }
    if (nclust < length(unique(labels))) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
                  "labels were provided but not enough sets of frequency bands were ",
                  "provided; ", nclust, " sets were given."))
    }
    if (any(unique(labels) < 1) | any(unique(labels) > nclust)) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
                  "some labels are outside the range 1:nclust."))
    }
  }

  # if no labels were provided but multiple sets of endpoints were provided,
  # compute collapsed measures using the first set and send a warning to user
  if (is.null(labels) & nclust == 1) {
    labels <- rep(1, nrep)
  }
  if (is.null(labels) & (nclust > 1)) {
    warning("Warning in computing collapsed measures by index; no labels were ",
            "provided, but multiple sets of endpoints were provided. Only the ",
            "first set of endpoints were used.")
    endpoints_index <- endpoints_index[1, , drop = FALSE]
    nclust <- 1
  }

  # check that the endpoints are in valid range
  if (any(endpoints_index > (nfreq + 1)) | any(endpoints_index < 1)) {
    stop(paste(c("Error in computing rep. collapsed measures by index; ",
                 "endpoints must be between 1 and nrow(spec) but\n",
                 endpoints_index, "\nwas provided."), collapse = " "))
  }

  # ensure the endpoints include 1 and (nfreq + 1) before continuing
  if (all(endpoints_index[, 1] != 1)) {
    endpoints_index <- cbind(rep(1, nclust), endpoints_index)
  } else if (any(endpoints_index[, 1] != 1)) {
    stop(paste(c("Error in computing rep. collapsed measures by index; ",
                 "ensure that the first provided endpoint of each cluster is 1 ",
                 "or the right endpoint of the first band. Got\n",
                 endpoints_index, "\n instead."), collapse = " "))
  }
  if (all(endpoints_index[, ncol(endpoints_index)] != (nfreq + 1))) {
    endpoints_index <- cbind(endpoints_index, rep(nfreq + 1, nclust))
  } else if (any(endpoints_index[, ncol(endpoints_index)] != (nfreq + 1))) {
    stop(paste(c("Error in computing rep. collapsed measures by index; ",
                 "ensure that the last provided endpoint of each cluster is ",
                 "nfreq + 1 or the left endpoint of the last band. Got\n",
                 endpoints_index, "\n instead."), collapse = " "))
  }
  for (j in 1:nclust) {
    # ensure endpoints are sorted
    endpoints_index[j,] <- sort(endpoints_index[j,])
  }
  nbands <- ncol(endpoints_index) - 1

  # compute and return collapsed measures
  collapsed <- matrix(nrow = nrep, ncol = nbands)
  for (j in 1:nclust) {
    for (b in 1:nbands) {
      left <- endpoints_index[j, b]
      right <- endpoints_index[j, b + 1] - 1
      band_spec <- spec[left:right, labels == j, drop = F]
      collapsed[labels == j, b] <- colMeans(band_spec)
    }
  }
  return(collapsed)
}

avg_collapsed_by_index <- function(spec, endpoints_index, labels = NULL) {
  spec <- as.matrix(spec)
  nrep <- dim(spec)[2]
  if (is.vector(endpoints_index)) {
    endpoints_index <- matrix(endpoints_index, nrow = 1, byrow = T)
  }
  nclust <- nrow(endpoints_index)

  rep_collapsed <- rep_collapsed_by_index(spec, endpoints_index, labels)
  avg_collapsed <- matrix(nrow = nclust, ncol = ncol(rep_collapsed))
  if (is.null(labels)) labels <- rep(1, nrep)
  for (j in 1:nclust) {
    clust_rep_collapsed <- rep_collapsed[labels == j, , drop = FALSE]
    avg_collapsed[j,] <- colMeans(clust_rep_collapsed)
  }
  return(avg_collapsed)
}

#### objective function and population evaluation ##############################
objective_function <- function(ch, labels, spec) {
  nclust <- dim(ch)[1]; nbands <- (dim(ch)[2] + 1) / 2
  nfreq <- dim(spec)[1]

  # set aside the endpoints from the solution
  endpoints_index <- ch[, 1:(nbands - 1), drop = FALSE]
  within_cluster_ss <- vector("double", nclust)
  for (j in 1:nclust) {
    # set aside cluster spectra and endpoints
    clust_spec <- spec[, labels == j, drop = FALSE]
    endpoints_index_clust <- endpoints_index[j, ]

    # get collapsed measures and widths in each band
    collapsed <- avg_collapsed_by_index(clust_spec, endpoints_index_clust)
    nfreq_per_band <- diff(c(1, endpoints_index_clust, nfreq + 1))

    # expand the collapsed measures and compute L2 distance
    center <- rep(collapsed, nfreq_per_band)
    within_cluster_ss[j] <- sum((clust_spec - center)^2) / (2 * (nfreq + 1))
  }
  return(sum(within_cluster_ss))
}

evaluate <- function(pop, spec) {
  popsize <- length(pop)
  fitness_values <- rep(0, popsize)
  for (i in seq_len(popsize)) {
    fitness_values[i] <- fitness(pop[[i]], spec)
  }
  return(fitness_values)
}

fitness <- function(ch, spec) {
  nclust <- dim(ch)[1]

  # assign time series to clusters and compute fitness accordingly
  labels <- apply(dist_mat(spec, ch), 1, which.min)
  if (length(unique(labels)) != nclust) {
    # penalize if there empty clusters
    return(1e-50)
  } else {
    # compute the objective function as usual
    return(1 / objective_function(ch, labels, spec))
  }
}

dist_mat <- function(spec, ch) {
  nclust <- dim(ch)[1]; nbands <- (dim(ch)[2] + 1) / 2
  nfreq <- dim(spec)[1]; nrep <- dim(spec)[2]

  # separate out endpoints and collapsed measures
  endpoints_index <- ch[, 1:(nbands - 1), drop = FALSE]
  collapsed <- ch[, nbands:dim(ch)[2], drop = FALSE]

  # compute distance matrix
  dist <- matrix(nrow = nrep, ncol = nclust)
  for (j in seq_len(nclust)) {
    # set aside cluster endpoints and collapsed measures
    endpoints_index_clust <- endpoints_index[j, ]
    collapsed_clust <- collapsed[j, ]

    # get widths of each band and expand collapsed measures
    nfreq_per_band <- diff(c(1, endpoints_index_clust, nfreq + 1))
    center <- rep(collapsed_clust, nfreq_per_band)

    # compute L2 distance between each spectrum and the jth cluster
    dist[, j] <- colSums((spec - center)^2) / (2 * (nfreq + 1))
  }
  return(dist)
}

#### population initialization #################################################
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

#### genetic operators #########################################################
parent_selection <- function(pop, nparents, fitness_values) {
  p <- cumsum(fitness_values / sum(fitness_values))
  s <- stochastic_universal_sampling(p, nparents)
  return(pop[s])
}

#' Stochastic Universal Sampling (SUS)
#'
#' This particular implementation is based on pseudocode in "Introduction to
#' Evolutionary Computing", 2e by Eiben and Smith (2015).
#'
#' @param cdf cumulative distribution function.
#' @param n number of items to select.
#'
#' @return integer indices of selected items.
stochastic_universal_sampling <- function(cdf, n) {
  cr <- i <- 1
  r <- runif(1, min = 0, max = 1 / n); s <- vector("integer", n)
  while (cr <= n) {
    while (r <= cdf[i]) {
      s[cr] <- i; r <- r + (1 / n); cr <- cr + 1
    }
    i <- i + 1
  }
  return(s)
}

mutate <- function(ch, pmutate, endpoint_range, collapsed_range, nfreq) {
  nclust <- dim(ch)[1]; nbands <- (dim(ch)[2] + 1) / 2

  # separate out endpoints and collapsed measures for this solution
  endpoints_index <- ch[, 1:(nbands - 1), drop = FALSE]
  collapsed <- ch[, nbands:dim(ch)[2], drop = FALSE]

  # endpoint mutation
  for (j in 1:nclust) {
    # get the endpoints for this cluster
    endpoints_index_clust <- endpoints_index[j, ]

    # pad the existing endpoints with 1 and `nfreq + 1` to get the
    # appropriate windows for mutation
    endpoints_index_new <- sort(c(1, endpoints_index_clust, nfreq + 1))

    # mutate each endpoint
    for (b in 2:(length(endpoints_index_new) - 1)) {
      # ensure that the current neighboring endpoints to the left and right
      # are not immediately to the left and right; otherwise, proceed with
      # mutation of this endpoint with probability pmutate
      if (!all(diff(endpoints_index_new[b + -1:1]) == 1) & runif(1) < pmutate) {
        # get the left and right endpoints for possible mutation
        # takes into account the existing endpoints so that the new endpoint
        # will not be selected outside these existing endpoints if they are
        # included in the range of mutation specified by endpoint_range
        left <- max(endpoints_index_new[b - 1] + 1,
                    endpoints_index_new[b] - endpoint_range)
        right <- min(endpoints_index_new[b + 1] - 1,
                     endpoints_index_new[b] + endpoint_range)

        # randomly sample the new endpoint from this window
        endpoints_index_new[b] <- sample(seq(left, right), 1)

        # I think this was included because I kept seeing an error and this
        # fixed it
        if(endpoints_index_new[b] == 1) {
          endpoints_index_new[b] <- 2
        }
      }
    }
    # replace existing (internal) endpoints with the new (internal) endpoints
    # sorting for a sanity check to ensure no errors later on in the algorithm
    new_internal <- endpoints_index_new[2:(length(endpoints_index_new) - 1)]
    endpoints_index[j, ] <- sort(new_internal)
  }

  # collapsed measure mutation
  # determine which of the existing collapsed measures should be mutated
  # each collapsed measure is mutated with probability pmutate independently
  # of the other measures
  s <- matrix(runif(nclust * ncol(collapsed)) < pmutate, nrow = nclust)
  collapsed[s] <- collapsed[s] +
    runif(sum(s), min = -(collapsed_range / 2), max = collapsed_range / 2)
  collapsed[collapsed < 0] <- 0 # ensure new measures are non-negative

  # return the mutated and re-combined solution
  return(cbind(endpoints_index, collapsed))
}

crossover <- function(p1, p2) {
  nclust <- dim(p1)[1]

  sl <- runif(nclust) < 0.5 # which clusters should be swapped?
  # swap the clusters as determined above between p1 and p2
  c1 <- rbind(p1[seq_len(nclust)[sl], ], p2[seq_len(nclust)[!sl], ])
  c2 <- rbind(p1[seq_len(nclust)[!sl], ], p2[seq_len(nclust)[sl], ])
  return(list(c1 = c1, c2 = c2))
}

#### selection index (joint) ###################################################
selection_index <- function(spec, labels, endpoints_index, collapsed) {
  nclust <- dim(endpoints_index)[1]
  band_index <- avg_global_band_similarity(spec, labels,
                                           endpoints_index, collapsed)
  if (nclust > 1) {
    clust_index <- avg_clust_similarity(spec, labels,
                                        endpoints_index, collapsed)
  } else {
    clust_index <- 0
  }
  return(clust_index + band_index)
}

#### selection index (no. of clusters) #########################################
clust_dist <- function(endpoints_index1, endpoints_index2, collapsed1,
                       collapsed2) {
  # ensure sorted endpoint indices
  # it is assumed these endpoints include 1 and nfreq + 1
  endpoints_index1 <- sort(endpoints_index1)
  endpoints_index2 <- sort(endpoints_index2)

  # expand the collapsed measures to full collapsed spectra
  # and compute L2 distance
  expanded1 <- rep(collapsed1, diff(endpoints_index1))
  expanded2 <- rep(collapsed2, diff(endpoints_index2))
  nfreq <- length(expanded1)
  sqrt(sum((expanded1 - expanded2)^2))
}

clust_sd <- function(clust_spec, endpoints_index, collapsed) {
  nfreq <- dim(clust_spec)[1]; nrep <- dim(clust_spec)[2]

  # ensure sorted endpoint indices
  # it is assumed these endpoints include 1 and nfreq + 1
  endpoints_index <- sort(endpoints_index)

  # expand the collapsed measures to full collapsed spectra
  # and compute the sd of the L2 distance around this object
  expanded <- rep(collapsed, diff(endpoints_index))
  sqrt(sum((clust_spec - expanded)^2 / nrep))
}

clust_similarity <- function(clust_spec1, clust_spec2,
                             endpoints_index1, endpoints_index2,
                             collapsed1, collapsed2) {
  # cluster standard deviations of distances to collapsed spectra
  sd1 <- clust_sd(clust_spec1, endpoints_index1, collapsed1)
  sd2 <- clust_sd(clust_spec2, endpoints_index2, collapsed2)

  # distance between cluster-specific collapsed spectra
  d <- clust_dist(endpoints_index1, endpoints_index2, collapsed1, collapsed2)

  # similarity is ratio of sum of std devs to distance
  return((sd1 + sd2) / d)
}

avg_clust_similarity <- function(spec, labels, endpoints_index, collapsed) {
  # allows for the possibility of empty clusters though this does not
  # happen often with the design of the GA
  ulabels <- unique(labels)
  n_nonempty_clust <- length(ulabels)
  nclust <- dim(endpoints_index)[1]; nfreq <- dim(spec)[1]

  # compute the similarity between each pair of clusters
  # the resulting matrix is symmetric so we only need to compute the lower
  # triangle of the matrix
  sim_mat <- matrix(0, nrow = nclust, ncol = nclust)
  for (i in ulabels) {
    # get the spectra, endpoint indices, and collapsed measures for cluster i
    # endpoints are padded with 1 and nfreq + 1 before computing similarity
    clust_spec_i <- spec[, labels == i, drop = F]
    endpoints_index_i <- c(1, endpoints_index[i, ], nfreq + 1)
    collapsed_i <- collapsed[i, ]

    for (j in ulabels) {
      # get the spectra, endpoint indices, and collapsed measures for cluster j
      # endpoints are padded with 1 and nfreq + 1 before computing similarity
      clust_spec_j <- spec[, labels == j, drop = F]
      endpoints_index_j <- c(1, endpoints_index[j, ], nfreq + 1)
      collapsed_j <- collapsed[j, ]

      # compute similarity if still in the lower triangular region of the
      # similarity matrix
      if (j < i) {
        sim_mat[i, j] <- clust_similarity(
          clust_spec_i, clust_spec_j,
          endpoints_index_i, endpoints_index_j,
          collapsed_i, collapsed_j
        )
      }
    }
  }
  # make matrix symmetric before computing maximum in each row
  sim_mat <- sim_mat + t(sim_mat)
  # cluster similarity index is average of the similarity of each cluster with
  # its most similar cluster
  return(sum(apply(sim_mat, 1, max)) / n_nonempty_clust)
}

#### selection index (no. of bands) ############################################
unioned_band_sd <- function(union_band_spec, nfreq) {
  union_collapsed <- mean(union_band_spec) # unioned collapsed measure of power
  sum((union_band_spec - union_collapsed)^2)
}

band_sd <- function(band_spec, band_collapsed, nfreq) {
  sum((band_spec - band_collapsed)^2)
}

band_similarity <- function(band_spec1, band_spec2, collapsed1, collapsed2,
                            nfreq) {
  # compute single band SDs
  sd1 <- band_sd(band_spec1, collapsed1, nfreq)
  sd2 <- band_sd(band_spec2, collapsed2, nfreq)

  # unioned band SD
  union_band_spec <- rbind(band_spec1, band_spec2)
  d <- unioned_band_sd(union_band_spec, nfreq)
  return((sd1 + sd2) / d)
}

avg_clust_band_similarity <- function(clust_spec, endpoints_index, collapsed) {
  nbands <- length(collapsed); nfreq <- dim(clust_spec)[1]

  # compute similarity of each band with its right neighboring band
  sim_vec <- rep(0, nbands - 1)
  for (b in 1:(nbands - 1)) {
    # endpoints, band spectra and collapsed measures for left band
    left_endpoint1 <- endpoints_index[b]
    right_endpoint1 <- endpoints_index[b + 1]
    band_spec1 <- clust_spec[left_endpoint1:(right_endpoint1 - 1), , drop = F]
    collapsed1 <- collapsed[b]

    # endpoints, band spectra and collapsed measures for right band
    left_endpoint2 <- endpoints_index[b + 1]
    right_endpoint2 <- endpoints_index[b + 2]
    band_spec2 <- clust_spec[left_endpoint2:(right_endpoint2 - 1), , drop = F]
    collapsed2 <- collapsed[b + 1]

    # similarity between left and right bands
    sim_vec[b] <- band_similarity(band_spec1, band_spec2,
                                  collapsed1, collapsed2,
                                  nfreq)
  }
  # compute average similarity over nbands
  return(sum(sim_vec) / (nbands - 1))
}

avg_global_band_similarity <- function(spec, labels, endpoints_index, collapsed) {
  ulabels <- unique(labels)
  n_nonnempty_clust <- length(ulabels)
  nclust <- dim(endpoints_index)[1]; nfreq <- dim(spec)[1]

  sim_vec <- rep(0, nclust)
  for (i in ulabels) {
    # set aside spectra, endpoints (which are then padded), and collapsed
    # measures for cluster i
    clust_spec <- spec[, labels == i, drop = F]
    endpoints_index_i <- c(1, endpoints_index[i, ], nfreq + 1)
    collapsed_i <- collapsed[i, ]

    # average band similarity for cluster i
    sim_vec[i] <- avg_clust_band_similarity(clust_spec,
                                            endpoints_index_i,
                                            collapsed_i)
  }
  return(mean(sim_vec))
}

#### multitaper estimation #####################################################
sine_multitaper <- function(X, ntapers = floor(sqrt(dim(as.matrix(X))[1]))) {
  if (is.vector(X)) X <- as.matrix(X)
  n <- dim(X)[1]
  if (n < dim(X)[2]) {
    warning("Wide matrix was provided. Ensure that the matrix is formatted such
            that columns (rather than rows) correspond to individual
            realizations.")
  }

  mtfreq <- seq(from = 1 / n, by = 1 / n, length = floor(n / 2) - 1)
  mtspec <- apply(X, 2, function(y) sine_multitaper_single(y, ntapers))

  return(list(
    mtfreq = mtfreq,
    mtspec = mtspec
  ))
}

sine_multitaper_single <- function(x, ntapers = floor(sqrt(length(x)))) {
  n <- length(x) # number of observations
  m <- floor(n / 2) - 1 # number of Fourier frequencies excluding 0 and n / 2

  sine_tapers <- outer(seq_len(ntapers), seq_len(n), function(x, y) {
    sqrt(2 / (n + 1)) * sin(pi * y * x / (n + 1))
  })
  direct_spec <- apply(sine_tapers, 1, function(r) Mod(fft(r * x)[2:(m + 1)])^2)
  return(rowMeans(direct_spec))
}

#### plotting functions ########################################################
plot_ga <- function(out) {
  numgen <- length(out$avgfit)
  df <- data.frame(gen = 1:numgen,
                   avgobjective = 1 / out$avgfit, minobjective = 1 / out$maxfit,
                   ninfeasible = out$ninfeasible)
  plot(df$gen, df$avgobjective,
       type = 'l', col = 1, lwd = 1.5,
       ylim = c(min(df$minobjective), max(df$avgobjective)),
       xaxt = 'n',
       xlab = "Generation Number", ylab = "Objective Function",
       main = "GA Evolution History")
  axis(1, at = floor(seq(0, numgen, length = 10)))
  grid()
  lines(df$gen, df$minobjective, type = 'l', col = 2, lwd = 1.5)
  legend("topright", legend = c("Average objective", "Minimum objective"),
         lwd = rep(1.5, 2), col = 1:2)
}

plot_ga_solution <- function(out) {
  n <- dim(out$X)[1]; nrep <- dim(out$X)[2]
  spec <- out$spec; nfreq <- dim(spec)[1]
  mtfreq <- seq(from = 1/n, by = 1/n, length = nfreq)

  nclust <- out$params$nclust; nbands <- out$params$nbands
  endpoints <- out$endpoints; collapsed <- out$collapsed
  endpoints_index <- endpoints * dim(out$X)[1]
  nendpoints <- nbands - 1

  if (nclust == 2) old.par <- par(mfrow = c(1, 2))
  if (nclust > 2) old.par <- par(mfrow = c(ceiling(nclust / 3), 3))
  for (j in 1:nclust) {
    matplot(mtfreq,
            spec[, out$labels == j, drop = F],
            type = "l", lty = 1, lwd = 0.2, col = "black",
            ylim = c(0, max(spec)),
            xlab = "Frequency", ylab = "Power",
            main = paste0("Cluster Label ", j)
    )
    abline(v = endpoints[j, ], col = 1, lwd = 2, lty = 1)
  }
  if (nclust >= 2) par(old.par)
}
