# command line arguments
# [1] model_name (e.g., "model1", "model2a", etc.)
# [2] nrep - number of replicate time series to use
# [3] len - length of the time series
# [4] nsim - number of simulations (e.g., 100)
# [5] RESULTS_DIR - where to save the simulation results
# [6] ncores - number of cores to use in parallelization
# [7] nclust - oracle number of clusters since no cluster validation is
# performed in this analysis
args <- commandArgs(trailingOnly = T)
model_name <- args[1]
nrep <- as.integer(args[2]); len <- as.integer(args[3]); nsim <- as.integer(args[4])
RESULTS_DIR <- args[5]
ncores <- as.integer(args[6])
nclust <- as.integer(args[7])

# seed for reproducibility
set.seed(541)
library(fda)

# confirm settings in the log
cat("Running BMARD simulation study with the following parameters:\n",
    "\tmodel_name\t\t=", model_name, "\n",
    "\tnrep\t\t\t=", nrep, "\n",
    "\tlen\t\t\t=", len, "\n",
    "\tnsim\t\t\t=", nsim, "\n",
    "\tRESULTS_DIR\t\t=", RESULTS_DIR, "\n",
    "\tncores\t\t\t=", ncores, "\n",
    "\tnclust\t\t\t=", nclust, "\n")

# source required functions
library(mclust)
Rcpp::sourceCpp("comparison-study/BMARD/CPPcode/BMARD_V112020.cpp")
source("comparison-study/BMARD/Rcode_auxiliary/ExtractionBMARDmaincomponentsmodes.R")
source('R/fbam.R') # functions for optimization and collapsed measures
source("R/generate_sim_data.R")
source("comparison-study/BMARD_comparison_functions.R")

# path to the results data file
# remove trailing separator if needed; only works on UNIX machine
if (substr(RESULTS_DIR, nchar(RESULTS_DIR), nchar(RESULTS_DIR)) == "/") {
  RESULTS_DIR <- substr(RESULTS_DIR, 1, nchar(RESULTS_DIR) - 1)
}
DATA_FILE_NAME <- file.path(RESULTS_DIR,
                            paste0(model_name, "_nrep=", nrep, "_len=", len,
                                   ".rda"))
cat("Results will be saved at", DATA_FILE_NAME, "\n")

# run simulation nsim times and save results to disk at each iteration
start <- Sys.time()
sim_data <- list()
for (n in 1:nsim) {
  tryCatch({
    # generate data
    cat("Simulation repetition n =", n, "of", nsim, "starting at",
        format(Sys.time(), usetz = TRUE), "\n\n")
    s <- Sys.time(); cat("Generating data...\n")
    data <- generate_simulation_data(model_name, nrep, len)
    cat("Data generation completed in", format(Sys.time() - s), "\n")
    nfreq <- length(data$mtfreq)

    # k-means clustering of cosine basis expansion coefficients
    cat("Running cluster -> BMARD -> band estimation procedure...\n")
    s <- Sys.time()
    kmlabels <- cluster_cosine_basis(data$mtfreq, data$mtspec, nclust,
                                     nbasis = 15)

    # BMARD for peak detection
    n_samples <- 500; n_chains <- 3
    len <- nrow(data$x); sample_freq <- 1

    ## periodograms
    pgrams <- list()
    for (i in 1:ncol(data$mtspec)) {
      x <- data$x[,i]; x <- (x - mean(x)) / sd(x)
      for (j in 1:n_chains) {
        pgrams[[n_chains * (i - 1) + j]] <- TSA::periodogram(x, plot = FALSE)$spec[-(len / 2)]
      }
    }

    ## sampling

    BDP <- parallel::mclapply(pgrams, function(x) {
      multirun(data$mtfreq, x, n_samples, sample_freq)
    }, mc.cores = ncores)

    ## process BMARD samples for each replicate
    cat("Processing samples...\n")
    # out <- matrix(nrow = ncol(data$mtspec), ncol = 3)
    out <- data.frame()
    for (i in 1:ncol(data$mtspec)) {
      start <- (n_chains * (i - 1) + 1)
      end <- start + n_chains - 1
      BDPout <- BDP[start:end]
      try({
        listmodes <- modes_extraction_gaussianmix(data = BDPout,
                                                  threshold_weight = .01,
                                                  chains = n_chains,
                                                  Nsamp = n_samples,
                                                  freq_rate = len,
                                                  quant = .9,
                                                  percsamp = .6,
                                                  qlow = .025,
                                                  qsup = .975)
      })
      peaks <- c(0, sort(listmodes$globalcenter[1,]), 0.5)
      endpoints <- peaks[1:(length(peaks) - 1)] + diff(peaks) / 2
      endpoints_index <- rep(0, length(endpoints))
      for(b in 1:length(endpoints)) {
        endpoints_index[b] <- which.min(abs(endpoints[b] - data$mtfreq))
      }
      if (endpoints_index[1] == 1) endpoints_index[1] <- 2
      if (endpoints_index[length(endpoints)] == nfreq + 1) {
        endpoints_index[length(endpoints)] <- nfreq
      }

      nfreq <- nrow(data$mtspec)
      clust_spec <- data$mtspec[, kmlabels == kmlabels[i], drop = FALSE]
      collapsed <- avg_collapsed_by_index(clust_spec, endpoints_index)
      expanded <- rep(collapsed, diff(c(1, endpoints_index, nfreq + 1)))
      clust_loss <- sum((clust_spec - expanded)^2) / (2 * (nfreq + 1))

      new.row <- data.frame(
        peaks = I(list(peaks)),
        endpoints = I(list(endpoints)),
        endpoints_index = I(list(endpoints_index)),
        clust_loss = clust_loss,
        label = kmlabels[i]
      )
      out <- rbind(out, new.row)
    }

    ## aggregate by cluster and compute average endpoints and bandwidths
    endpoints <- list()
    clust_loss <- rep(0, 3)
    out$n_endpoints <- sapply(out$endpoints, length)
    for (j in 1:3) {
      outj <- out[out$label == j, ]
      component_count <- table(outj$n_endpoints)
      mode <- which.max(component_count); mode <- as.integer(names(mode))
      outj <- outj[outj$n_endpoints == mode, ]
      min_loss <- outj[which.min(outj$clust_loss),]
      endpoints[j] <- min_loss$endpoints_index
      clust_loss[j] <- min_loss$clust_loss
    }
    bmard_time <- Sys.time() - s
    cat("Finished in ", bmard_time <- format(Sys.time() - s), "\n")

    s <- Sys.time(); cat("Running FBAM on generated data...\n")
    fit <- fbam(data$x, nclust_grid = nclust, nbands_grid = 2:6)
    cat("FBAM finished in", FBAM_time <- format(Sys.time() - s), "\n")

    cat("Saving results to disk...\n")
    sim_data[[n]] <- list(
      data = data,
      kmlabels = kmlabels,
      bmard_opt = endpoints,
      bmard_opt_l = sum(clust_loss),
      FBAM_fit = fit,
      bmard_time = bmard_time,
      FBAM_time = FBAM_time
    )
    save(sim_data, file = DATA_FILE_NAME)
  }, error = function(e) {
    cat("Simulation repetition n = ", n, "failed with error\n\n", str(e))
  })
}
cat("\n\nRun completed at", format(Sys.time(), usetz = TRUE), "\n")
cat("Total runtime:", format(Sys.time() - start, usetz = TRUE), "\n\n")

