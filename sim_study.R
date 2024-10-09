# study parameters from command line
# arguments 7-10 are optional. if not specified, default parameter values will
# be used in running the genetic algorithm.
# [1] model_name        string    e.g., "model1", "model2a"
# [2] nrep              integer   number of replicate time series (per cluster)
# [3] len               integer   length of the time series epochs
# [4] nsim              integer   number of simulation repetitions (e.g., 100)
# [5] results_dir       string    where to save the simulation results (.rda file)
# [6] ncores            integer   number of cluster cores to use
# [7] pcrossover        float     crossover probability (optional)
# [8] pmutate           float     mutation probability (optional)
# [9] endpoint_range    integer   endpoint mutation range (optional)
# [10] collapsed_range  float     collapsed measure mutation range (optional)
args <- commandArgs(trailingOnly = T)
model_name <- args[1]
nrep <- as.integer(args[2])
len <- as.integer(args[3])
nsim <- as.integer(args[4])
results_dir <- args[5]
ncores <- as.integer(args[6])
pcrossover <- as.double(args[7])
pmutate <- as.double(args[8])
endpoint_range <- as.integer(args[9])
collapsed_range <- as.double(args[10])

# make sure results directory exists
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# set seed for reproducibility
set.seed(451)
source("R/fbam.R") # functions for optimization
source("R/generate_sim_data.R") # functions for data generation

# print parameter settings to log
cat("Running simulation study with the following parameters:\n",
    "\tmodel_name\t\t=", model_name, "\n",
    "\tnrep\t\t\t=", nrep, "\n",
    "\tlen\t\t\t=", len, "\n",
    "\tnsim\t\t\t=", nsim, "\n",
    "\tRESULTS_DIR\t\t=", results_dir, "\n",
    "\tncores\t\t\t=", ncores, "\n",
    "\tpcrossover\t\t=", pcrossover, "\n",
    "\tpmutate\t\t\t=", pmutate, "\n",
    "\tendpoint_range\t\t=", endpoint_range, "\n",
    "\tcollapsed_range\t\t=", collapsed_range, "\n")

# path to results data file
# remove trailing separator if needed; this only works on UNIX machines
if (substr(results_dir, nchar(results_dir), nchar(results_dir)) == "/") {
  results_dir <- substr(results_dir, 1, nchar(results_dir) - 1)
}
data_fname <- file.path(results_dir,
                        paste0(model_name,
                               "_nrep=", nrep,
                               "_len=", len,
                               "_pcross=", pcrossover,
                               "_pmutate=", pmutate,
                               "_endpoint=", endpoint_range,
                               "_collapsed=", collapsed_range,
                               ".rda"))
cat("Results will be saved at", data_fname, "\n")

# create cluster for parallelization
nattempt <- 1
while (is.numeric(nattempt)) {
  nattempt <- tryCatch({
    cl <- parallelly::makeClusterPSOCK(ncores)
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function() {
      source('R/fbam.R'); source('R/generate_sim_data.R')
    })  # make needed functions available to each worker
    cat("Cluster creation succesful.\n")
  }, error = function(e) {
    message("attempt: ", nattempt, " ", conditionMessage(e))
    if (nattempt > 10) stop("too many failures", call.=FALSE)
    nattempt + 1
  })
}

# default parameter settings for optional command line parameters
if (is.na(pcrossover)) pcrossover <- 0.5
if (is.na(pmutate)) pmutate <- 0.05
if (is.na(endpoint_range)) endpoint_range <- NULL
if (is.na(collapsed_range)) collapsed_range <- NULL

# run simulation nsim times and save results to disk at each iteration
start <- Sys.time()
sim_data <- list()
for (n in 1:nsim) {
  tryCatch({
    cat("Simulation repetition n =", n, "of", nsim, "starting at",
        format(Sys.time(), usetz = TRUE), "\n\n")
    s <- Sys.time(); cat("Generating data...\n")
    data <- generate_simulation_data(model_name, nrep, len)
    cat("Data generation completed in", format(Sys.time() - s), "\n")

    s <- Sys.time(); cat("Running FBAM on generated data...\n")
    fit <- fbam(X = data$x, nclust_grid = 2:4, nbands_grid = 2:4,
                pcrossover = pcrossover, pmutate = pmutate,
                endpoint_range = endpoint_range,
                collapsed_range = collapsed_range)
    cat("FBAM finished in", runtime <- format(Sys.time() - s), "\n")

    cat("Saving results to disk...\n")
    sim_data[[n]] <- list(data = data, fit = fit, time = runtime)
    save(sim_data, file = data_fname)
  }, error = function(e) {
    cat("Simulation repetition n = ", n, "failed with error\n\n", str(e))
  })
}
cat("\n\nRun completed at", format(Sys.time(), usetz = TRUE), "\n")
cat("Total runtime:", format(Sys.time() - start, usetz = TRUE), "\n\n")
