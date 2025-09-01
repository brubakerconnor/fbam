source("fbam_Rfunctions.R")

# simulation parameters
# [1] model             string    e.g., "model1", "model2a"
# [2] nrep              integer   number of replicate time series (per subpop)
# [3] len               integer   length of the time series epochs
# [4] nsim              integer   number of simulation repetitions (e.g., 100)
# [5] ncores            integer   number of cores to use for parallelization
# [6] results           string    directory to save output (.rds file)
# [7] run_number        integer   index in sequence of runs
# [8] ntapers           double    (optional) exponential factor of length for
#                                 number of tapers (e.g, 0.75 = T^0.75).
#                                 defaults to 0.5
args <- commandArgs(trailingOnly = T)
model <- args[1]
nrep <- as.integer(args[2])
len <- as.integer(args[3])
nsim <- as.integer(args[4])
ncores <- as.integer(args[5])
results <- args[6]
run_number <- args[7]
taper_exp <- as.numeric(args[8])
if (is.na(taper_exp)) {
  ntapers <- floor(len ^ 0.5)
} else {
  ntapers <- floor(len ^ taper_exp)
}
cat("STUDY PARAMETERS:\n", "nsim: ", nsim, "\n", "ncores: ", ncores, "\n",
    "model_name: ", model, "\n", "nrep: ", nrep, "\n", "len: ", len, "\n",
    "results: ", results, "\n", "ntapers: ", ntapers, "\n")

# set up run
fname <- file.path(results, paste0(model, "_nrep=", nrep, "_len=",
                                   len, "_run=", run_number, ".rds"))
cat("results saved to: ", fname, "\n")
dat <- get(model)(nrep, len)

# fbam
run_start <- Sys.time()
fbam_out <- fbam(dat$x, nbands = 2:6, nsubpop = 2:6,
                 ncores = ncores, ntapers = ntapers, verbose = FALSE)
fbam_time <- as.numeric(Sys.time() - run_start, units = "secs")

# structure and save output
run <- list(fbam_out = fbam_out, time = fbam_time,
            model = model, nrep = nrep, len = len)
save(run, file = fname)
