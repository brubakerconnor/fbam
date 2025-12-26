## FBAM demo script using simulated data
source("fbam_Rfunctions.R")

## simulate data from one of the models defined in "fbam_Rfunctions.R"
## model1, model2a, model2b, model2c, or model3
model_name <- "model2c"
nrep <- 10 ## number of independent replicates per subpopulation
len <- 500 ## length of each replicate time series
dat <- get(model_name)(nrep, len)

## run with a fixed number of subpopulations and automatically determine 
## the number of frequency bands L = 2, ..., 6
## this function is run serially by default. specify ncores to run in parallel
fbam_out <- fbam(dat$x, nbands = 2:6, nsubpop = 3)

## automatically select number of subpopulations and frequeny bands
## from J = 2, ..., 6 and L = 2, ..., 6
fbam_out <- fbam(dat$x, nbands = 2:6, nsubpop = 2:6)

## if J = 1 is included in the vector of values for nsubpop, 
## two solutions are returned:
## (1) solution1: nsubpop = 1
## (2) solution: solution with automatically selection of J over values > 1
fbam_out <- fbam(dat$x, nbands = 2:6, nsubpop = 1, ncores = 4)

# structure and save output
run <- list(
  fbam_out = fbam_out,
  time = fbam_time,
  model = model,
  nrep = nrep,
  len = len
)
save(run, file = fname)
