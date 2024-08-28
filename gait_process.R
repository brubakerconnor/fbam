rm(list = ls())
source('R/fbam.R')
set.seed(9)

#### download and process gait data from physionet #############################
# download files from physionet directly if not already
print('Downloading data files from physionet...')
if (!dir.exists('data/')) {
  dir.create('data/gait-processed/')
  system('wget -r -N -c -np https://physionet.org/files/gaitndd/1.0.0/ -P data/')
}
DATA_DIR <- 'data/physionet.org/files/gaitndd/1.0.0'

# get the path to each of the ts files in the collection
ts_files <- list.files(DATA_DIR, pattern = "\\.ts$", full.names = T)
patient_id <- tools::file_path_sans_ext(basename(ts_files))
labels <- gsub("[0-9]", "", patient_id)

# note the following patients were removed due to too many artefacts in the
# stride interval series: als12, als4, park11
skip <- c("als12", "als4", "park11") # skip these
ts_files <- ts_files[!patient_id %in% skip]
labels <- labels[!patient_id %in% skip]
patient_id <- patient_id[!patient_id %in% skip]

# signal and estimation parameters
duration <- 210 # seconds
sample_rate <- 2 # samples per second
ts_length <- duration * sample_rate
mtfreq <- seq(from = 1 / ts_length, by = 1 / ts_length,
              length = floor(ts_length / 2) - 1)
ntapers <- floor(sqrt(ts_length))

# process each signal
xraw <- list()
x <- matrix(nrow = ts_length, ncol = length(ts_files))
mtspec <- matrix(nrow = length(mtfreq), ncol = length(ts_files))
for (i in 1:length(ts_files)) {
  print(paste0("Processing data contained in ", ts_files[i]))

  # read data and define variables
  dat <- readr::read_delim(ts_files[i], col_names = F, show_col_types = F)
  # left_stride <- dat[[2]]
  elapsed_time <- dat[[1]]
  left_stride <- diff(elapsed_time)
  elapsed_time <- elapsed_time[-1]
  xraw[[i]] <- list()
  xraw[[i]]$start <- dat[[1]][1]
  xraw[[i]]$elapsed <- dat[[1]]
  xraw[[i]]$left <- left_stride

  # outlier filtering
  # time points that lie outside the percentiles defined below are considered
  # outliers and are replaced by the median value of the series
  median_stride <- median(left_stride)
  quantiles <- quantile(left_stride, c(0.01, 0.96))
  replace <- left_stride < quantiles[1] | left_stride > quantiles[2]
  xraw[[i]]$flagged <- replace
  left_stride[replace] <- NA
  left_stride <- imputeTS::na_ma(left_stride, k = 6, weighting = 'exponential')

  # linear interpolation to make regularly sampled series
  start_time <- min(elapsed_time); end_time <- start_time + duration
  time_grid <- seq(from = start_time, to = end_time, by = 1 / sample_rate)
  time_grid <- time_grid[1:ts_length]
  left_stride <- approx(x = elapsed_time, y = left_stride, xout = time_grid,
                        method = "linear")$y

  # remove trend
  left_stride <- as.vector(gsignal::detrend(left_stride, p = 1))

  # standardize by standard deviation of series
  left_stride <- left_stride / sd(left_stride)
  x[, i] <- left_stride # add to ts data

  # multitaper estimate
  mtspec[, i] <- as.vector(sine_multitaper(left_stride)$mtspec)
}

#### descriptors file ##########################################################
# read and clean up the subject-level descriptors file
print('Processing subject descriptors file...')
desc <- read.delim("data/physionet.org/files/gaitndd/1.0.0/subject-description.txt",
                   row.names = NULL, na.strings = "MISSING")
names(desc) <- c("patient_id", "group", "age", "height", "weight", "gender",
                 "gait_speed", "duration_severity")
desc <- desc[desc$patient_id != "", ] # remove that weird mostly empty row
desc$group[desc$group == "subjects"] <- "als"

# impute those missing als gait speed values with group mean gait speed
als_mean_speed <- mean(desc[desc$group == "als", 'gait_speed'], na.rm = TRUE)
desc$gait_speed[is.na(desc$gait_speed) & desc$group == 'als'] <- als_mean_speed

# The duration/severity column contains condition-specific measures of disease
# duration or severity. For control patients, this value is 0. For Parkinson's
# patients, this is the Hohn and Yahr score which ranges from 1 to 5 with higher
# values indicating more advanced disease. For Huntington's patients, this is
# the total functional capacity measure (lower scores mean more advanced
# functional impairment) which ranges from 0 to 13. Finally, for ALS patients,
# this is simply the number of months since diagnosis.

# in order to provide a standardized scale, the duration or severity within
# each group is re-calculated in the following way:
# PD: 0.25 * score - 0.25
# HD: (1/13) * (13 - score) - 1/13
# ALS: (1/max_duration) * score - (1/max_duration)
# Standardization is done in this way to put all values between 0 and 1 where
# a value of 0 is "best" and 1 is "worst".
desc$duration_severity_std <- desc$duration_severity
desc$duration_severity_std[desc$group == "park"] <- 0.25 *
  desc$duration_severity[desc$group == "park"] - 0.25
desc$duration_severity_std[desc$group == "hunt"] <- (1/13) *
  (13 - desc$duration_severity[desc$group == "hunt"] - (1/13))

max_duration <- max(desc$duration_severity_std[desc$group == "als"])
desc$duration_severity_std[desc$group == "als"] <- (1 / max_duration) *
  desc$duration_severity[desc$group == "als"] - (1 / max_duration)

#### self similarity parameter #################################################
print('Computing self similarity parameters...')
ssp <- rep(0, ncol(x))
for (i in 1:ncol(x)) {
  print(patient_id[i])
  ssp[i] <- DFA::SSP(x[,i])
}

#### FBAM no clustering and clustering solutions ###############################
print('Running FBAM...')
## setting (a)
noclust <- fbam(x, 1, 2:6, parallel = TRUE)$selected
names(noclust$rep_collapsed) <- c("LF", "HF")

## setting (b)
clust <- fbam(x, 2:6, 2:6, parallel = TRUE)$selected
names(clust$rep_collapsed) <- c("LF", "HF")

x <- data.frame(x); names(x) <- patient_id
mtspec <- data.frame(mtspec); names(mtspec) <- patient_id
gait <- list(
  xraw = xraw,
  x = x,
  mtfreq = mtfreq,
  mtspec = mtspec,
  labels = labels,
  patient_id = patient_id,
  noclust = noclust,
  clust = clust,
  ssp = ssp
)

#### csv file ##################################################################
# create simplified csv file
# add the LF, HF, and LF:HF measures from the no clustering and clustering
# solutions along with cluster labels from the clustering solution
noclust_df <- data.frame(patient_id = gait$patient_id,
                         noclust_LF = noclust$rep_collapsed[,1],
                         noclust_HF = noclust$rep_collapsed[,2],
                         noclust_ratio = noclust$rep_collapsed[,1] / noclust$rep_collapsed[,2])

clust_df <- data.frame(patient_id = gait$patient_id,
                       clust_LF = clust$rep_collapsed[,1],
                       clust_HF = clust$rep_collapsed[,2],
                       clust_ratio = clust$rep_collapsed[,1] / clust$rep_collapsed[,2],
                       clust_label = gait$clust$labels)

gait_df <- merge(desc, noclust_df, by = "patient_id")
gait_df <- merge(gait_df, clust_df, by = "patient_id")

# add endpoint columns and ssp
gait_df$noclust_endpoint <- noclust$endpoints[1]
gait_df$clust_endpoint <- clust$endpoints[gait_df$clust_label,]
gait_df$ssp <- gait$ssp

#### save data #################################################################
# create data object and save
print('Saving data...')
save(gait, file = "data/gait-processed/gait.rda")

# save results to csv
write.csv(gait_df, file = "data/gait-processed/gait.csv", row.names = FALSE)

