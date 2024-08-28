args <- commandArgs(trailingOnly = T)
DATA_DIR <- args[1]  # directory where .rda files from comparison study are saved

# utility functions
# convert run time to seconds if not already
convert_to_seconds <- function(s) {
  unit <- sub("^[0-9.]+ ", "", s)
  if (unit == "secs") {
    time <- sub(" secs$", "", s)
    return(as.double(time) * 60)
  }
  if (unit == "mins") {
    time <- sub(" mins$", "", s)
    return(as.double(time) * 60)
  }
}

# function to get mode of a vector
get_mode <- function(v) {
  uniq_vals <- unique(v)
  uniq_vals[which.max(tabulate(match(v, uniq_vals)))]
}

# collect all file names to include in the report
fnames <- list.files(DATA_DIR, pattern = "\\.rda$", full.names = TRUE)
report <- data.frame()

for (fname in fnames) {
  print(fname)
  load(fname)

  # get parameters of the run
  model <- strsplit(tail(strsplit(fname, "/")[[1]], 1), "_")[[1]][1]
  nrep <- dim(sim_data[[1]]$data$x)[2] / length(unique(sim_data[[1]]$data$labels))
  len <- dim(sim_data[[1]]$data$x)[1]
  nsim <- length(sim_data)
  pcrossover <- sim_data[[1]]$FBAM_fit$selected$params$pcrossover
  pmutate <- sim_data[[1]]$FBAM_fit$selected$params$pmutate
  endpoint_range <- sim_data[[1]]$FBAM_fit$selected$params$endpoint_range[1]
  collapsed_range <- sim_data[[1]]$FBAM_fit$selected$params$collapsed_range[1]

  # collect data
  FBAM_ari <- BMARD_ari <- rep(0, nsim)
  FBAM_nbands <- BMARD_nbands <- rep(0, nsim)
  FBAM_obj <- BMARD_obj <- FBAM_time <- BMARD_time <- rep(0, nsim)
  for (n in 1:nsim) {
    FBAM_ari[n] <- mclust::adjustedRandIndex(
      sim_data[[n]]$FBAM_fit$selected$labels,
      sim_data[[n]]$data$labels
    )
    BMARD_ari[n] <- mclust::adjustedRandIndex(
      sim_data[[n]]$kmlabels,
      sim_data[[n]]$data$labels
    )
    FBAM_nbands[n] <- sim_data[[n]]$FBAM_fit$selected$params$nbands
    BMARD_nbands[n] <- mean(unlist(lapply(sim_data[[n]]$bmard_opt, length)))
    FBAM_obj[n] <- sim_data[[n]]$FBAM_fit$selected$objective
    BMARD_obj[n] <- sim_data[[n]]$bmard_opt_l
    FBAM_time[n] <- convert_to_seconds(sim_data[[n]]$FBAM_time)
    BMARD_time[n] <- convert_to_seconds(sim_data[[n]]$bmard_time)
  }
  x <- data.frame(
    model = model,
    nrep = nrep,
    len = len,
    pcrossover = pcrossover,
    pmutate = pmutate,
    endpoint_range = endpoint_range,
    collapsed_range = collapsed_range,
    mean_FBAM_ari = mean(FBAM_ari, na.rm = T),
    se_FBAM_ari = sd(FBAM_ari, na.rm = T),
    mean_BMARD_ari = mean(BMARD_ari, na.rm = T),
    se_BMARD_ari = sd(BMARD_ari, na.rm = T),
    mode_FBAM_nbands = get_mode(FBAM_nbands),
    percent_FBAM_nbands = sum(FBAM_nbands == get_mode(FBAM_nbands)),
    mode_BMARD_nbands = get_mode(BMARD_nbands),
    percent_BMARD_nbands = sum(BMARD_nbands == get_mode(BMARD_nbands)),
    mean_FBAM_obj = mean(FBAM_obj),
    median_FBAM_obj = median(FBAM_obj),
    se_FBAM_obj = sd(FBAM_obj),
    mad_FBAM_obj = mad(FBAM_obj),
    mean_BMARD_obj = mean(BMARD_obj),
    median_BMARD_obj = median(BMARD_obj),
    se_BMARD_obj = sd(BMARD_obj),
    mad_BMARD_obj = mad(BMARD_obj),
    mean_FBAM_time_secs = mean(FBAM_time),
    se_FBAM_time_secs = sd(FBAM_time),
    mean_BMARD_time_secs = mean(BMARD_time),
    se_BMARD_time_secs = sd(BMARD_time)
  )
  report <- rbind(report, x)
  rm(sim_data)
}
write.csv(report, paste0(DATA_DIR, "/report.csv"), row.names = FALSE)
