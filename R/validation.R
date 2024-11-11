validation_index <- function(spec, labels, endpoints_index, avg_summary) {
  nsubpop <- dim(endpoints_index)[1]
  band_index <- avg_global_band_sim(spec, labels, endpoints_index, avg_summary)
  subpop_index <- 0
  if (nsubpop > 1) subpop_index <- avg_subpop_sim(spec, labels, endpoints_index, avg_summary)
  return(subpop_index + band_index)
}

subpop_dist <- function(endpoints_index1, endpoints_index2, avg_summary1, avg_summary2) {
  # ensure sorted endpoint indices
  # it is assumed these endpoints_index include 1 and nfreq + 1
  endpoints_index1 <- sort(endpoints_index1)
  endpoints_index2 <- sort(endpoints_index2)

  # expand the collapsed measures to full collapsed spectra
  # and compute L2 distance
  expanded1 <- rep(avg_summary1, diff(endpoints_index1))
  expanded2 <- rep(avg_summary2, diff(endpoints_index2))
  nfreq <- length(expanded1)
  sqrt(sum((expanded1 - expanded2)^2))
}

subpop_sd <- function(subpop_spec, endpoints_index, avg_summary) {
  nrep <- dim(subpop_spec)[2]
  endpoints_index <- sort(endpoints_index)
  expanded <- rep(avg_summary, diff(endpoints_index))
  sqrt(sum((subpop_spec - expanded)^2 / nrep))
}

subpop_sim <- function(subpop_spec1, subpop_spec2, endpoints_index1,
                       endpoints_index2, avg_summary1, avg_summary2) {
  sd1 <- subpop_sd(subpop_spec1, endpoints_index1, avg_summary1)
  sd2 <- subpop_sd(subpop_spec2, endpoints_index2, avg_summary2)
  d <- subpop_dist(endpoints_index1, endpoints_index2, avg_summary1, avg_summary2)
  return((sd1 + sd2) / d)
}

avg_subpop_sim <- function(spec, labels, endpoints_index, avg_summary) {
  # allows for the possibility of empty subpopulations
  ulabels <- unique(labels)
  n_nonempty_subpop <- length(ulabels)
  nfreq <- dim(spec)[1]; nsubpop <- dim(endpoints_index)[1]

  # compute the similarity between each pair of subpopulations
  # the resulting matrix is symmetric so we only need the lower triangle
  sim_mat <- matrix(0, nrow = nsubpop, ncol = nsubpop)
  for (i in ulabels) {
    subpop_speci <- spec[, labels == i, drop = F]
    endpoints_indexi <- c(1, endpoints_index[i,], nfreq + 1)
    avg_summaryi <- avg_summary[i,]

    for (j in ulabels) {
      if (j < i) {
        subpop_specj <- spec[, labels == j, drop = F]
        endpoints_indexj <- c(1, endpoints_index[j,], nfreq + 1)
        avg_summaryj <- avg_summary[j,]
        sim_mat[i, j] <- subpop_sim(
          subpop_speci, subpop_specj,
          endpoints_indexi, endpoints_indexj,
          avg_summaryi, avg_summaryj
        )
      }
    }
  }
  sim_mat <- sim_mat + t(sim_mat) # make symmetric
  return(sum(apply(sim_mat, 1, max)) / n_nonempty_subpop)
}

band_sim <- function(band_spec1, band_spec2, avg_summary1, avg_summary2) {
  sd1 <- sqrt(sum((band_spec1 - avg_summary1)^2))
  sd2 <- sqrt(sum((band_spec2 - avg_summary2)^2))
  union_band_spec <- rbind(band_spec1, band_spec2)
  union_collapsed <- mean(union_band_spec) # unioned collapsed measure of power
  d <- sqrt(sum((union_band_spec - union_collapsed)^2))
  return((sd1 + sd2) / d)
}

avg_subpop_band_sim <- function(subpop_spec, endpoints_index, avg_summary) {
  nbands <- length(avg_summary); nfreq <- dim(subpop_spec)[1]

  # compute similarity of each band with its right neighboring band
  sim_vec <- rep(0, nbands - 1)
  for (b in 1:(nbands - 1)) {
    # endpoints, band spectra and collapsed measures for left band
    left_endpoint1 <- endpoints_index[b]
    right_endpoint1 <- endpoints_index[b + 1]
    band_spec1 <- subpop_spec[left_endpoint1:(right_endpoint1 - 1), , drop = F]
    avg_summary1 <- avg_summary[b]

    # endpoints, band spectra and collapsed measures for right band
    left_endpoint2 <- endpoints_index[b + 1]
    right_endpoint2 <- endpoints_index[b + 2]
    band_spec2 <- subpop_spec[left_endpoint2:(right_endpoint2 - 1), , drop = F]
    avg_summary2 <- avg_summary[b + 1]

    # similarity between left and right bands
    sim_vec[b] <- band_sim(band_spec1, band_spec2, avg_summary1, avg_summary2)
  }
  # return(sum(sim_vec) / (nbands - 1))
  return(sum(sim_vec) / nbands)
}

avg_global_band_sim <- function(spec, labels, endpoints_index, collapsed) {
  ulabels <- unique(labels)
  n_nonnempty_subpop <- length(ulabels)
  nsubpop <- dim(endpoints_index)[1]; nfreq <- dim(spec)[1]
  sim_vec <- rep(0, nsubpop)
  for (i in ulabels) {
    subpop_speci <- spec[, labels == i, drop = F]
    endpoints_indexi <- c(1, endpoints_index[i,], nfreq + 1)
    avg_summaryi <- collapsed[i,]
    sim_vec[i] <- avg_subpop_band_sim(subpop_speci, endpoints_indexi, avg_summaryi)
  }
  return(mean(sim_vec))
}
