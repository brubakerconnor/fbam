#' Unioned band standard deviation
#'
#' Used in the denominator of band redundancy, unioned band SD is the
#' standard deviation of the distances of each spectra constrained to the
#' union of two neighboring frequency bands to the collapsed measure of
#' power computed over those two bands.
#'
#' @param union_band_spec A matrix of shape `(union_nfreq, nrep_clust)` where
#' `union_nfreq` is the number of frequencies across both bands being compared
#' and `nrep_clust` is the number of replicates in the cluster.
#' @param nfreq Integer. The number of Fourier frequencies used in estimation.
#' This is to provide proper scaling for integral approximation.
#'
#' @return Double.
unioned_band_sd <- function(union_band_spec, nfreq) {
  union_collapsed <- mean(union_band_spec) # unioned collapsed measure of power
  sqrt(sum((union_band_spec - union_collapsed)^2) / (2 * (nfreq + 1)))
}

#' (Single) Band standard deviation
#'
#' Used in the numerator of band redundancy, band standard deviation is defined
#' as the standard deviation of the distances of each spectra constrained to a
#' frequency band to the collapsed measure of power computed within that band.
#'
#' @param band_spec A matrix of shape `(nfreq, nrep_clust)` where
#' `nfreq` is the number of frequencies in the current band and `nrep_clust`
#' is the number of replicates in the cluster.
#' @param band_collapsed Double. Collapsed measure of power over the current
#' band (as estimated by the GA, for example).
#' @param nfreq Integer. The number of Fourier frequencies used in estimation.
#' This is to provide proper scaling for integral approximation.
#'
#' @return Double.
band_sd <- function(band_spec, band_collapsed, nfreq) {
  sqrt(sum((band_spec - band_collapsed)^2) / (2 * (nfreq + 1)))
}

#' Band redundancy
#'
#' Band redundancy is the ratio of within-band standard deviation of distances
#' to the same standard deviation computed if the bands were union. Higher
#' redundancy means that these two quantities are similar in terms of loss
#' reduction, suggesting that these two bands could be union without significant
#' loss in approximation quality.
#'
#' @param band_spec1 A matrix of shape `(nfreq1, nrep_clust)` where
#' `nfreq1` is the number of frequencies in the first band and `nrep_clust`
#' is the number of replicates in the cluster.
#' @param band_spec2 A matrix of shape `(nfreq2, nrep_clust)` where
#' `nfreq2` is the number of frequencies in the second band and `nrep_clust`
#' is the number of replicates in the cluster.
#' @param collapsed1 Double. Collapsed measure of power over the first
#' band (as estimated by the GA, for example).
#' @param collapsed2 Double. Collapsed measure of power over the second
#' band (as estimated by the GA, for example).
#' @param nfreq Integer. The number of Fourier frequencies used in estimation.
#' This is to provide proper scaling for integral approximation.
#'
#' @return Double.
band_redundancy <- function(band_spec1, band_spec2, collapsed1, collapsed2,
                           nfreq) {
  # compute single band SDs
  sd1 <- band_sd(band_spec1, collapsed1, nfreq)
  sd2 <- band_sd(band_spec2, collapsed2, nfreq)

  # unioned band SD
  union_band_spec <- rbind(band_spec1, band_spec2)
  d <- unioned_band_sd(union_band_spec, nfreq)
  return((sd1 + sd2) / d)
}

#' Average within-cluster band redundancy
#'
#' Average the redundancy of each band with its right neighbor over all
#' bands associated to a single cluster.
#'
#' @param clust_spec A matrix of shape (`nfreq`, `nrep_clust`) where
#' `nrep_clust` is the number of spectra clustered into the current cluster.
#' @param endpoints_index Vector of Fourier indices representing the
#' frequency bands associated to the current cluster. It is assumed these
#' include `1` and `nfreq + 1`.
#' @param collapsed Vector of collapsed measures for this cluster.
#'
#' @return Double.
avg_clust_band_redundancy <- function(clust_spec, endpoints_index, collapsed) {
  nbands <- length(collapsed); nfreq <- dim(clust_spec)[1]

  # compute redundancy of each band with its right neighboring band
  red_vec <- rep(0, nbands - 1)
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

    # redundancy between left and right bands
    red_vec[b] <- band_redundancy(band_spec1, band_spec2,
                                  collapsed1, collapsed2,
                                  nfreq)
  }
  # compute average redundancy over nbands
  return(sum(red_vec) / nbands)
}

#' Average band redundancy across all clusters
#'
#' @param spec A matrix of shape `(nfreq, nrep)` containing all spectra across
#' all clusters being comapred.
#' @param labels Vector of labels.
#' @param endpoints_index A matrix of shape `(nclust, nbands - 1)` with the
#' Fourier indices of the frequency band boundaries for each cluster. It is
#' assumed these endpoints do NOT contain `1` and `nfreq + 1` since these are
#' intended to come directly from the output of the GA which does not carry
#' these values in its representation.
#' @param collapsed A matrix of shape `(nclust, nbands)` with the collapsed
#' measures of each cluster.
#'
#' @return Double.
avg_global_band_redundancy <- function(spec, labels, endpoints_index, collapsed) {
  ulabels <- unique(labels)
  n_nonnempty_clust <- length(ulabels)
  nclust <- dim(endpoints_index)[1]; nfreq <- dim(spec)[1]

  red_vec <- rep(0, nclust)
  for (i in ulabels) {
    # set aside spectra, endpoints (which are then padded), and collapsed
    # measures for cluster i
    clust_spec <- spec[, labels == i, drop = F]
    endpoints_index_i <- c(1, endpoints_index[i, ], nfreq + 1)
    collapsed_i <- collapsed[i, ]

    # average band redundancy for cluster i
    red_vec[i] <- avg_clust_band_redundancy(clust_spec,
                                            endpoints_index_i,
                                            collapsed_i)
  }
  return(mean(red_vec))
}
