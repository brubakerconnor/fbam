#' Cluster/band quantity selection index
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
selection_index <- function(spec, labels, endpoints_index, collapsed) {
  nclust <- dim(endpoints_index)[1]
  band_index <- avg_global_band_redundancy(spec, labels,
                                           endpoints_index, collapsed)
  if (nclust > 1) {
    clust_index <- avg_clust_similarity(spec, labels,
                                        endpoints_index, collapsed)
  } else {
    clust_index <- 0
  }
  return(clust_index + band_index)
}
