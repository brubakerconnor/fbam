#' Distance between two clusters
#'
#' The distance between two clusters is defined as the L2 distance between
#' their collapsed spectra.
#'
#' @param endpoints_index1 Vector of Fourier indices representing the
#' frequency bands associated to the first cluster being compared. It is
#' assumed these include `1` and `nfreq + 1`.
#' @param endpoints_index2 Vector of Fourier indices representing the
#' frequency bands associated to the second cluster being compared. It is
#' assumed these include `1` and `nfreq + 1`.
#' @param collapsed1 Vector of collapsed measures for the first cluster being
#' compared.
#' @param collapsed2 Vector of collapsed measures for the second cluster being
#' compared.
#'
#' @return Double.
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
  sqrt(sum((expanded1 - expanded2)^2) / (2 * (nfreq + 1)))
}

#' Cluster variance
#'
#' Cluster variance (used in the numerator of similarity) is defined as the
#' standard deviation of the L2 distances between each spectra in the cluster
#' and the cluster-average collapsed spectrum.
#'
#' @param clust_spec A matrix of shape (`nfreq`, `nrep_clust`) where
#' `nrep_clust` is the number of spectra clustered into the current cluster.
#' @param endpoints_index Vector of Fourier indices representing the
#' frequency bands associated to this cluster. It is assumed these include `1`
#' and `nfreq + 1`.
#' @param collapsed Vector of collapsed measures for this cluster.
#'
#' @return Double.
clust_sd <- function(clust_spec, endpoints_index, collapsed) {
  nfreq <- dim(clust_spec)[1]; nrep <- dim(clust_spec)[2]

  # ensure sorted endpoint indices
  # it is assumed these endpoints include 1 and nfreq + 1
  endpoints_index <- sort(endpoints_index)

  # expand the collapsed measures to full collapsed spectra
  # and compute the sd of the L2 distance around this object
  expanded <- rep(collapsed, diff(endpoints_index))
  sqrt(sum((clust_spec - expanded)^2) / (2 * (nfreq + 1) * nrep))
}

#' Cluster similarity
#'
#' Cluster similarity is defined as the ratio of the sum of standard deviations
#' of the distances of each cluster spectra to their respective cluster centers
#' (collapsed spectra) to the distance between the collapsed spectra of each
#' cluster.
#'
#' @param clust_spec1 A matrix of shape (`nfreq`, `nrep_clust1`) where
#' `nrep_clust1` is the number of spectra clustered into the current cluster.
#' @param clust_spec2 A matrix of shape (`nfreq`, `nrep_clust2`) where
#' `nrep_clust2` is the number of spectra clustered into the current cluster.
#' @param endpoints_index1 Vector of Fourier indices representing the
#' frequency bands associated to the first cluster being compared. It is
#' assumed these include `1` and `nfreq + 1`.
#' @param endpoints_index2 Vector of Fourier indices representing the
#' frequency bands associated to the second cluster being compared. It is
#' assumed these include `1` and `nfreq + 1`.
#' @param collapsed1 Vector of collapsed measures for the first cluster being
#' compared.
#' @param collapsed2 Vector of collapsed measures for the second cluster being
#' compared.
#'
#' @return Double.
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

#' Cluster similarity index
#'
#' The index of cluster similarity is the average similarity of each cluster
#' with its most similar cluster.
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


