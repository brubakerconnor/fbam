#' Collapsed measures from frequency boundaries
#'
#' Computes replicate-specific or cluster-average collapsed measures of power
#' (band average power) using frequency-valued band boundaries for possible
#' more than one cluster.
#'
#' @param spec Vector (one replicate) or matrix of shape (`nfreq`, `nrep`) where
#' `nfreq` is the number of Fourier frequencies used in estimation
#' (typically `floor(T / 2) - 1`) and `nrep` is the number of replicates.
#' @param mtfreq Vector of frequencies corresponding to each row of `spec`.
#' @param endpoints Vector of length `nbands - 1` or `nbands + 1`
#' (in the case of one cluster) or matrix (in the case of more than one cluster)
#' of shape (`nclust`, `nbands - 1`) or (`nclust`, `nbands + 1`).
#' See below for more details.
#' @param labels If computing collapsed measures for more than one cluster,
#' `labels` should be a vector of length `nrep` where the `k`th value points
#' to the row of `endpoints_index` that should be used to collapse the `k`th
#' spectrum (column) in `spec`. Default `NULL`. If labels are not provided but
#' endpoints_index includes more than one row, the first row will be used in
#' collapsing all replicates.
#'
#' @details
#' The vector or matrix `endpoints_index` can be either of length `nbands - 1`
#' or shape (`nclust`, `nbands - 1`) or of length `nbands + 1` or shape
#' (`nclust`, `nbands + 1`). The former case is appropriate when the first and
#' last frequencies are NOT included, and the latter case for when they are. The
#' first frequency is always `0` and the last frequency is always `0.5`. These
#' are optional, but if they are included, they must be included for ALL clusters,
#' that is, in each row of `endpoints`.
#'
#' @return A matrix of shape (`nrep`, `nbands`) for replicate-specific collapsed
#' measures or a matrix of shape (`nclust`, `nbands`) for cluster-average
#' collapsed measures.
#' @export
#'
#' @examples
#' # 10 replicates of white noise where T = 400
#' # spec is simulated using the asymptotic distribution of the periodogram
#' # for a white noise process with variance 1
#' length <- 400; nfreq <- floor(length / 2) - 1
#' nrep <- 10
#' mtfreq <- seq(1 / length, by = 1 / length, length = nfreq)
#' spec <- matrix(rchisq(nfreq * nrep, df = 2) / 2, nrow = nfreq)
#' labels <- rep(1:2, each = nrep / 2) # equally divided clusters
#'
#' # first cluster has frequency band endpoints at 50 / T and 100 / T
#' # second cluster has frequency band endpoints at 70 / T and 120 / T
#' endpoints <- matrix(c(0.125, 0.25, 0.175, 0.3), nrow = 2, byrow = T)
#' rep_collapsed(spec, mtfreq, endpoints, labels)
#' avg_collapsed(spec, mtfreq, endpoints, labels)
rep_collapsed <- function(spec, mtfreq, endpoints, labels = NULL) {
  if (is.vector(endpoints)) endpoints <- matrix(endpoints, nrow = 1)
  # convert endpoints to indices and send off to
  # rep_collapsed_by_index to handle the rest
  endpoints_index <- apply(endpoints, c(1, 2), g, ref_vector = mtfreq)
  return(rep_collapsed_by_index(spec, endpoints_index, labels))
}

#' @rdname rep_collapsed_by_index
#' @export
avg_collapsed <- function(spec, mtfreq, endpoints, labels = NULL) {
  if (is.vector(endpoints)) endpoints <- matrix(endpoints, nrow = 1)
  # convert endpoints to indices and send off to
  # rep_collapsed_by_index to handle the rest
  endpoints_index <- apply(endpoints, c(1, 2), g, ref_vector = mtfreq)
  return(avg_collapsed_by_index(spec, endpoints_index, labels))
}

# determine the index of the value in ref_vector that is closest to value in
# absolute distance
g <- function(value, ref_vector) {
  which.min(abs(ref_vector - value))
}
