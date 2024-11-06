#' @useDynLib fbam
#' @importFrom Rcpp sourceCpp
evaluate <- function(pop, spec, nsubpop, nbands) {
  apply(pop, MARGIN = 1, FUN = function(x) {
    p <- matrix(x, nrow = nsubpop, byrow = TRUE) # organize by subpopulation

    # assign replicates to subpopulations using minimum l2 distance to
    # subpopulation-specific average collapsed spectra
    nsubpop <- dim(p)[1]; nbands <- (dim(p)[2] + 1) / 2
    nfreq <- dim(spec)[1]; nrep <- dim(spec)[2]
    endpoints <- p[, 1:(nbands - 1), drop = FALSE]
    collapsed <- p[, nbands:dim(p)[2], drop = FALSE]
    if (nsubpop > 1) {
      labels <- l2assign(spec, p)
    } else {
      labels <- rep(1, nrep)
    }

    # if there are empty subpopulations, assign fitness penalty
    if (length(unique(labels)) != nsubpop) return(1e-50)

    # otherwise, compute within cluster sum of squares
    wcss <- vector("double", nsubpop)
    for (j in 1:nsubpop) {
      clust_spec <- spec[, labels == j, drop = FALSE]
      # get collapsed measures and widths in each band
      # collapsed <- avg_summary(clust_spec, endpoints[j,])
      widths <- diff(c(1, endpoints[j,], nfreq + 1))
      # expand the collapsed measures and compute L2 distance
      # center <- rep(collapsed, widths)
      center <- rep(collapsed[j,], widths)
      wcss[j] <- sum((clust_spec - center)^2) / (2 * (nfreq + 1))
    }
    return(1 / sum(wcss))
  })
}

# l2_assign <- function(spec, p) {
#   nsubpop <- dim(p)[1]; nbands <- (dim(p)[2] + 1) / 2
#   nfreq <- dim(spec)[1]; nrep <- dim(spec)[2]
#
#   # separate out endpoints and collapsed measures
#   endpoints <- p[, 1:(nbands - 1), drop = FALSE]
#   collapsed <- p[, nbands:dim(p)[2], drop = FALSE]
#
#   # compute distance matrix
#   dist <- matrix(nrow = nrep, ncol = nsubpop)
#   for (j in seq_len(nsubpop)) {
#     # get widths of each band and expand collapsed measures
#     widths <- diff(c(1, endpoints[j, ], nfreq + 1))
#     center <- rep(collapsed[j, ], widths)
#
#     # compute l2 distance between each spectrum and the jth cluster
#     dist[,j] <- colSums((spec - center)^2) / (2 * (nfreq + 1))
#   }
#   return(apply(dist, 1, which.min))
# }
