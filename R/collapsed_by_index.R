#' Collapsed measures from Fourier indices
#'
#' Computes replicate-specific or cluster-average collapsed measures from
#' Fourier indices for possibly more than one cluster.
#'
#' @param spec Vector (one replicate) or matrix of shape (`nfreq`, `nrep`) where
#' `nfreq` is the number of Fourier frequencies used in estimation
#' (typically `floor(T / 2) - 1`) and `nrep` is the number of replicates.
#' @param endpoints_index Vector of length `nbands - 1` or `nbands + 1`
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
#' last indices are NOT included, and the latter case for when they are. The
#' first index is always `1` and the last index is always `nfreq + 1`. These are
#' optional, but if they are included, they must be included for ALL clusters,
#' that is, in each row of `endpoints_index`.
#'
#' @return A matrix of shape (`nrep`, `nbands`) for replicate-specific collapsed
#' measures or a matrix of shape (`nclust`, `nbands`) for cluster-average
#' collapsed measures.
#'
#' @examples
#' # 10 replicates of white noise where T = 400
#' # spec is simulated using the asymptotic distribution of the periodogram
#' # for a white noise process with variance 1
#' length <- 400; nfreq <- floor(length / 2) - 1
#' nrep <- 10
#' spec <- matrix(rchisq(nfreq * nrep, df = 2) / 2, nrow = nfreq)
#' labels <- rep(1:2, each = nrep / 2) # equally divided clusters
#'
#' # first cluster has frequency band endpoints at 50 / T and 100 / T
#' # second cluster has frequency band endpoints at 70 / T and 120 / T
#' endpoints_index <- matrix(c(51, 101, 71, 121), nrow = 2, byrow = T)
#' rep_collapse_by_index(spec, endpoints_index, labels)
rep_collapsed_by_index <- function(spec, endpoints_index, labels = NULL) {
  spec <- as.matrix(spec)
  nfreq <- dim(spec)[1]
  nrep <- dim(spec)[2]
  if (is.vector(endpoints_index)) {
    endpoints_index <- matrix(endpoints_index, nrow = 1, byrow = T)
  }
  nclust <- nrow(endpoints_index)

  # if labels were provided, ensure that an appropriate matrix of endpoints
  # was provided, too
  if (!is.null(labels)) {
    if (length(labels) < nrep) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
                  "not enough labels provided."))
    }
    if (length(labels) > nrep) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
                  "too many labels provided."))
    }
    if (nclust < length(unique(labels))) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
           "labels were provided but not enough sets of frequency bands were ",
           "provided; ", nclust, " sets were given."))
    }
    if (any(unique(labels) < 1) | any(unique(labels) > nclust)) {
      stop(paste0("Error in computing rep. collapsed measures by index; ",
                  "some labels are outside the range 1:nclust."))
    }
  }

  # if no labels were provided but multiple sets of endpoints were provided,
  # compute collapsed measures using the first set and send a warning to user
  if (is.null(labels) & nclust == 1) {
    labels <- rep(1, nrep)
  }
  if (is.null(labels) & (nclust > 1)) {
    warning("Warning in computing collapsed measures by index; no labels were ",
            "provided, but multiple sets of endpoints were provided. Only the ",
            "first set of endpoints were used.")
    endpoints_index <- endpoints_index[1, , drop = FALSE]
    nclust <- 1
  }

  # check that the endpoints are in valid range
  if (any(endpoints_index > (nfreq + 1)) | any(endpoints_index < 1)) {
    stop(paste(c("Error in computing rep. collapsed measures by index; ",
                "endpoints must be between 1 and nrow(spec) but\n",
                endpoints_index, "\nwas provided."), collapse = " "))
  }

  # ensure the endpoints include 1 and (nfreq + 1) before continuing
  if (all(endpoints_index[, 1] != 1)) {
    endpoints_index <- cbind(rep(1, nclust), endpoints_index)
  } else if (any(endpoints_index[, 1] != 1)) {
    stop(paste(c("Error in computing rep. collapsed measures by index; ",
                "ensure that the first provided endpoint of each cluster is 1 ",
                "or the right endpoint of the first band. Got\n",
                endpoints_index, "\n instead."), collapse = " "))
  }
  if (all(endpoints_index[, ncol(endpoints_index)] != (nfreq + 1))) {
    endpoints_index <- cbind(endpoints_index, rep(nfreq + 1, nclust))
  } else if (any(endpoints_index[, ncol(endpoints_index)] != (nfreq + 1))) {
    stop(paste(c("Error in computing rep. collapsed measures by index; ",
                "ensure that the last provided endpoint of each cluster is ",
                "nfreq + 1 or the left endpoint of the last band. Got\n",
                endpoints_index, "\n instead."), collapse = " "))
  }
  for (j in 1:nclust) {
    # ensure endpoints are sorted
    endpoints_index[j,] <- sort(endpoints_index[j,])
  }
  nbands <- ncol(endpoints_index) - 1

  # compute and return collapsed measures
  collapsed <- matrix(nrow = nrep, ncol = nbands)
  for (j in 1:nclust) {
    for (b in 1:nbands) {
      left <- endpoints_index[j, b]
      right <- endpoints_index[j, b + 1] - 1
      band_spec <- spec[left:right, labels == j, drop = F]
      collapsed[labels == j, b] <- colMeans(band_spec)
    }
  }
  return(collapsed)
}

#' @rdname rep_collapsed_by_index
avg_collapsed_by_index <- function(spec, endpoints_index, labels = NULL) {
  spec <- as.matrix(spec)
  nrep <- dim(spec)[2]
  if (is.vector(endpoints_index)) {
    endpoints_index <- matrix(endpoints_index, nrow = 1, byrow = T)
  }
  nclust <- nrow(endpoints_index)

  rep_collapsed <- rep_collapsed_by_index(spec, endpoints_index, labels)
  avg_collapsed <- matrix(nrow = nclust, ncol = ncol(rep_collapsed))
  if (is.null(labels)) labels <- rep(1, nrep)
  for (j in 1:nclust) {
    clust_rep_collapsed <- rep_collapsed[labels == j, , drop = FALSE]
    avg_collapsed[j,] <- colMeans(clust_rep_collapsed)
  }
  return(avg_collapsed)
}

