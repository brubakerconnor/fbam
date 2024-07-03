#' Mutation of a single solution
#'
#' Mutation of a single solution proceeds in two stages since each cluster
#' center encoded by a single solution consists of Fourier indices for the
#' endpoints and collapsed measures.
#'
#' Fourier indices (endpoints_index) are each mutated with probability pmutate,
#' independently of the rest. If an endpoint is selected for mutation, it is
#' randomly set to a new endpoint that is within a window no wider than
#' `2 * endpoint_range` with the existing endpoint at the center. The
#' existing endpoints are taken into account so that the new endpoint is not
#' sampled outside of them resulting in an unordered solution (though this
#' could probably be relaxed).
#'
#' Collapsed measures are also each mutated with probability pmutate,
#' independently of the rest. If a collapsed measure is selected for mutation,
#' it is randomly shifted by some amount sampled according to a uniform
#' distribution on `-(collapsed_range / 2)` to `collapsed_range / 2`. Any
#' mutated collapsed measures that are mutated to a negative value are
#' set to zero to ensure validity of the mutated solution.
#'
#' @param ch A matrix valued solution of shape (`nclust`, `2 * nbands - 1`).
#' @param pmutate Double; probability of mutation of any single Fourier index
#' or collapsed measure.
#' @param endpoint_range Integer; the total maximum width of the window in which
#' any new endpoint is selected.
#' @param collapsed_range Double; the total width of the window in which any
#' new collapsed measure is shifted.
#' @param nfreq Number of frequencies used in estimation. This is used to
#' ensure new endpoints do not leave the allowed range of Fourier indices.
#'
#' @return A matrix of the same shape as `ch` with the new values of
#' any mutated endpoints or collapsed measures.
mutate <- function(ch, pmutate, endpoint_range, collapsed_range, nfreq) {
  nclust <- dim(ch)[1]; nbands <- (dim(ch)[2] + 1) / 2

  # separate out endpoints and collapsed measures for this solution
  endpoints_index <- ch[, 1:(nbands - 1), drop = FALSE]
  collapsed <- ch[, nbands:dim(ch)[2], drop = FALSE]

  # endpoint mutation
  for (j in 1:nclust) {
    # get the endpoints for this cluster
    endpoints_index_clust <- endpoints_index[j, ]

    # pad the existing endpoints with 1 and `nfreq + 1` to get the
    # appropriate windows for mutation
    endpoints_index_new <- sort(c(1, endpoints_index_clust, nfreq + 1))

    # mutate each endpoint
    for (b in 2:(length(endpoints_index_new) - 1)) {
      # ensure that the current neighboring endpoints to the left and right
      # are not immediately to the left and right; otherwise, proceed with
      # mutation of this endpoint with probability pmutate
      if (!all(diff(endpoints_index_new[b + -1:1]) == 1) & runif(1) < pmutate) {
        # get the left and right endpoints for possible mutation
        # takes into account the existing endpoints so that the new endpoint
        # will not be selected outside these existing endpoints if they are
        # included in the range of mutation specified by endpoint_range
        left <- max(endpoints_index_new[b - 1] + 1,
                    endpoints_index_new[b] - endpoint_range)
        right <- min(endpoints_index_new[b + 1] - 1,
                     endpoints_index_new[b] + endpoint_range)

        # randomly sample the new endpoint from this window
        endpoints_index_new[b] <- sample(seq(left, right), 1)

        # I think this was included because I kept seeing an error and this
        # fixed it
        if(endpoints_index_new[b] == 1) {
          endpoints_index_new[b] <- 2
        }
      }
    }
    # replace existing (internal) endpoints with the new (internal) endpoints
    # sorting for a sanity check to ensure no errors later on in the algorithm
    new_internal <- endpoints_index_new[2:(length(endpoints_index_new) - 1)]
    endpoints_index[j, ] <- sort(new_internal)
  }

  # collapsed measure mutation
  # determine which of the existing collapsed measures should be mutated
  # each collapsed measure is mutated with probability pmutate independently
  # of the other measures
  s <- matrix(runif(nclust * ncol(collapsed)) < pmutate, nrow = nclust)
  collapsed[s] <- collapsed[s] +
    runif(sum(s), min = -(collapsed_range / 2), max = collapsed_range / 2)
  collapsed[collapsed < 0] <- 0 # ensure new measures are non-negative

  # return the mutated and re-combined solution
  return(cbind(endpoints_index, collapsed))
}
