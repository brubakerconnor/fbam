#' @importFrom stats runif
mutate <- function(p, pmutate, spec) {
  nfreq <- dim(spec)[1]
  nsubpop <- dim(p)[1]; nbands <- (dim(p)[2] + 1) / 2
  endpoints <- p[, 1:(nbands - 1), drop = FALSE]
  collapsed <- p[, nbands:dim(p)[2], drop = FALSE]

  # endpoints mutation with truncated normal distribution
  for (j in 1:nsubpop) {
    new_endpoints <- c(1, endpoints[j,], nfreq + 1)
    for (l in 2:nbands) {
      # ensure that the current neighboring endpoints to the left and right
      # are not immediately to the left and right; otherwise, proceed with
      # mutation of this endpoint with probability pmutate
      if (!all(diff(new_endpoints[l + -1:1]) == 1) & runif(1) < pmutate) {
        # update widths since the endpoints may have changed in last iteration
        widths <- diff(new_endpoints) - 1
        s <- truncnorm::rtruncnorm(1, a = -widths[l - 1], b = widths[l],
                                   mean = 0, sd = max(widths[(l - 1):l]) / 3)
        new_endpoints[l] <- new_endpoints[l] + round(s)
      }
    }
    endpoints[j,] <- new_endpoints[-c(1, length(new_endpoints))]
  }

  # summary measure mutation with truncated normal distribution
  lottery <- matrix(runif(nsubpop * ncol(collapsed)) < pmutate, nrow = nsubpop)
  if (sum(lottery) > 0) {
    collapsed[lottery] <- truncnorm::rtruncnorm(sum(lottery),
                                                a = 0, b = max(spec),
                                                mean = collapsed[lottery],
                                                sd = max(spec) / 8)
  }

  # return modified individual
  return(as.vector(t(cbind(endpoints, collapsed))))
}
