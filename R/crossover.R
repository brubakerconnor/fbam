#' Crossover
#'
#' Perform crossover operation between two parents. Crossover is performed by
#' randomly swapping complimentary subsets of cluster solutions between the
#' two parents.
#'
#' @param p1 Parent 1 - single solution from population
#' @param p2 Parent 2 - single solution from population
#'
#' @return A list of containing the byproducts of crossover (children) with
#' names `c1` and `c2`.
crossover <- function(p1, p2) {
  nclust <- dim(p1)[1]

  sl <- runif(nclust) < 0.5 # which clusters should be swapped?
  # swap the clusters as determined above between p1 and p2
  c1 <- rbind(p1[seq_len(nclust)[sl], ], p2[seq_len(nclust)[!sl], ])
  c2 <- rbind(p1[seq_len(nclust)[!sl], ], p2[seq_len(nclust)[sl], ])
  return(list(c1 = c1, c2 = c2))
}
