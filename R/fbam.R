#' Frequency band analysis for multiple stationary time series
#'
#' Simultaneous segmentation and frequency band estimation for a collection of
#' independent time series realizations. When either or both of the number of
#' frequency bands or subpopulations is unknown, the final solution is chosen
#' by minimization of validation criteria.
#'
#' @param X Column-wise matrix of replicate time series. A vector is treated as a single replicate.
#' @param nbands Possible values for the number of frequencies bands. Must all be greater than or equal to 2.
#' @param nsubpop Possible values for the number of subpopulations.
#' If more than one value is supplied, they all must be greater than or equal to 2.
#' Default value is \code{1}.
#' @param popsize Population size. Default is \code{50}.
#' @param pmutate Probability of mutation. Default value is \code{0.1}.
#' @param maxgen Maximum number of generations to run GA. Default value is \code{500}.
#' @param maxrun Maxmimum number of generations without improvement (defined below) before the GA is
#' terminated. Default value is \code{100}.
#' @param tol Tolerance used in determining improvement. Default value is \code{1e-2} which corresponds to 5% improvement.
#' @param ntapers Number of tapers used in multitaper estimates from \code{X}.
#' Default value is \code{floor(sqrt(nrow(as.matrix(X))))}, the floor of the square root
#' of the length of the time series.
#' @param parallel Number of cores to use in parallelization. Default value is \code{1} (no parallelization).
#'
#'
#' @return A list with the following components:
#' \itemize{
#'  \item \code{selected_solution}: The selected solution based on minimum validation criterion which is
#'  itself a list with the following components:
#'  \itemize{
#'    \item \code{spec}: Column-wise multitaper estimates in correspondence with the columns of \code{X}.
#'    \item \code{labels}: Subpopulation assignments.
#'    \item \code{endpoints}: Row-wise matrix of subpopulation-specific frequency-valued band boundaries.
#'    \item \code{avg_summary}: Row-wise matrix of subpopulation-specific average summary measures.
#'    \item \code{rep_summary}: Matrix of replicate-specific summary measures whose rows are in correspondence with
#'    the columns of \code{X}.
#'    \item \code{objective}: The objective function value of the final solution.
#'    \item \code{validation}: The value of the validation criteria.
#'    \item \code{avgfit}, \code{maxfit}, \code{ninfeasible}: The average fitness value, maximum fitness value, and
#'    number of infeasible solutions at the end of each generation. Vectors of length equal to the number of generations
#'    completed before termination of the GA.
#'    \item \code{params}: List of parameters sent to the GA.
#'  }
#'  \item \code{all_solutions}: A list of solutions of the same list structure as \code{selected_solution.}
#' }
#' In the case that \code{length(nsubpop) == 1} and \code{length(nbands) == 1}, \code{selected_solution}
#' and \code{all_solutions} are the same object. This is done intentionally for consistent output structure.
#' @export
#'
#' @examples
#' X <- matrix(nrow = 500, ncol = 20)
#' for (i in 1:20) X[,i] <- arima.sim(list(ar = runif(1, 0.2, 0.8)), n = 500)
#' sine_mt(X)
fbam <- function(X, nbands, nsubpop = 1, popsize = 50, pmutate = 0.1,
                 maxgen = 500, maxrun = 100, tol = 1e-2,
                 ntapers = floor(sqrt(nrow(as.matrix(X)))),
                 parallel = 1) {
  param_grid <- expand.grid(nbands = nbands, nsubpop = nsubpop)
  all_solutions <- parallel::mclapply(1:nrow(param_grid), function(i) {
    nsubpop <- param_grid$nsubpop[i]
    nbands <- param_grid$nbands[i]
    ga(X, nbands, nsubpop, popsize, pmutate, maxgen, maxrun, tol, ntapers,
       parallel == 1)
  }, mc.cores = parallel)
  validation <- unlist(lapply(all_solutions, function(x) x$validation))
  selected <- all_solutions[[which.min(validation)]]
  return(list(selected_solution = selected, all_solutions = all_solutions))
}
