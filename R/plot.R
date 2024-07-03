#' Plot GA evolution history
#'
#' `plot_ga()` displays the evolution history in terms of average and minimum
#' values of the loss function within each generation of a GA run.
#' `plot_ga_solution()` displays a plot of the spectra colored according to
#' the label they were assigned along with the estimated frequency bands in
#' each cluster.
#'
#' @param out The list object that it outputted by \code{\link{ga}}
#'
#' @export
#'
#' @examples
#' data(model1)
#' out <- ga(model1$x, 3, 3)
#' plot_ga(out)
plot_ga <- function(out) {
  numgen <- length(out$avgfit)
  df <- data.frame(gen = 1:numgen,
                   avgloss = 1 / out$avgfit, minloss = 1 / out$maxfit,
                   ninfeasible = out$ninfeasible)
  plot(df$gen, df$avgloss,
       type = 'l', col = 1, lwd = 1.5,
       ylim = c(min(df$minloss), max(df$avgloss)),
       xaxt = 'n',
       xlab = "Generation Number", ylab = "Loss Function",
       main = "GA Evolution History")
  axis(1, at = floor(seq(0, numgen, length = 10)))
  grid()
  lines(df$gen, df$minloss, type = 'l', col = 2, lwd = 1.5)
  legend("topright", legend = c("Average Loss", "Minimum Loss"),
         lwd = rep(1.5, 2), col = 1:2)
}

#' @rdname plot_ga
#' @export
plot_ga_solution <- function(out) {
  n <- dim(out$X)[1]; nrep <- dim(out$X)[2]
  spec <- out$spec; nfreq <- dim(spec)[1]
  mtfreq <- seq(from = 1/n, by = 1/n, length = nfreq)

  nclust <- out$params$nclust; nbands <- out$params$nbands
  endpoints <- out$endpoints; collapsed <- out$collapsed
  endpoints_index <- endpoints * dim(out$X)[1]
  nendpoints <- nbands - 1

  if (nclust == 2) old.par <- par(mfrow = c(1, 2))
  if (nclust > 2) old.par <- par(mfrow = c(ceiling(nclust / 3), 3))
  for (j in 1:nclust) {
    matplot(mtfreq,
            spec[, out$labels == j, drop = F],
            type = "l", lty = 1, lwd = 0.2, col = "black",
            ylim = c(0, max(spec)),
            xlab = "Frequency", ylab = "Power",
            main = paste0("Cluster Label ", j)
    )
    abline(v = endpoints[j, ], col = 1, lwd = 2, lty = 1)
  }
  if (nclust >= 2) par(old.par)
}

#' #' Visualize spectra
#' #'
#' #' @param spec
#' #' @param freq
#' #' @param endpoints
#' #' @param labels
#' #'
#' #' @return
#' #'
#' #' @examples
#' plot_spec <- function(spec, freq = NULL, endpoints = NULL, labels = NULL) {
#'   # [ TODO ] change color palette?
#'   if (is.vector(spec)) spec <- as.matrix(spec)
#'   if (is.null(freq)) {
#'     n <- 2 * dim(spec)[1] + 1
#'     freq <- seq(from = 1/n, by = 1/n, length = dim(spec)[1])
#'   }
#'   if (is.null(labels)) labels <- rep(1, dim(spec)[2])
#'   if (is.vector(endpoints)) endpoints <- matrix(endpoints, nrow = 1)
#'
#'   matplot(freq, spec, type = 'l', lty = 1, col = labels,
#'           xlab = "Frequency", ylab = "Power")
#'   if (!is.null(endpoints)) {
#'     for (j in 1:dim(endpoints)[1]) {
#'       abline(v = endpoints[j,], col = j)
#'     }
#'   }
#' }
