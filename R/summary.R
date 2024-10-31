rep_summary <- function(spec, endpoints) {
  nfreq <- dim(spec)[1]; nrep <- dim(spec)[2]
  nbands <- length(endpoints) + 1
  endpoints <- c(1, endpoints, nfreq + 1)
  summary <- matrix(nrow = nrep, ncol = nbands)
  for (l in 1:nbands) {
    band_spec <- spec[endpoints[l]:(endpoints[l + 1] - 1), , drop = F]
    summary[,l] <- colMeans(band_spec)
  }
  return(summary)
}

avg_summary <- function(spec, endpoints) {
  rep <- rep_summary(spec, endpoints)
  return(colMeans(rep))
}
