generate_simulation_data <- function(model_name, nrep, len) {
  tryCatch({
    return(get(model_name)(nrep, len))
  }, error = function() {
    print("Provided model_name not valid.")
  })
}

#### model 1 ###################################################################
# single step function
# x - frequencies
# breaks - vector of breakpoints; assumed between min(x) and max(x)
# values - function values on the sub-intervals defined by breaks
step_function <- function(x, breaks, values) {
  xmin <- min(x); xmax <- max(x)
  breaks <- sort(c(xmin, xmax, breaks))
  y <- rep(values[1], length(x))
  for(i in 2:(length(breaks) - 1)) {
    y[x >= breaks[i]] <- values[i]
  }
  return(y)
}

# cluster of step functions such that the breaks and values have replicate
# specific uniformly distributed deviations from a marginal value
# nrep - number of time series to generate
# x - frequencies
# breaks - vector of (marginal) breakpoints; assumed between min(x) and max(x)
# values - (marginal) function values on the sub-intervals defined by breaks
# breaks_width - double; controls range of the uniformly distributed noise
#   added to each of the marginal breakpoints defined in breaks
# values_width - double; controls range of the uniformly distributed noise
#   added to each of the marginal breakpoints defined in breaks
clust_step_function <- function(nrep, x, breaks, values, breaks_width,
                                values_width) {
  spec <- matrix(nrow = length(x), ncol = nrep)
  for (k in 1:nrep) {
    breaks_k <- breaks + runif(length(breaks),
                               min = -breaks_width / 2,
                               max = breaks_width / 2)
    values_k <- values + runif(length(values),
                               min = -values_width / 2,
                               max = values_width / 2)
    spec[, k] <- step_function(x, breaks_k, values_k)
  }
  return(spec)
}

# generate data from model 1
# nrep - number of replicates per cluster
# len - length of realizations
model1 <- function(nrep, len) {
  # model parameters
  # locations of marginal breakpoints
  breaks1 <- c(0.1, 0.25) # lower
  breaks2 <- c(0.2, 0.3) # middle
  breaks3 <- c(0.25, 0.4) # upper
  breaks_width <- 0.02

  # marginal values of the constant segments
  values1 <- c(15, 7.5, 1) # lower
  values2 <- c(25, 12.5, 1) # middle
  values3 <- c(35, 17.5, 1) # upper
  values_widths <- c(0.5, 0.75, 0.25)

  # generate theoretical spectra
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = floor(len / 2))
  spec1 <- clust_step_function(nrep, freq, breaks1, values1,
                               breaks_width, values_widths[1])
  spec2 <- clust_step_function(nrep, freq, breaks2, values2,
                               breaks_width, values_widths[2])
  spec3 <- clust_step_function(nrep, freq, breaks3, values3,
                               breaks_width, values_widths[3])
  spec <- cbind(spec1, spec2, spec3)

  # time series realizations
  x <- spec_sim(rbind(spec, spec[dim(spec)[1]:1, ]))

  # MT spectral estimates
  mtout <- sine_multitaper(x, ntapers = floor(sqrt(len / 2)))

  return(list(x = x, labels = labels, freq = freq, spec = spec,
              mtfreq = mtout$mtfreq, mtspec = mtout$mtspec))
}

#### model 2 ###################################################################
# AR(2) theoretical spectrum
# freq - frequencies to evaluate on
# phi1, phi2 - doubles; autoregressive parameters
# sd - positive double; standard deviation of the innovation process
ar2_spec <- function(freq, phi1, phi2, sd) {
  sd^2 / (1 + phi1^2 + phi2^2 -
            2 * phi1 * (1 - phi2) * cos(2 * pi * freq) -
            2 * phi2 * cos(4 * pi * freq))
}

# master function that is used to simulate from the submodels (a), (b), and (c)
# nrep -  number of replicates per cluster
# len - length of time series realizations
# peaks - vector of marginal peak locations in each cluster
# bandwidths - vector of corresponding marginal bandwidths
# sd - marginal innovation standard deviation
# peak_ranges - vector of doubles; controls range of the uniformly distributed
#   noise added to each of the marginal peaks
# bw_ranges - vector of doubles; controls range of the uniformly distributed
#   noise added to each of the marginal bandwidths
# sd_ranges - vector of doubles; controls range of the uniformly distributed
#   noise added to each of the marginal innovation SDs
model2 <- function(nrep, len, peaks, bandwidths, sd, peak_ranges, bw_ranges, sd_ranges) {
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = 250) # true frequencies
  spec <- matrix(nrow = 250, ncol = 3 * nrep) # true spectra
  x <- matrix(nrow = len, ncol = 3 * nrep) # time series data

  for(i in 1:nrep) {
    # realizations of parameters using peak/bw parameterization from
    # Granados-Garcia et al. (2022)
    peaks_ <- peaks + runif(length(peaks), -peak_ranges, peak_ranges)
    bw_ <- bandwidths + runif(length(bandwidths), -bw_ranges, bw_ranges)
    sd_ <- sd + runif(length(sd), -sd_ranges, sd_ranges)
    phi1 <- 2 * cos(2 * pi * peaks_) * exp(-bw_)
    phi2 <- -exp(-2 * bw_)

    # theoretical spectra
    spec[, i] <- ar2_spec(freq, phi1[1], phi2[1], sd_[1])
    spec[, nrep + i] <- ar2_spec(freq, phi1[2], phi2[2], sd_[2])
    spec[, 2*nrep + i] <- ar2_spec(freq, phi1[3], phi2[3], sd_[3])

    # time series data
    x[, i] <- arima.sim(list(ar = c(phi1[1], phi2[1])), n = len, sd = sd_[1])
    x[, nrep + i] <- arima.sim(list(ar = c(phi1[2], phi2[2])), n = len, sd = sd_[2])
    x[, 2 * nrep + i] <- arima.sim(list(ar = c(phi1[3], phi2[3])), n = len, sd = sd_[3])
  }

  # multitaper spectral estimates with ntapers = floor(sqrt(len))
  mtout <- sine_multitaper(x)

  return(list(x = x, labels = labels, freq = freq, spec = spec,
              mtfreq = mtout$mtfreq, mtspec = mtout$mtspec))
}

# model2a - most overlap between peaks
model2a <- function(nrep, len) {
  peaks <- c(0.23, 0.25, 0.27)
  bws <- rep(0.15, 3)
  sd <- rep(2.25, 3)
  peak_ranges <- rep(0.01, 3)
  bw_ranges <- rep(0.005, 3)
  sd_ranges <- rep(0, 3)

  return(model2(nrep = nrep, len = len, peaks = peaks, bandwidths = bws,
                sd = sd, peak_ranges = peak_ranges, bw_ranges = bw_ranges,
                sd_ranges = sd_ranges
  ))
}

# model2b - medium overlap between peaks
model2b <- function(nrep, len) {
  peaks <- c(0.22, 0.25, 0.28)
  bws <- c(0.15, 0.165, 0.18)
  sd <- rep(2.25, 3)
  peak_ranges <- rep(0.01, 3)
  bw_ranges <- rep(0.005, 3)
  sd_ranges <- rep(0, 3)

  return(model2(nrep = nrep, len = len, peaks = peaks, bandwidths = bws,
                sd = sd, peak_ranges = peak_ranges, bw_ranges = bw_ranges,
                sd_ranges = sd_ranges
  ))
}

# model2c - least overlap between peaks
model2c <- function(nrep, len) {
  peaks <- c(0.21, 0.25, 0.29)
  bws <- c(0.15, 0.17, 0.20)
  sd <- rep(2.25, 3)
  peak_ranges <- rep(0.01, 3)
  bw_ranges <- rep(0.005, 3)
  sd_ranges <- rep(0, 3)

  return(model2(nrep = nrep, len = len, peaks = peaks, bandwidths = bws,
                sd = sd, peak_ranges = peak_ranges, bw_ranges = bw_ranges,
                sd_ranges = sd_ranges
  ))
}


#### model 3 ###################################################################
# AR(1) theoretical spectrum
# freq - frequencies to evaluate on
# phi1 - double; autoregressive parameter
# sd - positive double; standard deviation of the innovation process
ar1_spec <- function(x, phi, sd) {
  sd^2 / (1 + phi^2 - 2 * phi * cos(2 * pi * x))
}

# model3 - three clusters of AR(1) processes that mimic gait data
model3 <- function(nrep, len) {
  # model parameters
  phi1_lower <- c(0.78, 0.36, 0.52); phi1_upper <- c(0.82, 0.44, 0.58)
  sd1_lower <- sqrt(c(1.2, 3.8, 1.2)); sd1_upper <- sqrt(c(1.3, 4.5, 1.4))

  # theoretical spectra/ts data
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = 250)
  spec <- matrix(nrow = 250, ncol = 3 * nrep)
  x <- matrix(nrow = len, ncol = 3 * nrep)

  for(i in 1:nrep) {
    # draw realizations of each of the parameters
    phi1_ <- runif(length(phi1_lower), phi1_lower, phi1_upper)
    sd_ <- runif(length(sd1_lower), sd1_lower, sd1_upper)

    # theoretical spectra
    spec[, i] <- ar1_spec(freq, phi1_[1], sd_[1])
    spec[, nrep + i] <- ar1_spec(freq, phi1_[2], sd_[2])
    spec[, 2*nrep + i] <- ar1_spec(freq, phi1_[3], sd_[3])

    # generate time series realization of length len
    x[, i] <- arima.sim(list(ar = phi1_[1]), n = len, sd = sd_[1])
    x[, nrep + i] <- arima.sim(list(ar = phi1_[2]), n = len, sd = sd_[2])
    x[, 2 * nrep + i] <- arima.sim(list(ar = phi1_[3]), n = len, sd = sd_[3])
  }

  ## multitaper spectral estimates with n_tapers = floor(sqrt(len))
  mtout <- sine_multitaper(x)

  return(list(x = x, labels = labels, freq = freq, spec = spec,
              mtfreq = mtout$mtfreq, mtspec = mtout$mtspec))
}

#### simulate from spectrum ####################################################
# simulate time series realization from a single theoretical spectrum
# spec - vector of spectrum values
spec_sim_single <- function(spec) {
  nx <- length(spec)
  sd <- sqrt(1 / (2 * nx))
  z <- vector("complex", nx)
  y <- vector("complex")
  i <- complex(imaginary = 1)

  for (j in 1:nx) {
    if (j / nx == 0.5 || j / nx == 1) {
      z[j] <- rnorm(1, sd = sd)
    } else if (j < floor(nx / 2) + 1) {
      z[j] <- complex(real = rnorm(1, sd = sd), imaginary = rnorm(1, sd = sd))
    } else {
      z[j] <- Conj(z[nx - j])
    }
  }

  for (t in 1:nx) {
    y[t] <- sum(sqrt(spec) * exp(2 * pi * i * (1:nx) * t / nx) * z)
  }

  return(Re(y))
}

# simulate a single time series realization from a each of many theoretical
# spectra
# spec - matrix of spectra (column-wise)
spec_sim <- function(spec) {
  spec <- as.matrix(spec)
  x <- apply(spec, 2, spec_sim_single)
  x <- matrix(x, nrow = nrow(spec), ncol = ncol(spec))
  return(if (ncol(x) == 1) as.vector(x) else x)
}

