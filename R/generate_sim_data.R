generate_simulation_data <- function(model_name, nrep, len) {
  tryCatch({
    return(get(model_name)(nrep, len))
  }, error = function() {
    print("Provided model_name not valid.")
  })
}

#### model 1 ###################################################################
# piecewise smooth step function
# x         inputs to at which to evaluate function
# values    vector of length L
# breaks    vector of length L-1
# delta     double, spacing around breaks for cubic spline
pws_step <- function(x, values, breaks, delta = 0.025) {
  coefs <- matrix(nrow = length(breaks), ncol = 4)
  for (i in 1:length(breaks)) {
    left <- (breaks[i] - delta/2); right <- (breaks[i] + delta/2)
    A <- matrix(c(
      left^3, left^2, left, 1,
      right^3, right^2, right, 1,
      3 * left^2, 2 * left, 1, 0,
      3 * right^2, 2 * right, 1, 0
    ), nrow = 4, ncol = 4, byrow = TRUE)
    coefs[i, ] <- solve(A, c(values[i], values[i + 1], 0, 0))
  }
  y <- rep(values[length(values)], length(x))
  marks <- sort(c(0, breaks - (delta/2), breaks + (delta/2), 0.5))
  for (i in 1:length(breaks)) {
    ind <- 2 * (i - 1) + 1
    y[x >= marks[ind] & x < marks[ind + 1]] <- values[i]
    x_gap <- x[x >= marks[ind + 1] & x < marks[ind + 2]]
    y[x %in% x_gap] <- matrix(c(
      x_gap^3, x_gap^2, x_gap, rep(1, length(x_gap))
    ), ncol = 4) %*% coefs[i, ]
  }
  return(y)
}

# generate data from model 1
# nrep - number of replicates per cluster
# len - length of realizations
model1 <- function(nrep, len) {
  # model parameters
  # locations of marginal breakpoints
  breaks <- matrix(c(
    0.1, 0.25, #lower
    0.2, 0.3, #middle
    0.25, 0.4 #upper
  ), nrow = 3, byrow = TRUE)
  delta1 <- 0.025  # spacing for cubic spline

  # marginal values of the constant segments
  values <- matrix(c(
    15, 7.5, 2, #lower
    25, 12.5, 4, #middle
    35, 17.5, 6 #upper
  ), nrow = 3, byrow = TRUE)
  delta2 <- 2

  # generate theoretical spectra
  labels <- rep(1:3, each = nrep)
  freq <- seq(0, 0.5, length = floor(len / 2))
  # spec1 <- clust_step_function(nrep, freq, breaks1, values1,
  #                              breaks_width, values_widths[1])
  # spec2 <- clust_step_function(nrep, freq, breaks2, values2,
  #                              breaks_width, values_widths[2])
  # spec3 <- clust_step_function(nrep, freq, breaks3, values3,
  #                              breaks_width, values_widths[3])
  # spec <- cbind(spec1, spec2, spec3)
  spec <- matrix(nrow = floor(len / 2), ncol = 3 * nrep)
  for (j in 1:3) {
    for (k in 1:nrep) {
      spec[, nrep * (j - 1) + k] <- pws_step(
        freq,
        values[j, ] + runif(ncol(values), min = -delta2, max = delta2),
        breaks[j, ],
        delta1
      )
    }
  }

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
  peaks <- c(0.21, 0.25, 0.29)
  bws <- c(0.155, 0.15, 0.155)
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
  peaks <- c(0.19, 0.25, 0.31)
  bws <- c(0.165, 0.15, 0.165)
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

