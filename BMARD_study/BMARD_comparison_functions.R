# construct a cosine basis using nbasis functions for each estimated spectrum 
# in spec and use these coefficients for clustering in k-means; k-means is
# run with 30 initial starting points
# freq - the frequencies used in estimation
# spec - matrix of spectra (column-wise)
# nclust - number of clusters
# nbasis - number of cosine basis functions
cluster_cosine_basis <- function(freq, spec, nclust, nbasis = 10) {
  # since we are using create.fourier.basis which includes cosine and sine
  # functions plus the constant function
  nbasis <- 2 * nbasis 
  
  # basis expansion on the spectra reflected across the origin 
  # from -0.5 to 0.5
  reflected_mtfreq <- c(-rev(freq), freq)
  reflected_mtspec <- rbind(spec[nrow(spec):1,], spec)
  
  # create cosine basis
  basis <- fda::create.fourier.basis(rangeval = c(-0.5, 0.5), nbasis = nbasis)
  fdPar <- fda::smooth.basis(argvals = seq(-0.5, 0.5, length.out = nrow(reflected_mtspec)),
                        y = reflected_mtspec, fdParobj = basis)
  cos_basis <- fdPar$y2cMap[seq(from = 1, to = nbasis + 1, by = 2),]
  
  # for each spectrum in spec, use least squares to find coefficient of the
  # cosine basis
  cos_coefs <- matrix(nrow = ncol(reflected_mtspec), ncol = (nbasis / 2) + 1)
  # reconstructed <- matrix(nrow = nrow(reflected_mtspec), ncol = ncol(reflected_mtspec))
  for (i in 1:ncol(spec)) {
    df <- data.frame(t(cos_basis))
    df$y <- reflected_mtspec[, i]
    cos_coefs[i,] <- lm(y ~ . - 1, data = df)$coefficients
    # reconstructed[, i] <- colSums(cos_coefs[i, ] * cos_basis)
  }
  # matplot(reflected_mtfreq, reconstructed, type = 'l')
  # plot(fdPar)
  km <- kmeans(cos_coefs, nclust, nstart = 30)
  return(km$cluster)
}

# BMARD sampling
# freqs - Fourier frequencies
# P - list of periodograms
# n_samples - number of samples to obtain from sampler
# sample_freq - length of time series realizations
multirun <- function(freqs, P, n_samples, sample_freq) {
  trunc <- sample(20:30, 1); # random initial level of truncation
  out <- SpectralBDP(Tsize = 500, 
                     omegas = freqs, 
                     Pgram = P, 
                     Nsamp = n_samples, 
                     L = trunc,
                     epsilon = rep(runif(1, .01, .05), trunc),
                     epsilon_BW = runif(1, 0.01, 0.05) , 
                     SamplingRate = sample_freq,
                     MaxFreq = 0.5, 
                     Sup = n_samples, 
                     Lambda_init = sample(1:20, 1), 
                     Lambda_prior = "flat",
                     Lambda_max = 50L, 
                     Lambda_prop_dist = "up1down1", 
                     alpha_const = F,
                     alpha_init = 1, 
                     a_alpha_L = .1, 
                     b_alpha_L = .1, 
                     q = 1)
  return(out)
  # return(list(out = out, ID = P$ID, label = P$label))
}
