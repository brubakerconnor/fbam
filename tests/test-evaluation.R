# informal testing on the accuracy of fitness evaluation
nrep <- 10; n <- 1000
spec_band1 <- matrix(rchisq((n / 2) * nrep, df = 2) / 2, ncol = nrep)
spec_band2 <- matrix(rchisq((n / 2) * nrep, df = 2) / 2, ncol = nrep) + 10
spec <- rbind(spec_band1, spec_band2)
par(mfrow = c(2, 1))
matplot(spec, type = 'l', col = 1, lwd = 0.5)

# plot loss over a grid of two bands and verify that the minimizer occurs at
# lambda_star = 500
lambda_star <- seq(2, dim(spec)[1] - 1, by = 1)
loss <- rep(0, length(lambda_star))
for (i in 1:length(lambda_star)) {
  ch <- matrix(c(lambda_star[i], NA, NA), nrow = 1)
  loss[i] <- clustFBE:::loss_function(ch, rep(1, nrep), spec)
}
plot(lambda_star, loss, type = 'l')
print(which.min(loss))

# do the same thing but now with two clusters and the loss function
# depends on proper assignment of spectra
nrep <- 10; n <- 100
spec_band1 <- matrix(rchisq((3 * n / 4) * nrep, df = 2) / 2, ncol = nrep)
spec_band2 <- matrix(rchisq((n / 4) * nrep, df = 2) / 2, ncol = nrep) + 10
spec1 <- rbind(spec_band1, spec_band2)

spec_band1 <- matrix(rchisq((n / 4) * nrep, df = 2) / 2, ncol = nrep) + 5
spec_band2 <- matrix(rchisq((3 * n / 4) * nrep, df = 2) / 2, ncol = nrep) + 15
spec2 <- rbind(spec_band1, spec_band2)
spec <- cbind(spec1, spec2)

par(mfrow = c(1, 1))
freq <- seq(from = 1/n, by = 1/n, length = dim(spec)[1])
matplot(freq, cbind(spec1, spec2), type = 'l', col = rep(c(2, 4), each = nrep), lwd = 0.5)
abline(v = c(0.25, 0.75), col = c(4, 2), lwd = 2)

# plot loss over a grid of two bands and verify that the minimizer occurs at
# lambda_star = 26 and lambda_star = 76
lambda_star <- seq(2, dim(spec)[1] - 1, by = 1)
loss <- matrix(0, length(lambda_star), length(lambda_star))
for (i in 1:length(lambda_star)) {
  for (j in 1:length(lambda_star)) {
    collapsed1 <- clustFBE:::avg_collapsed_by_index(spec[, 1:nrep], lambda_star[i])
    ch1 <- matrix(c(lambda_star[i], collapsed1), nrow = 1)

    collapsed2 <- clustFBE:::avg_collapsed_by_index(spec[, (nrep + 1):dim(spec)[2]],
                                         lambda_star[j])
    ch2 <- matrix(c(lambda_star[j], collapsed2), nrow = 1)

    loss[i, j] <- 1 / clustFBE:::fitness(rbind(ch1, ch2), spec)
  }
}
heatmap(loss, Colv = NA, Rowv = NA, scale = 'column')
par(mfrow = c(1, 1))
freq <- seq(from = 1/n, by = 1/n, length = dim(spec)[1])
matplot(freq, cbind(spec1, spec2), type = 'l', col = rep(c(2, 4), each = nrep), lwd = 0.5)
abline(v = lambda_star[which(loss == min(loss), arr.ind = T)] / n, col = c(2, 4), lwd = 2)
print(lambda_star[which(loss == min(loss), arr.ind = T)] / n)
