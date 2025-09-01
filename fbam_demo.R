rm(list = ls())
source("fbam_Rfunctions.R")
library(jcolors)
library(tidyverse)
set.seed(5)
cols <- jcolors::jcolors("pal5")
theme_set(theme_light())

# simulate data from a pre-defined simulation model
dat <- model2c(20, 500)

# run GA for a single combination of number of bands/subpopulations
# this could be used when these parameters are known
ga_out <- ga(dat$mtspec, nbands = 3, nsubpop = 3) # input spectra instead of time series
plot(ga_out, which = 1) # show evolution history (avg loss and min loss)
plot(ga_out, which = 2) # show infeasible solutions

# grid search over combination of (nbands, nsubpop)
fbam_out <- fbam(dat$x, nbands = 2:3, nsubpop = 3:4, ncores = 4)
sol <- fbam_out$solution
print(sel_nsubpop <- sol$params$nsubpop)
print(sel_nbands <- sol$params$nbands)

# visualize solutions
par(mfrow = c(1, sel_nsubpop))
for (j in seq(sel_nsubpop)) {
  matplot(dat$mtfreq, dat$mtspec[, sol$labels == j],
          type = 'l', lty = 1, col = cols[j],
          xlab = "frequency", ylab = "power", main = paste0("subpopulation ", j))
  abline(v = dat$mtfreq[sol$cuts[j, -c(1, sel_nbands + 1)]],
         col = cols[j], lwd = 2)
}

# get summary measures
s <- do.call(rbind, lapply(seq(sel_nsubpop), function(j) {
  sm <- summarize_power(dat$mtspec[, sol$labels == j], sol$cuts[j,])
  sm <- data.frame(sm$rep)
  names(sm) <- as.character(seq(sel_nsubpop))
  sm %>%
    pivot_longer(everything(), names_to = "band", values_to = "summary") %>%
    dplyr::mutate(band = factor(band), subpop = paste0("subpop ", j))
}))
s %>%
  ggplot(aes(x = band, y = summary)) +
  geom_boxplot(width = 0.5) +
  facet_wrap(vars(subpop))



