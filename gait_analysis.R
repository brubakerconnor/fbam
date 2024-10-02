source('R/fbam.R')
load('data/gait-processed/gait.rda')
gait_df <- read.csv('data/gait-processed/gait.csv')
gait_df <- gait_df[complete.cases(gait_df),]
set.seed(9)
library(dplyr)

#### time series and spectra plots #############################################
if (!dir.exists('figures')) dir.create('figures')
## original time series against elapsed time
pdf("figures/all_stride_ts_raw.pdf", width = 16, height = 64)
par(mfrow = c(16, 4))
for (i in 1:ncol(gait$x)) {
  plot(gait$xraw[[i]]$elapsed[-1], gait$xraw[[i]]$left, type = 'o',
       main = gait$patient_id[i],
       col = as.numeric(gait$xraw[[i]]$flagged) + 1)
}
dev.off()

## processed time series against elapsed time
pdf("figures/all_stride_ts.pdf", width = 16, height = 64)
par(mfrow = c(16, 4))
for (i in 1:ncol(gait$x)) {
  plot(seq(1, by = 0.5, length = nrow(gait$x)), gait$x[,i],
       type = 'l',
       main = gait$patient_id[i])
}
dev.off()

# Figure 1
pdf('figures/Figure_1.pdf', width = 12, height = 6)
par(mfrow = c(2, 4))
groups <- c("control", "als", "hunt", "park")
group_names <- c("Control", "ALS", "Huntington's", "Parkinson's")
for (i in 1:length(groups)) {
  plot(seq(1, by = 0.5, length = dim(gait$x)[1]),
       gait$x[, which.max(gait$labels == groups[i])],
       type = 'l',
       ylim = c(-5, 8),
       xlab = "seconds", ylab = "seconds",
       main = group_names[i])
  grid()
}
for (i in 1:length(groups)) {
  matplot(gait$mtfreq, gait$mtspec[, gait$labels == groups[i]],
          type = 'l', col = 1, lty = 1, lwd = 0.5,
          ylim = c(0, max(gait$mtspec)),
          xlab = "frequency", ylab = "relative power")
  grid()
}
dev.off()

#### results plot (figure 3) ###################################################
pdf('figures/Figure_3.pdf', width = 10, height = 4)
par(mfrow = c(1, 2))
# matplot(gait$mtfreq, gait$mtspec, type = 'l', col = 1, lty = 1, lwd = 0.3,
#         xlab = "frequency", ylab = "relative power",
#         main = "Setting A")
# abline(v = gait$noclust$endpoints[1])

for (j in 1:gait$clust$params$nclust) {
  matplot(gait$mtfreq, gait$mtspec[, gait$clust$labels == j],
          type = 'l', col = 'grey70', lty = 1, lwd = 0.3,
          ylim = c(0, max(gait$mtspec)),
          xlab = "Frequency", ylab = "Relative Power",
          main = paste0("Subpopulation ", j))

  # vertical line for boundary
  abline(v = gait$clust$endpoints[j, 1], lwd = 2, lty = 2, col = 'red')

  # piecewise constant approximation
  endpointsj <- c(1, gait$clust$endpoints_index[j, ], nrow(gait$mtspec) + 1)
  w <- diff(endpointsj)
  for (l in 1:(ncol(gait$clust$endpoints_index) + 1)) {
    print(length(gait$mtfreq[endpointsj[l]:(endpointsj[l+1] - 1)]))
    print(w[l])
    lines(gait$mtfreq[endpointsj[l]:(endpointsj[l + 1] - 1)],
          rep(gait$clust$collapsed[j, l], w[l]), col = 'black', lwd = 2)
  }

}
dev.off()

#### cluster membership by condition ###########################################
cat('CLUSTER MEMBERSHIP BY CONDITION:\n')
print(table(gait$labels, gait$clust$labels))

#### cluster analysis ##########################################################
cat('\nSUBJECT-LEVEL COVARIATES BY CLUSTER\n')
cat('Cluster 1\n')
cat('#######################################################################\n')
print(summary(gait_df[gait_df$clust_label == 1,
                      c('age', 'height', 'weight', 'gait_speed',
                        'duration_severity_std', 'ssp')]))

cat('Cluster 2\n')
cat('#######################################################################\n')
print(summary(gait_df[gait_df$clust_label == 2,
                      c('age', 'height', 'weight', 'gait_speed',
                        'duration_severity_std', 'ssp')]))

#### wilcoxon sum rank tests ###################################################
cat('\nWILCOXON SUM RANK TESTS\n')
cols <- c('age', 'height', 'weight', 'gait_speed', 'duration_severity_std', 'ssp')
for (col in cols) {
  wt <- wilcox.test(
    gait_df[gait_df$clust_label == 1, col],
    gait_df[gait_df$clust_label == 2, col],
    exact = FALSE
  )
  cat("Variable: ", col, "\t p-value:", round(wt$p.value, 5), "\n")
}

#### kruskal-wallis tests (figure 4) ###########################################
pdf('figures/Figure_4.pdf', width = 10, height = 4)
par(mfrow = c(1, 2))
boxplot(gait_df$noclust_LF ~ gait_df$group,
        xlab = 'group', ylab = 'low frequency power',
        main = 'LF Collapsed Power by Group')
# low frequency power
w <- seq(2, dim(gait$mtspec)[1]) # candidate boundaries
lf_p_values <- rep(0, length(w)) # p-values
for (i in 1:length(w)) {
  # compute replicate specific BAPs computed using current boundary
  rep_bap <- rep_collapsed_by_index(
    gait$mtspec, c(1, w[i], dim(gait$mtspec)[1] + 1)
  )
  # compute Kruskal-Wallis p-value of LF across group
  lf_p_values[i] <- kruskal.test(rep_bap[,1] ~ gait$labels)$p.value
}

plot(gait$mtfreq[w], log(lf_p_values), type = 'l',
     xlab = "boundary", ylab = "logarithm of p-value", main = "Kruskal-Wallis p-values")
abline(v = gait$noclust$endpoints[1], col = 'black', lwd = 2)
abline(h = lf_p_values[gait$noclust$endpoints_index[1]], col = 4)
abline(h = log(0.05), col = "purple", lty = 2, lwd = 2)
dev.off()

#### logistic regression analysis ##############################################
# # binary label: 1 = healthy, 0 otherwise
gait_df$class <- as.factor(ifelse(
  gait_df$group == "hunt", 1, 0
))
gait_df$gender <- as.factor(gait_df$gender)

glm_base <- glm(class ~ age + height + gait_speed, data = gait_df, family = binomial())
glm_cond1 <- glm(class ~ age + height + gait_speed + noclust_LF, data = gait_df, family = binomial())
base_acc <- mean(as.numeric(glm_base$fitted.values > 0.5) == gait_df$class)
cond1_acc <- mean(as.numeric(glm_cond1$fitted.values > 0.5) == gait_df$class)
cat('\nLOGISTIC REGRESSION ANALYSIS\n')
cat(paste0('Accuracy for baseline model: ', round(base_acc, 4), '\n'))
cat(paste0('Accuracy with condition 1 LF measures: ', round(cond1_acc, 4), '\n'))

#### regression analysis #######################################################
cat('\nLINEAR REGRESSION ANALYSIS\n')
gait_df_hd <- gait_df[gait_df$group == 'hunt', ]
lm <- lm(duration_severity ~ age + height + weight + clust_ratio + clust_LF,
         data = gait_df_hd)
print(olsrr::ols_step_best_subset(lm))
print(summary(lm(duration_severity ~ age + clust_ratio,
                 data = gait_df_hd)))

#### variability ###############################################################
mean_mtspec <- mean(as.matrix(gait$mtspec))
totalss <- sum((gait$mtspec - mean_mtspec)^2)

nfreq <- dim(gait$mtspec)[1]
dfreq <- diff(c(1, gait$noclust$endpoints_index[1, ], nfreq + 1))
mean_collapsed <- rep(gait$noclust$collapsed[1, ], dfreq)
mean_collapsed <- matrix(rep(mean_collapsed, ncol(gait$mtspec)),
                       ncol = ncol(gait$mtspec))
noclustss <- sum((gait$mtspec - mean_collapsed)^2)
print(paste0('% variability captured for noclust solution: ',
      round(noclustss / totalss * 100, 4)))

clustss <- 0
for (j in 1:gait$clust$params$nclust) {
  dfreq <- diff(c(1, gait$clust$endpoints_index[j, ], nfreq + 1))
  mean_collapsed <- rep(gait$clust$collapsed[j, ], dfreq)
  mean_collapsed <- matrix(rep(mean_collapsed, ncol(gait$mtspec)),
                         ncol = ncol(gait$mtspec))
  clustss <- clustss + sum((gait$mtspec - mean_collapsed)^2)
}
print(paste0('% variability captured for clust solution: ',
             round(clustss / totalss * 100, 4)))










