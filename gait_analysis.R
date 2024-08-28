source('R/fbam.R')
load('data/gait-processed/gait.rda')
gait_df <- read.csv('data/gait-processed/gait.csv')
gait_df <- gait_df[complete.cases(gait_df),]
set.seed(9)

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

#### results plot (figure 2) ###################################################
pdf('figures/Figure_2.pdf', width = 12, height = 4)
par(mfrow = c(1, 3))
matplot(gait$mtfreq, gait$mtspec, type = 'l', col = 1, lty = 1, lwd = 0.3,
        xlab = "frequency", ylab = "relative power",
        main = "Setting A")
abline(v = gait$noclust$endpoints[1])

for (j in 1:gait$clust$params$nclust) {
  matplot(gait$mtfreq, gait$mtspec[, gait$clust$labels == j],
          type = 'l', col = 1, lty = 1, lwd = 0.3,
          ylim = c(0, max(gait$mtspec)),
          xlab = "frequency", ylab = "relative power",
          main = paste0("Setting B; Cluster ", j))
  abline(v = gait$clust$endpoints[j, 1])
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

#### kruskal-wallis tests (figure 3) ###########################################
pdf('figures/Figure_3.pdf', width = 12, height = 4)
par(mfrow = c(1, 3))
boxplot(gait_df$noclust_LF ~ gait_df$group,
        xlab = 'Group', ylab = 'LF Collapsed Power',
        main = 'LF Collapsed Power by Group (Setting A)')
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

plot(gait$mtfreq[w], lf_p_values, type = 'l', ylim = c(0, 1),
     xlab = "Boundary", ylab = "p-value", main = "LF Power")
abline(v = gait$noclust$endpoints[1], col = 4)
abline(h = lf_p_values[gait$noclust$endpoints_index[1]], col = 4)
abline(h = 0.05, col = "purple", lty = 2)

# high frequency power
hf_p_values <- rep(0, length(w)) # p-values
for (i in 1:length(w)) {
  # compute replicate specific BAPs computed using current boundary
  rep_bap <- rep_collapsed_by_index(
    gait$mtspec, c(1, w[i], dim(gait$mtspec)[1] + 1)
  )
  # compute Kruskal-Wallis p-value of LF across group
  hf_p_values[i] <- kruskal.test(rep_bap[,2] ~ gait$labels)$p.value
}

plot(gait$mtfreq[w], hf_p_values, type = 'l', ylim = c(0, 1),
     xlab = "Boundary", ylab = "p-value", main = "HF Power")
# abline(v = gait$mtfreq[which.min(hf_p_values)], col = 2)
# abline(h = min(hf_p_values), col = 2)
abline(v = gait$noclust$endpoints[1], col = 4)
abline(h = hf_p_values[gait$noclust$endpoints_index[1]], col = 4)
abline(h = 0.05, col = "purple", lty = 2)
dev.off()
