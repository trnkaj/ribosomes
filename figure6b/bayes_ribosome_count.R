library(rethinking)
library(tidyverse)
library(readxl)

# Load the data
ribosome_counts <- read_excel("./figure6b/ribosome_counts.xlsx")
#select only relevant columns
ribosome_counts <- subset(ribosome_counts, select = c("control", "sirna", "mix_control", "mix_sirna"))

#exploratory statistics
summary(ribosome_counts)
#draw histograms for each treatment in one chart
par(mfrow = c(2, 2))
hist(ribosome_counts$control, main = "Histogram of control ribosome counts", xlab = "Ribosome counts", ylab = "Frequency", col = "#60AF1A", breaks = 100)
hist(ribosome_counts$sirna, col = "#0B7D9E", breaks = 100)
hist(ribosome_counts$mix_control, col = "#A8DC79", breaks = 100)
hist(ribosome_counts$mix_sirna, col = "#00C0BA", breaks = 100)

#reset to one plot
par(mfrow = c(1, 1)) 

#prepare data for analysis
# Convert to two-column format with one column "counts" and the other "treatment" containing the existing column labels
c_df <- data.frame(counts = ribosome_counts$control, treatment = rep("control", length(ribosome_counts$control)))
s_df <- data.frame(counts = ribosome_counts$sirna, treatment = rep("sirna", length(ribosome_counts$sirna)))
mc_df <- data.frame(counts = ribosome_counts$mix_control, treatment = rep("mix_control", length(ribosome_counts$mix_control)))
ms_df <- data.frame(counts = ribosome_counts$mix_sirna, treatment = rep("mix_sirna", length(ribosome_counts$mix_sirna)))
#bind dataframes to one df
ribosome_counts <- rbind(c_df, s_df, mc_df, ms_df)
#clean NA values
ribosome_counts <- ribosome_counts[!is.na(ribosome_counts$counts), ]
#create treatment index column, 1 for control, 2 for sirna, 3 for mix_control, 4 for mix_sirna
ribosome_counts$treatment_id <- ifelse(ribosome_counts$treatment == "control", 1,
                                      ifelse(ribosome_counts$treatment == "sirna", 2,
                                             ifelse(ribosome_counts$treatment == "mix_control", 3, 4)))
#standardize counts
ribosome_counts$std_counts <- standardize(ribosome_counts$counts)

#model
ribosome_counts_m1 <- ulam(
  alist(
    std_counts ~ dnorm(mu, sigma),
    mu <- a[treatment_id],
    a[treatment_id] ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = ribosome_counts,
  chains = 4, cores = 4, iter = 2000
)

precis(ribosome_counts_m1, depth = 2)

post_ribosome_counts_m1 <- extract.samples(ribosome_counts_m1)
#destandardise samples
post_ribosome_counts_m1$a <- post_ribosome_counts_m1$a * sd(ribosome_counts$counts) + mean(ribosome_counts$counts)

#plot posterior distributions for treatments
#export plot to a pdf
pdf("./figure6b/ribosome_counts_plot_posterior.pdf", width = 8, height = 6)
#control
dens(post_ribosome_counts_m1$a[, 1], col = "#60AF1A", xlim = c(60, 240), #main = "Posterior distribution of ribosome counts", xlab = "Number of ribosomes per cell", ylab = "Density", 
     lwd = 2)
#sirna
dens(post_ribosome_counts_m1$a[, 2], col = "#0B7D9E", add = TRUE, lwd = 2)
#mix_control
dens(post_ribosome_counts_m1$a[, 3], col =  "#A8DC79", add = TRUE, lwd = 2)
#mix_sirna
dens(post_ribosome_counts_m1$a[, 4], col = "#00C0BA", add = TRUE, lwd = 2)
dev.off()

#export plot to a pdf
pdf("./figure6b/ribosome_counts_plot_differences.pdf", width = 8, height = 6)
#plot posterior differences
dens(post_ribosome_counts_m1$a[, 3] - post_ribosome_counts_m1$a[, 1], col = "#A8DC79", xlim = c(-70, 60),# main = "Posterior distribution of differences", xlab = "Difference in ribosome counts (mixed - non-mixed)", ylab = "Density", 
     lwd = 2)
dens(post_ribosome_counts_m1$a[, 4] - post_ribosome_counts_m1$a[, 2], col = "#00C0BA", add = TRUE, lwd = 2)
abline(v = 0, col = "black", lwd = 2, lty = 2)
#close pdf device
dev.off()
