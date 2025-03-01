require(mixtools)
source("code/Bayes_mix_normal.R")
set.seed(123)
n <- 500
ymix_norm <- rnormmix(n = n, lambda = c(0.6, rep(0.05, 8), 0.1),
mu = seq(1, 20, by = 2), sigma = rep(1, 10))
y <- ymix_norm
plot(ecdf(y))
grid <- seq(min(y) - 1, max(y) + 1, len = 200)
ngrid <- length(grid)
# Components number (K = 20, 40, 50)
K <- 20
# Alpha (0.1, 0.5, 0.7, 1)
alpha = 1
fit_gmm <- gmm(y = y, grid = grid,  K = K, amu = 5, b2mu = 0.01, asigma2 = 1, 
               bsigma2 = 1, alpha_dirichlet = rep(alpha, K), nsim = 5000, nburn = 1000)

gmm_samples <- fit_gmm$samples
gmm_dens <- fit_gmm$density

# Effective sample size
ess_values_norm <- apply(gmm_dens, 2, coda::effectiveSize)

plot(grid, ess_values_norm, type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)",
     main = "ESS across Grid Points for GMM (alpha = 1) on Data ")
points(grid, ess_values_norm, pch = 16, col = "black")

# Traceplot
set.seed(123)
random_index <- sample(1:ncol(gmm_dens), 1)

plot(gmm_dens[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")

# Weight components plot
fit_gmm_samples <- as.matrix(fit_gmm$samples)
w_samples_20 <- fit_gmm_samples[, grep("^w\\[", colnames(fit_gmm_samples))]
set.seed(2333)  # For reproducibility
random_sample <- sample(1:nrow(w_samples_20), 1)  # random index

weights <- w_samples_20[random_sample, ]  
plot(1:20, weights, type = "o", col = "blue", pch = 16, lwd = 2,
     main = paste(""),
     xlab = "Component", ylab = "Weight")

# Estimated density
dens3m <- apply(gmm_dens, 2, mean)
dens3l <- apply(gmm_dens, 2, quantile, prob = 0.025)
dens3h <- apply(gmm_dens, 2, quantile, prob = 0.975)

dfhist <- data.frame(y = y)
dfdens3 <- data.frame(dm = dens3m, dl = dens3l, dh = dens3h,
                      seqgrid = grid)
ggplot(dfdens3, aes(x = seqgrid, y = dm)) + geom_line(size = 1, colour = "blue") +
  geom_ribbon(data = dfdens3, aes(x = seqgrid, ymin = dl, ymax = dh),
              alpha = 0.3, fill = "dodgerblue1") + xlab("y") + ylab("Density") + geom_histogram(data = dfhist, aes(x = y,
                                                    y = after_stat(density)), alpha = 0.2, bins = 40, inherit.aes = FALSE,
                                 fill = "gray", colour = "black")+ theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))
