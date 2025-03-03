source("code/Bayes_mix_poisson.R")
source("code/Bayes_mix_normal.R")
require(ggplot2)
require(coda)
require(rjags)

# Underdispersed data
set.seed(123)
counts_under <- rbinom(500, 10, 0.8)
obs_data <- counts_under

# -------------------------------------------------
# Poisson Mixture Models
# -------------------------------------------------
# Change to desired K components
K = 20
# Change to desired alpha parameter 
alpha_p = 1
grid <- seq(min(obs_data) - 1, max(obs_data) + 1, by = 1)
res_pois <- pmm(y = obs_data, grid = grid, K = K, alpha = 1, beta = 1, alpha_dir = rep(alpha_p,K))
density_pmm <- res_pois$density 

# ESS
ess_values <- apply(density_pmm, 2, coda::effectiveSize)
plot(grid, ess_values, type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)")
points(grid, ess_values, pch = 16, col = "black")

# Traceplot
set.seed(123)  
random_index <- sample(1:ncol(density_pmm), 1)

# Plot the trace of the density at this selected grid point
plot(density_pmm[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")

# Estimated density
mean_density <- apply(density_pmm, 2, mean) 
lower_density <- apply(density_pmm, 2, quantile, prob = 0.025)
upper_density <- apply(density_pmm, 2, quantile, prob = 0.975)

density_df <- data.frame(
  grid = grid,
  mean_density = mean_density,
  lower_density = lower_density,
  upper_density = upper_density
)

# plot
offset <- 0.1
obs_density <- table(factor(obs_data, levels = grid)) / length(obs_data) 
y_max <- max(c(mean_density, obs_density)) * 1
plot(grid, mean_density, type = "h", lwd = 2, col = "gray",
     xlab = "y", ylab = "Density", ylim = c(0, y_max),
     main = paste("Poisson Mixture Model Density for K =", K))

# Add shaded area for the 95% credible interval
polygon(c(grid, rev(grid)),
        c(lower_density, rev(upper_density)),
        col = rgb(30/255, 144/255, 255/255, 0.3), border = NA)
points(grid + offset, obs_density, type = "h", lwd = 2, col = "black", lty = 3)  
legend("topright", legend = c("Mean Density", "95% CI", "Observed Data"),
       col = c("gray", rgb(0, 0, 1, 0.2), "black"), lty = c(1, NA, 3),
       lwd = c(2, NA, 2), fill = c(NA, rgb(30/255, 144/255, 255/255, 0.3), NA), border = NA)


# -------------------------------------------------
# Square root Gaussian Mixture Models
# -------------------------------------------------
grid_norm <- seq(min(obs_data) - 1, max(obs_data) + 1, len = 2000)
K = 20
fit_norm <- gmm(y = sqrt(obs_data), grid = grid_norm, K = K, amu = 0, b2mu = 100, asigma2 = 0.1,
            bsigma2 = 1, alpha = rep(alpha_p, K), nsim = 5000, nburn = 1000)
density_sqrt_norm <- fit_norm$density

# ESS
ess_values_sqrt_norm <- apply(density_sqrt_norm, 2, coda::effectiveSize)
plot(grid_norm, ess_values_sqrt_norm, type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)")
points(grid_norm, ess_values_sqrt_norm, pch = 16, col = "black")

# Traceplot
set.seed(123)
random_index <- sample(1:ncol(fit_norm$density), 1)

plot(fit_norm$density[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")

dens3m <- apply(fit_norm$density, 2, mean)
dens3l <- apply(fit_norm$density, 2, quantile, prob = 0.025)
dens3h <- apply(fit_norm$density, 2, quantile, prob = 0.975)

dfhist <- data.frame(y = sqrt(obs_data))
dfdens3 <- data.frame(dm = dens3m, dl = dens3l, dh = dens3h,
                      seqgrid = grid_norm)
ggplot(dfdens3, aes(x = seqgrid, y = dm)) + geom_line(size = 1, colour = "blue") +
  geom_ribbon(data = dfdens3, aes(x = seqgrid, ymin = dl, ymax = dh),
              alpha = 0.3, fill = "dodgerblue1") + xlab("y") + ylab("Density") + geom_histogram(data = dfhist, aes(x = y,
                                                    y = after_stat(density)), alpha = 0.2, bins = 40, inherit.aes = FALSE,
                                 fill = "gray", colour = "black") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))
