source("code/Bayes_mix_normal.R")
source("code/Bayes_mix_poisson.R")
source("code/Bayes_mix_student_t.R")
require(ggplot2)

# -------------------------------------------------
# DP Gaussian Mixture Models
# -------------------------------------------------
set.seed(123)
n <- 500
ymix_norm <- rnormmix(n = n, lambda = c(0.6, rep(0.05, 8), 0.1),
mu = seq(1, 20, by = 2), sigma = rep(1, 10))
y <- ymix_norm
plot(ecdf(y))
grid <- seq(min(y) - 1, max(y) + 1, len = 200)
ngrid <- length(grid)


L <- 50
fit_dp_gmm <- dp_gmm(y = y,grid = grid,L = L, amu =  0,
                     b2mu = 0.01,asigma2 = 1,bsigma2 = 1,
                     aalpha = 1,balpha =  1,n_iter = 5000,nburn = 1000)
dp_gmm_dens <- fit_dp_gmm$density
dp_gmm_dens <- fit_dp_gmm$density


# traceplot
set.seed(123)
random_index <- sample(1:ncol(dp_gmm_dens), 1)

plot(dp_gmm_dens[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")

# ess
ess_values_norm <- apply(dp_gmm_dens, 2, coda::effectiveSize)
plot(grid, ess_values_norm, type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)")

# Estimated density
dens3m <- apply(dp_gmm_dens, 2, mean)
dens3l <- apply(dp_gmm_dens, 2, quantile, prob = 0.025)
dens3h <- apply(dp_gmm_dens, 2, quantile, prob = 0.975)

dfhist <- data.frame(y = y)
dfdens3 <- data.frame(dm = dens3m, dl = dens3l, dh = dens3h,
                      seqgrid = grid)
ggplot(dfdens3, aes(x = seqgrid, y = dm)) + geom_line(size = 1, colour = "blue") +
  geom_ribbon(data = dfdens3, aes(x = seqgrid, ymin = dl, ymax = dh),
              alpha = 0.3, fill = "dodgerblue1") + xlab("y") + ylab("Density") + geom_histogram(data = dfhist, aes(x = y,
                                                    y = after_stat(density)), alpha = 0.2, bins = 40, inherit.aes = FALSE,
                                 fill = "gray", colour = "black") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))

# -------------------------------------------------
# DP Poisson Mixture Models - over-dispersed data
# -------------------------------------------------
set.seed(123)
n <- 500
y_van <- numeric(n)
u <- runif(n)
for(i in 1:n){
if(u[i] < 0.4) {y_van[i] <- rpois(1, lambda = 2)} 
else
  if(u[i] < 0.8) {y_van[i] <- rpois(1, lambda = 6)}
   else {y_van[i] <- rpois(1, lambda = 10)}
}
obs_data <- y_van
L = 50
grid <- seq(min(obs_data) - 1, max(obs_data) + 1, by = 1)
res_pois <- dp_pmm_nimble(y = obs_data, grid = grid, L = L, 
                      alambda = 1, blambda = 1, aalpha = 1, 
                      balpha =1, n_iter = 5000, nburn = 1000)
density_pmm <- res_pois$density 

# traceplot
random_index <- sample(1:ncol(density_pmm), 1)
plot(density_pmm[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")

# ESS
ess_values <- apply(density_pmm, 2, coda::effectiveSize)
plot(grid[grid >=0], ess_values[grid >=0], type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)", cex.main = 1.1)
points(grid[grid >=0], ess_values[grid >=0], pch = 16, col = "black")

mean_density <- apply(density_pmm, 2, mean) # Taking mean over each column (1x200)
lower_density <- apply(density_pmm, 2, quantile, prob = 0.025)
upper_density <- apply(density_pmm, 2, quantile, prob = 0.975)

# Estimated density
obs_density <- table(factor(obs_data, levels = grid)) / length(obs_data)  
y_max <- max(c(mean_density, obs_density)) * 1  

plot(grid[grid >=0], mean_density[grid >=0], type = "h", lwd = 2, col = "gray",
     xlab = "y", ylab = "Density", ylim = c(0, y_max),
     main = "DP Poisson Mix Model on Over-dispersed Data")
polygon(c(grid[grid >=0], rev(grid[grid >=0])),
        c(lower_density[grid >=0], rev(upper_density[grid >=0])),
        col = rgb(30/255, 144/255, 255/255, 0.3), border = NA)
offset <- 0.1  
obs_density <- table(factor(obs_data, levels = grid)) / length(obs_data)  
points(grid + offset, obs_density, type = "h", lwd = 2, col = "black", lty = 3)  
legend("topright", legend = c("Mean Density", "95% CI", "Observed Data"),
       col = c("gray", rgb(0, 0, 1, 0.2), "black"), lty = c(1, NA, 3),
       lwd = c(2, NA, 2), fill = c(NA, rgb(30/255, 144/255, 255/255, 0.3), NA), border = NA)

# -------------------------------------------------
# DP Student-t Mixture Models
# -------------------------------------------------
set.seed(123)
# Parameters for the main normal distribution
N <- 100
mu <- 3
sigma <- 1.15

# Generate main normal data
y <- rnorm(N, mu, sigma)

# Proportion of outliers and amount of shift for outliers
prop_outliers <- 0.1  # Higher proportion of outliers for heavier tails
shift <- 20           # Larger shift for more extreme outliers
t_df <- 3             # Degrees of freedom for t-distribution (heavier tails)

# Determine number of outliers
n_outliers <- round(N * prop_outliers)

# Select indices for outliers
ind_outliers <- sample(1:N, size = n_outliers)

# Initialize outlier-adjusted data
y_outliers <- numeric(N)

# Assign non-outlier data
y_outliers[1:(N - length(ind_outliers))] <- y[-ind_outliers]

# Assign outliers using a t-distribution for heavy tails and add shift
y_outliers[(N - length(ind_outliers) + 1):N] <- rt(n_outliers, df = t_df) * sigma + mu + shift

# Final data with outliers and heavier tails
y_data <- y_outliers

L = 50
grid <- seq(min(y_data) - 1, max(y_data) + 1, length.out = 200)

tmix_sample <- DP_tmix_nimble(y_data, grid = grid, L = L, amu = 0, b2mu = 0.01,
                           asigma2 = 1, bsigma2 = 1, rate_nu = 0.1, aalpha = 1, balpha = 1,
                           n_iter = 5000, nburn = 1000)
samples_t_mixture <- tmix_sample$samples
sample_t_dens <- tmix_sample$density

# traceplot
set.seed(123)
random_index <- sample(1:ncol(sample_t_dens), 1)
plot(sample_t_dens[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")

# ess
ess_values <- apply(sample_t_dens, 2, coda::effectiveSize)
plot(grid, ess_values, type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)")
points(grid, ess_values, pch = 16, col = "black")

density <- sample_t_dens
mean_density <- apply(density, 2, mean)
lower_density <- apply(density, 2, quantile, prob = 0.025)
upper_density <- apply(density, 2, quantile, prob = 0.975)
density_df <- data.frame(grid = grid, 
                         mean_density = mean_density, 
                         lower_density = lower_density, 
                         upper_density = upper_density)
df_hist <- data.frame(y = y_data)

# estimated density plot
ggplot() + 
  geom_histogram(data = df_hist, aes(x = y, y = after_stat(density)), bins = 40, alpha = 0.3, fill = "gray", color = "black") +
  geom_line(data = density_df, aes(x = grid, y = mean_density), color = "blue", size = 1) + 
  geom_ribbon(data = density_df, aes(x = grid, ymin = lower_density, ymax = upper_density), alpha = 0.3, fill = "dodgerblue") +
  geom_rug(data = df_hist, aes(x = y), sides = "b", color = "darkblue", alpha = 0.7) + 
  xlab("y") + 
  ylab("Density") + 
  theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))

