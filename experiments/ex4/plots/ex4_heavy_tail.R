source("code/Bayes_mix_student_t.R")
source("code/Bayes_mix_normal.R")
require(ggplot2)
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


# -------------------------------------------------
# Student-t Mixture Models
# -------------------------------------------------
K =50
alpha = 0.1
grid <- seq(min(y_data) - 1, max(y_data) + 1, length.out = 200)
tmix_sample <- tmix_nimble(y_data, grid = grid, K = K, amu = 0, b2mu = 0.01, 
                           asigma2 = 2, bsigma2 = 2, rate_nu = 0.1, alpha_dirichlet = rep(alpha,K),
                            nsim = 5000, nburn = 1000)
samples_t_mixture <- tmix_sample$samples
sample_t_dens <- tmix_sample$density

# Traceplot
set.seed(123)
random_index <- sample(1:ncol(sample_t_dens), 1)
random_index
plot(sample_t_dens[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")

# ESS
ess_values <- apply(sample_t_dens, 2, coda::effectiveSize)
plot(grid, ess_values, type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)")
points(grid, ess_values, pch = 16, col = "black")


density <- sample_t_dens
mean_density <- apply(density, 2, mean)
lower_density <- apply(density, 2, quantile, prob = 0.025)
upper_density <- apply(density, 2, quantile, prob = 0.975)

# Create a data frame for the densities
density_df <- data.frame(grid = grid, 
                         mean_density = mean_density, 
                         lower_density = lower_density, 
                         upper_density = upper_density)
df_hist <- data.frame(y = y_data)

# Estimated density plot
ggplot() + 
  geom_histogram(data = df_hist, aes(x = y, y = after_stat(density)), bins = 40, alpha = 0.3, fill = "gray", color = "black") +
  geom_line(data = density_df, aes(x = grid, y = mean_density), color = "blue", size = 1) + 
  geom_ribbon(data = density_df, aes(x = grid, ymin = lower_density, ymax = upper_density), alpha = 0.3, fill = "dodgerblue") +
  geom_rug(data = df_hist, aes(x = y), sides = "b", color = "darkblue", alpha = 0.7) + 
  xlab("y") + 
  ylab("Density") + 
  theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 0.43)


# -------------------------------------------------
# Gaussian Mixture Models
# -------------------------------------------------
K =50
alpha = 0.1
fit_norm <- gmm(y = y_data, grid = grid, K = 20, amu = 0, b2mu = 0.01, asigma2 = 0.1,
            bsigma2 = 1, alpha = rep(alpha, 20), nsim = 5000, nburn = 1000)

# traceplot
plot(fit_norm$density[, random_index], type = "l",
     xlab = "Iteration", ylab = "Density")
# ess
ess_values_gmm <- apply(fit_norm$density, 2, coda::effectiveSize)
plot(grid, ess_values_gmm, type = "l", col = "black", lwd = 2,
     xlab = "Grid Points", ylab = "Effective Sample Size (ESS)")
points(grid, ess_values_gmm, pch = 16, col = "black")

dens3m <- apply(fit_norm$density, 2, mean)
dens3l <- apply(fit_norm$density, 2, quantile, prob = 0.025)
dens3h <- apply(fit_norm$density, 2, quantile, prob = 0.975)

df_hist <- data.frame(y = y_data)
density_df <- data.frame(grid = grid, 
                         mean_density = dens3m, 
                         lower_density = dens3l, 
                         upper_density = dens3h)
ggplot() + 
  geom_histogram(data = df_hist, aes(x = y, y = after_stat(density)), bins = 40, alpha = 0.3, fill = "gray", color = "black") +
  geom_line(data = density_df, aes(x = grid, y = mean_density), color = "blue", size = 1) + 
  geom_ribbon(data = density_df, aes(x = grid, ymin = lower_density, ymax = upper_density), alpha = 0.3, fill = "dodgerblue") +
  geom_rug(data = df_hist, aes(x = y), sides = "b", color = "darkblue", alpha = 0.7) + 
  xlab("y") + 
  ylab("Density")+ 
  theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))+ ylim(NA, 0.43)
