source("code/Bayes_mix_student_t.R")
source("code/Bayes_mix_normal.R")
# repeat experiment 20 times
set.seed(123)
# Parameters for the main normal distribution
N <- 100
mu <- 3
sigma <- 1.15
# Prop outliers: 0.01, 0.05, 0.10
prop_outliers <- 0.1  # Higher proportion of outliers for heavier tails
shift <- 20           # Larger shift for more extreme outliers
t_df <- 3             # Degrees of freedom for t-distribution (heavier tails)
num_datasets <- 20    # Number of datasets to generate

# Initialize a list to store each dataset
datasets <- vector("list", num_datasets)

# Generate 20 datasets with outliers and store them in the list
for (i in 1:num_datasets) {
  # Generate the main normal data
  y <- rnorm(N, mu, sigma)
  
  # Determine the number of outliers
  n_outliers <- round(N * prop_outliers)
  
  # Select indices for outliers randomly
  ind_outliers <- sample(1:N, size = n_outliers)
  
  # Initialize outlier-adjusted data
  y_outliers <- numeric(N)
  
  # Assign non-outlier data
  y_outliers[1:(N - length(ind_outliers))] <- y[-ind_outliers]
  
  # Assign outliers using a t-distribution for heavy tails and add shift
  y_outliers[(N - length(ind_outliers) + 1):N] <- rt(n_outliers, df = t_df) * sigma + mu + shift
  
  # Store the generated dataset in the list
  datasets[[i]] <- y_outliers
}

# looping 20 times for t mix model
set.seed(123)
# Set parameters
num_experiments <- length(datasets)  # Number of datasets (20 in this case)
K <- 20  # Number of mixture components; K = 20, 50

# Initialize vectors to store WAIC and LPML values for each experiment
waic_values <- numeric(num_experiments)
lpml_values <- numeric(num_experiments)

# Loop through each dataset to fit the t-mixture model and calculate good fit criteria
for (i in 1:num_experiments) {
  grid <- seq(min(datasets[[i]]) - 1, max(datasets[[i]]) + 1, len = 200) 
  # Fit the t-mixture model on the current dataset
  tmix_sample <- tmix_nimble(datasets[[i]], grid = grid, K = K, amu = 0, b2mu = 0.01,
                             asigma2 = 1, bsigma2 = 1, rate_nu = 0.1, alpha_dirichlet = rep(1, K),
                             nsim = 5000, nburn = 1000)
  # Extract posterior samples
  samples_t_mixture <- tmix_sample$samples
  samples_matrix_t_mixture <- as.matrix(samples_t_mixture)
  
  # Extract weights, means, variances, and degrees of freedom samples
  w_samples_t <- samples_matrix_t_mixture[, grep("w", colnames(samples_matrix_t_mixture))]
  mu_samples_t <- samples_matrix_t_mixture[, grep("mu", colnames(samples_matrix_t_mixture))]
  tau_samples_t <- samples_matrix_t_mixture[, grep("tau", colnames(samples_matrix_t_mixture))]
  df_samples_t <- samples_matrix_t_mixture[, grep("nu", colnames(samples_matrix_t_mixture))]
  sigma2_samples_t <- 1 / tau_samples_t
  
  # Calculate WAIC and LPML for the current experiment
  ghc_t <- good_fit_criteria_t(y = datasets[[i]], P = w_samples_t, Mu = mu_samples_t,
                               Sigma2 = sigma2_samples_t, df_samples = df_samples_t, nburn = 1000)
  
  # Store the WAIC and LPML values
  waic_values[i] <- ghc_t$WAIC
  lpml_values[i] <- ghc_t$LPML
}

# Calculate the average WAIC and LPML across all experiments
average_waic <- mean(waic_values)
average_lpml <- mean(lpml_values)
sd_lpml <- sd(lpml_values)
sd_waic <- sd(waic_values)

# Df results
results <- data.frame(dataset_id = 1:num_experiments, WAIC = waic_values, LPML = lpml_values)
print(results)
write.csv(results, file = "experiments/t_mix/prop01/experiment_heavy_t_mix_results_k20_alpha1.csv", row.names = FALSE)

# Df summary
summary <- data.frame(Average_WAIC = average_waic, Average_LPML = average_lpml,
                      SD_WAIC = sd_waic, SD_LPML <- sd_lpml)
write.csv(summary, file = "experiments/t_mix/prop01/summary_heavy_t_mix_k20_alpha1.csv", row.names = FALSE)

#------------------------------------------------------------------------------#
# Looping for gmm on heavy tail data
set.seed(123)
num_experiments <- length(datasets)  # Number of datasets (20 in this case)
K <- 20  # Number of mixture components
# Initialize vectors to store WAIC and LPML values for each experiment
waic_values2 <- numeric(num_experiments)
lpml_values2 <- numeric(num_experiments)

# Loop through each dataset to fit the Gaussian mixture model and calculate good fit criteria
for (i in 1:num_experiments) {
  grid <- seq(min(datasets[[i]]) - 1, max(datasets[[i]]) + 1, len = 200)
  # Fit the Gaussian mixture model on the current dataset
  fit_norm <- gmm(y = datasets[[i]], grid = grid, K = K, amu = 0, b2mu = 0.01, 
                  asigma2 = 1, bsigma2 = 1, alpha_dirichlet = rep(1, K), nsim = 5000, nburn = 1000)
  
  # Extract posterior samples for weights, means, and variances
  samples_norm <- fit_norm$samples
  samples_matrix_norm <- as.matrix(samples_norm)
  
  # Extract weights, means, and variances
  w_samples <- samples_matrix_norm[, grep("w", colnames(samples_matrix_norm))]
  mu_samples <- samples_matrix_norm[, grep("mu", colnames(samples_matrix_norm))]
  tau_samples <- samples_matrix_norm[, grep("tau", colnames(samples_matrix_norm))]
  sigma2_samples <- 1 / tau_samples
  
  # Calculate WAIC and LPML for the current experiment
  ghc_norm <- good_fit_criteria(y = datasets[[i]], P = w_samples, Mu = mu_samples, 
                                Sigma2 = sigma2_samples, nburn = 1000)
  
  # Store the WAIC and LPML values
  waic_values2[i] <- ghc_norm$WAIC
  lpml_values2[i] <- ghc_norm$LPML
}

# Calculate the average WAIC and LPML across all experiments
average_waic2 <- mean(waic_values2)
average_lpml2 <- mean(lpml_values2)
sd_lpml2 <- sd(lpml_values2)
sd_waic2 <- sd(waic_values2)

# Df results
results_gmm <- data.frame(dataset_id = 1:num_experiments, WAIC = waic_values2, LPML = lpml_values2)
print(results_gmm)
write.csv(results_gmm, file = "experiments/t_mix/prop01/experiment_heavy_gmm_results_k50_alpha1.csv", row.names = FALSE)

# Df summary
summary_gmm <- data.frame(Average_WAIC = average_waic2, Average_LPML = average_lpml2,
                      SD_WAIC = sd_waic2, SD_LPML <- sd_lpml2)
write.csv(summary_gmm, file = "experiments/t_mix/prop01/summary_heavy_gmm_k50_alpha1.csv", row.names = FALSE)


