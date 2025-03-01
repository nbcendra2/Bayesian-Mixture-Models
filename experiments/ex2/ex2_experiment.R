source("code/Bayes_mix_normal.R")
source("code/Bayes_mix_poisson.R")

# Repeating experiments 20 times Overdispersed Poisson Mixture

# Define the number of experiments
set.seed(123)
num_datasets <- 20
# Initialize a list to store each dataset
datasets_overdispersed <- vector("list", num_datasets)

# Set Parameters
n <- 500
K <- 20
for (i in 1:num_datasets) {
  # Generate synthetic over-dispersed data (if you already have obs_data, skip this part)
  obs_data <- numeric(n)
  u <- runif(n)
  for (j in 1:n) {
    if (u[j] < 0.4) {
      obs_data[j] <- rpois(1, lambda = 2)
    } else if (u[j] < 0.8) {
      obs_data[j] <- rpois(1, lambda = 6)
    } else {
      obs_data[j] <- rpois(1, lambda = 10)
    }
  }
  
  datasets_overdispersed[[i]] <- obs_data
}

# Examining datasets


# Number of datasets to be iterated
num_experiments <- length(datasets_overdispersed)

# Results df for LPML and WAIC
results <- data.frame(dataset_id = integer(num_experiments), WAIC = numeric(num_experiments),
                       LPML = numeric(num_experiments))

# Loop to perform experiment 20 times for overdispersed

for (i in 1:num_experiments) {
  # alpha and beta prior maybe can use method of moment to determine the best results?
  # Fit Poisson mixture model with current dataset
  res_pois <- pmm(y = datasets_overdispersed[[i]], 
                  grid = seq(min(datasets_overdispersed[[i]]) - 1, 
                             max(datasets_overdispersed[[i]]) + 1, by = 1), 
                  K = K, alpha = 1, beta = 1, alpha_dir = rep(0.1, K))
  
  # Extract samples from the fitted model
  samples_mat_pois <- as.matrix(res_pois$samples)
  lambda_samples <- samples_mat_pois[, 1:K]
  w_samples <- samples_mat_pois[, (K + 1):(2 * K)]
  
  # Calculate good fit criteria
  gh3_pois <- good_fit_criteria_pois(datasets_overdispersed[[i]], w_samples, lambda_samples, nburn = 1000)
  
  # Store results in a data frame
  results$dataset_id[i] <- i
  results$LPML[i] <- gh3_pois$LPML
  results$WAIC[i] <- gh3_pois$WAIC
}

# Calculate average LPML and WAIC
average_LPML <- mean(results$LPML)
average_WAIC <- mean(results$WAIC)
sd_LPML <- sd(results$LPML)
sd_WAIC <- sd(results$WAIC)

# Display the results
print(results)
write.csv(results, file = "experiments/Poisson_mix/overdispersed_dataset/experiment_pois_overdispersed_results_k20.csv", row.names = FALSE)

# Summary
summary_over_pois <- data.frame(Average_WAIC = average_WAIC, Average_LPML = average_LPML,
           SD_WAIC = sd_WAIC, SD_lpml = sd_LPML)

write.csv(summary_over_pois, file = "experiments/Poisson_mix/overdispersed_dataset/summary_over_pois_k20.csv", row.names = FALSE)
#---------------------------------------------------------------------------------#
# Repeating experiments 20 times Overdispersed Square Root Normal Mixture

# Results df for LPML and WAIC
results <- data.frame(dataset_id = integer(num_experiments), WAIC = numeric(num_experiments),
                       LPML = numeric(num_experiments))
for (i in 1:num_experiments) {
  # Fit the square root normal mixture model
  grid <- seq(min(datasets_overdispersed[[i]]) - 1, max(datasets_overdispersed[[i]]) + 1, len = 200)
  fit_norm <- gmm(y = sqrt(datasets_overdispersed[[i]]), grid = grid, K = K, amu = 0, b2mu = 0.01, 
                  asigma2 = 1, bsigma2 = 1, alpha = rep(0.1, K), nsim = 5000, nburn = 1000)
  
  # Extract posterior samples for weights, means, and variances
  samples_50 <- fit_norm$samples
  samples_matrix_50 <- as.matrix(samples_50)
  w_samples <- samples_matrix_50[, grep("w", colnames(samples_matrix_50))]
  mu_samples <- samples_matrix_50[, grep("mu", colnames(samples_matrix_50))]
  tau_samples <- samples_matrix_50[, grep("tau", colnames(samples_matrix_50))]
  sigma2_samples <- 1 / tau_samples
  
  # Calculate WAIC and LPML for each experiment
  ghc_norm <- good_fit_criteria(y = sqrt(datasets_overdispersed[[i]]), P = w_samples, Mu = mu_samples, 
                                Sigma2 = sigma2_samples, nburn = 1000)
  
  # Store results in a data frame
  results$dataset_id[i] <- i
  results$LPML[i] <- ghc_norm$LPML
  results$WAIC[i] <- ghc_norm$WAIC
}

# Calculate the average WAIC and LPML across all experiments
average_waic <- mean(results$WAIC)
average_lpml <- mean(results$LPML)
sd_lpml <- sd(results$LPML)
sd_waic <- sd(results$WAIC)

# Display the results for sqrt normal for overdispersed data
print(results)
write.csv(results, file = "experiments/Poisson_mix/overdispersed_dataset/experiment_sqrtNorm_overdispersed_results_k20.csv", row.names = FALSE)

# Summary
summary_over_sqrtnorm <- data.frame(Average_WAIC = average_waic, Average_LPML = average_lpml,
           SD_WAIC = sd_waic, SD_lpml = sd_lpml)

write.csv(summary_over_sqrtnorm, file = "experiments/Poisson_mix/overdispersed_dataset/summary_over_sqrtnorm_k20.csv", row.names = FALSE)

#------------------------------------------------------------------------------------#
# Repeating experiments for underdispersed data for poisson
set.seed(123)
num_datasets <- 20
# Initialize a list to store each dataset
datasets_underdispersed <- vector("list", num_datasets)

#Note underdispersed: Mean Data > Var Data
# Set Parameters
n <- 500
K <- 50
for (i in 1:num_datasets) {
  # Generate synthetic over-dispersed data (if you already have obs_data, skip this part)
  datasets_underdispersed[[i]] <- rbinom(500, 10, 0.8)
}
# Number of datasets to be iterated
num_experiments <- length(datasets_underdispersed)

# Results df for LPML and WAIC
results <- data.frame(dataset_id = integer(num_experiments), WAIC = numeric(num_experiments),
                       LPML = numeric(num_experiments))

# Loop to perform experiment 20 times for overdispersed

for (i in 1:num_experiments) {
  # alpha and beta prior maybe can use method of moment to determine the best results?
  # Fit Poisson mixture model with current dataset
  res_pois <- pmm(y = datasets_underdispersed[[i]], 
                  grid = seq(min(datasets_underdispersed[[i]]) - 1, 
                             max(datasets_underdispersed[[i]]) + 1, by = 1), 
                  K = K, alpha = 1, beta = 1, alpha_dir = rep(0.1, K))
  
  # Extract samples from the fitted model
  samples_mat_pois <- as.matrix(res_pois$samples)
  lambda_samples <- samples_mat_pois[, 1:K]
  w_samples <- samples_mat_pois[, (K + 1):(2 * K)]
  
  # Calculate good fit criteria
  gh3_pois <- good_fit_criteria_pois(datasets_underdispersed[[i]], w_samples, lambda_samples, nburn = 1000)
  
  # Store results in a data frame
  results$dataset_id[i] <- i
  results$LPML[i] <- gh3_pois$LPML
  results$WAIC[i] <- gh3_pois$WAIC
}

# Calculate average LPML and WAIC
average_LPML <- mean(results$LPML)
average_WAIC <- mean(results$WAIC)
sd_LPML <- sd(results$LPML)
sd_WAIC <- sd(results$WAIC)

# Display the results
print(results)
write.csv(results, file = "experiments/Poisson_mix/underdispersed_dataset/experiment_pois_underdispersed_results_k20.csv", row.names = FALSE)

# Summary
summary_over_pois <- data.frame(Average_WAIC = average_WAIC, Average_LPML = average_LPML,
           SD_WAIC = sd_WAIC, SD_lpml = sd_LPML)

write.csv(summary_over_pois, file = "experiments/Poisson_mix/underdispersed_dataset/summary_under_pois_k20.csv", row.names = FALSE)

#------------------------------------------------------------------------------------#
# Repeating experiments 20 times underdispersed Square Root Normal Mixture

# Results df for LPML and WAIC
results <- data.frame(dataset_id = integer(num_experiments), WAIC = numeric(num_experiments),
                       LPML = numeric(num_experiments))
for (i in 1:num_experiments) {
  # Fit the square root normal mixture model
  grid <- seq(min(datasets_underdispersed[[i]]) - 1, max(datasets_underdispersed[[i]]) + 1, len = 200)
  fit_norm <- gmm(y = sqrt(datasets_underdispersed[[i]]), grid = grid, K = K, amu = 0, b2mu = 0.01, 
                  asigma2 = 1, bsigma2 = 1, alpha = rep(0.1, K), nsim = 5000, nburn = 1000)
  
  # Extract posterior samples for weights, means, and variances
  samples_50 <- fit_norm$samples
  samples_matrix_50 <- as.matrix(samples_50)
  w_samples <- samples_matrix_50[, grep("w", colnames(samples_matrix_50))]
  mu_samples <- samples_matrix_50[, grep("mu", colnames(samples_matrix_50))]
  tau_samples <- samples_matrix_50[, grep("tau", colnames(samples_matrix_50))]
  sigma2_samples <- 1 / tau_samples
  
  # Calculate WAIC and LPML for each experiment
  ghc_norm <- good_fit_criteria(y = sqrt(datasets_underdispersed[[i]]), P = w_samples, Mu = mu_samples, 
                                Sigma2 = sigma2_samples, nburn = 1000)
  
  # Store results in a data frame
  results$dataset_id[i] <- i
  results$LPML[i] <- ghc_norm$LPML
  results$WAIC[i] <- ghc_norm$WAIC
}

# Calculate the average WAIC and LPML across all experiments
average_waic <- mean(results$WAIC)
average_lpml <- mean(results$LPML)
sd_lpml <- sd(results$LPML)
sd_waic <- sd(results$WAIC)

# Display the results for sqrt normal for overdispersed data
print(results)
write.csv(results, file = "experiments/Poisson_mix/underdispersed_dataset/experiment_sqrtNorm_underdispersed_results_k20.csv", row.names = FALSE)

# Summary
summary_over_sqrtnorm <- data.frame(Average_WAIC = average_waic, Average_LPML = average_lpml,
           SD_WAIC = sd_waic, SD_LPML = sd_lpml)

write.csv(summary_over_sqrtnorm, file = "experiments/Poisson_mix/underdispersed_dataset/summary_under_sqrtnorm_k20.csv", row.names = FALSE)

