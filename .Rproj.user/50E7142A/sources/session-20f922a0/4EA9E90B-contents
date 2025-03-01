source("code/Bayes_mix_normal.R")
require(mixtools)
require(coda)
require(rjags)
library(ggplot2)

# Repeating experiments 20 times Overfitted Mixture Normal

# Define the number of experiments
set.seed(123)
num_datasets <- 20
# Initialize a list to store each dataset
datasets <- vector("list", num_datasets)
n <- 500
for (i in 1:num_datasets) {
  # Generate the main normal data
  y <- rnormmix(n = n, lambda = c(0.6, rep(0.05, 8), 0.1), 
                               mu = seq(1, 20, by = 2), sigma = rep(1, 10))
  datasets[[i]] <- y
}

num_experiments <- length(datasets)  # Number of datasets (20 in this case)
K <- 20  # Number of mixture components, K = 20, 40, 50
# Initialize vectors to store WAIC and LPML values for each experiment
waic_values2 <- numeric(num_experiments)
lpml_values2 <- numeric(num_experiments)

# Loop through each dataset to fit the Gaussian mixture model and calculate good fit criteria
for (i in 1:num_experiments) {
  grid <- seq(min(datasets[[i]]) - 1, max(datasets[[i]]) + 1, len = 200)
  # Fit the Gaussian mixture model on the current dataset
  fit_norm <- gmm_nimble(y = datasets[[i]], grid = grid, K = K, amu = 5, b2mu = 0.010, 
                  asigma2 = 1, bsigma2 = 1, alpha_dirichlet = rep(0.1, K), nsim = 5000, nburn = 1000)
    # fit_norm <- gmm_nimble(y = datasets[[i]], grid = grid, K = K, amu = 5, b2mu = 0.010, 
    #               asigma2 = 1, bsigma2 = 1, alpha = rep(1, K), nsim = 5000, nburn = 1000)
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
sd_waic2 <- sd(waic_values2)
sd_lpml2 <- sd(lpml_values2)

# Combine WAIC and LPML values into a data frame
results_df <- data.frame(
  dataset_id = seq_len(num_experiments),
  WAIC       = waic_values2,
  LPML       = lpml_values2
)

# Inspect the results
print(results_df)

write.csv(results_df, file = "experiments/Overfitted_mix/alpha1/experiment_overfitted_results_k50.csv", row.names = FALSE)

# Print results
summary <- data.frame(Average_WAIC = average_waic2, Average_LPML = average_lpml2,
           SD_WAIC = sd_waic2, SD_lpml2 = sd_lpml2)

# Inspect the results
print(summary)

write.csv(summary, file = "experiments/Overfitted_mix/alpha1/experiment_overfitted_summary_k50.csv", row.names = FALSE)
