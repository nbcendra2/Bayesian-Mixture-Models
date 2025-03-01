require(mixtools)
require(coda)
require(rjags)
source("code/Bayes_mix_normal.R")
source("code/Bayes_mix_poisson.R")
source("code/Bayes_mix_student_t.R")

# Repeating experiments 20 times Overfitted Mixture Normal

# Define the number of experiments
set.seed(123)
num_datasets <- 20
n <- 500
# Initialize a list to store each dataset
datasets <- vector("list", num_datasets)
for (i in 1:num_datasets) {
  # Generate the main normal data
  y <- rnormmix(n = n, lambda = c(0.6, rep(0.05, 8), 0.1), 
                               mu = seq(1, 20, by = 2), sigma = rep(1, 10))
  datasets[[i]] <- y
}

# Number of datasets (20 in this case)
num_experiments <- length(datasets)

# Upperbound of mixture components
L <- 50  

# Initialize vectors to store WAIC and LPML values for each experiment
waic_values2 <- numeric(num_experiments)
lpml_values2 <- numeric(num_experiments)

# Loop through each dataset to fit the Gaussian mixture model and calculate good fit criteria
for (i in 1:num_experiments) {
  grid <- seq(min(datasets[[i]]) - 1, max(datasets[[i]]) + 1, len = 200)
  # Fit the Gaussian mixture model on the current dataset
  fit_norm <- dp_gmm(y = datasets[[i]],grid = grid,L = L, amu =  5, 
                     b2mu = 0.01,asigma2 = 1,bsigma2 = 1, 
                     aalpha = 1,balpha =  1,n_iter = 5000,nburn = 1000)
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

write.csv(results_df, file = "experiments/Overfitted_mix/dp/experiment_overfitted_dp_results.csv", row.names = FALSE)

# Print results
summary <- data.frame(Average_WAIC = average_waic2, Average_LPML = average_lpml2,
           SD_WAIC = sd_waic2, SD_lpml2 = sd_lpml2)

# Inspect the results
print(summary)

write.csv(summary, file = "experiments/Overfitted_mix/dp/summary_overfitted_dp.csv", row.names = FALSE)

#---------------------------------------------------------------------------------#
# Repeat experiment 20 times for DP Pois mixture
set.seed(123)
num_datasets2 <- 20    # Number of datasets to generate
n <- 500
# Initialize a list to store each dataset
datasets_pois <- vector("list", num_datasets2)

# Generate 20 datasets with outliers and store them in the list
for (l in 1:num_datasets2) {
  y_van <- numeric(n)
  u <- runif(n)
  for(i in 1:n){
    if(u[i] < 0.4) {y_van[i] <- rpois(1, lambda = 2)} 
    else
      if(u[i] < 0.8) {y_van[i] <- rpois(1, lambda = 6)}
        else {y_van[i] <- rpois(1, lambda = 10)}
  }
  datasets_pois[[l]] <- y_van
}

num_experiments_pois <- length(datasets_pois)  # Number of datasets (20 in this case)
L <- 50  # Number of mixture components
# Initialize vectors to store WAIC and LPML values for each experiment
waic_values_pois <- numeric(num_experiments_pois)
lpml_values_pois <- numeric(num_experiments_pois)

results_pois <- data.frame(dataset_id = integer(num_experiments_pois), WAIC = numeric(num_experiments_pois),
                       LPML = numeric(num_experiments_pois))

# Loop through each dataset to fit the Gaussian mixture model and calculate good fit criteria
for (i in 1:num_experiments_pois) {
  grid_pois <- seq(min(datasets_pois[[i]]) - 1, max(datasets_pois[[i]]) + 1, by = 1)
  # Fit the Gaussian mixture model on the current dataset
  fit_dp_pmm <- dp_pmm_nimble(y = datasets_pois[[i]], grid = grid_pois, L = L, 
                      alambda = 1, blambda = 1, aalpha = 1, 
                      balpha =1, n_iter = 5000, nburn = 1000)
  # Extract posterior samples for weights, means, and variances
  samples_pois <- fit_dp_pmm$samples
  samples_matrix_pois <- as.matrix(samples_pois)
  
  # Extract weights, lambda
  w_samples_pois <- samples_matrix_pois[, grep("w", colnames(samples_matrix_pois))]
  lambda_samples <- samples_matrix_pois[, grep("lambda", colnames(samples_matrix_pois))]
  
  # Calculate WAIC and LPML for the current experiment
  ghc_pois <- good_fit_criteria_pois(y = datasets_pois[[i]], P = w_samples_pois, Lambda = lambda_samples, 
                                nburn = 1000)
  
  # Store the WAIC and LPML values
  results_pois$dataset_id[i] <- i
  results_pois$WAIC[i] <- ghc_pois$WAIC
  results_pois$LPML[i] <- ghc_pois$LPML
}

# Calculate the average WAIC and LPML across all experiments
average_waic_pois <- mean(results_pois$WAIC)
average_lpml_pois <- mean(results_pois$LPML)
sd_LPML_pois <- sd(results_pois$LPML)
sd_WAIC_pois <- sd(results_pois$WAIC)

# Display the results
print(results_pois)
write.csv(results_pois, file = "experiments/Poisson_mix/dp/experiment_pois_overdispersed_dp_results.csv", row.names = FALSE)

# Summary
summary_over_pois <- data.frame(Average_WAIC = average_waic_pois, Average_LPML = average_lpml_pois,
           SD_WAIC = sd_WAIC_pois, SD_lpml = sd_LPML_pois)

write.csv(summary_over_pois, file = "experiments/Poisson_mix/dp/summary_over_pois_dp.csv", row.names = FALSE)

#---------------------------------------------------------------------------------#
# Repeat experiment 20 times for DP t mix
set.seed(123)
# Parameters for the main normal distribution
N <- 100
mu <- 3
sigma <- 1.15
# For DP use 10% outlier
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

# looping 20 times.
# Set parameters
num_experiments <- length(datasets)  # Number of datasets (20 in this case)
L <- 50  # Number of mixture components

# head(datasets)
# Initialize vectors to store WAIC and LPML values for each experiment
waic_values <- numeric(num_experiments)
lpml_values <- numeric(num_experiments)

# Loop through each dataset to fit the t-mixture model and calculate good fit criteria
for (i in 1:num_experiments) {
  grid <- seq(min(datasets[[i]]) - 1, max(datasets[[i]]) + 1, len = 200) 
  # Fit the t-mixture model on the current dataset
  tmix_sample <- DP_tmix_nimble(datasets[[i]], grid = grid, L = L, amu = 0, b2mu = 0.01,
                           asigma2 = 1, bsigma2 = 1, rate_nu = 0.1, aalpha = 1, balpha = 1,
                           n_iter = 5000, nburn = 1000)
  
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

# Print results
# Df results
results <- data.frame(dataset_id = 1:num_experiments, WAIC = waic_values, LPML = lpml_values)
write.csv(results, file = "experiments/t_mix/dp/experiment_heavy_t_mix_dp_results.csv", row.names = FALSE)

# Df summary
summary <- data.frame(Average_WAIC = average_waic, Average_LPML = average_lpml,
                      SD_WAIC = sd_waic, SD_LPML <- sd_lpml)
write.csv(summary, file = "experiments/t_mix/dp/summary_heavy_t_mix_dp.csv", row.names = FALSE)