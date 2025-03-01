require(MCMCpack)
require(coda)
require(rjags)
require(nimble)
#----------------------------------------------------------------------------------------#
#' Poisson Mixture Model using Jags
#----------------------------------------------------------------------------------------#
#' This function implements a Poisson mixture model with K components using MCMC for inference
#'
#' @param y Observed data, a vector of length n
#' @param grid A set of grid points where the density is evaluated (min & max)
#' @param K Number of components in the mixture model
#' @param alpha Shape parameter for prior gamma distribution for each lambda_k
#' @param beta Rate parameter for prior gamma distribution for each lambda_k
#' @param alpha_dir Dirichlet prior parameters for the mixture weights w
#' @param n_iter Number of MCMC iterations for sampling
#' @param nburn Number of burn-in iterations for MCMC
#' 
#' @return A list containing the following:
#' \item{samples}{Posterior samples of the parameters from the MCMC}
#' \item{density}{Matrix of densities for each component at each grid point for all MCMC iterations}
#' \item{grid}{Grid points used for evaluating the densities}
pmm <- function(y, grid, K, alpha, beta, alpha_dir, n_iter = 5000, nburn =1000){
  
  # pmm hierarchy
  pois_mix_mod <- "model{
    for (i in 1:n){
      # latent var
      z[i] ~ dcat(w[])
      
      # likelihood
      y[i] ~ dpois(lambda[z[i]])
    }
    
    # priors
    for (k in 1:K){
      lambda[k] ~ dgamma(alpha,beta)
    }
    w ~ ddirch(alpha_dirichlet[])
  
  }"
  # --------------
  data_list <- list(
  y = y,
  n = length(y),
  K = K,
  alpha = alpha,
  beta = beta,
  alpha_dirichlet = alpha_dir
)
  inits <- function(){
  list(
    z = sample(1:K, size = length(y), replace = TRUE),
    lambda = rgamma(K, alpha, beta), 
    w = rep(1/K, K)
  )
  }
  params_pois <- c("lambda", "w")
  
  # model
  poisModel <- jags.model(textConnection(pois_mix_mod), data = data_list, inits = inits, n.chains = 1, n.adapt = 1000)
  
  # removes burn in
  update(poisModel, nburn)
  
  # sampling
  samples_pois <- coda.samples(poisModel, variable.names = params_pois, n.iter = n_iter)
  samples_mat_pois <- as.matrix(samples_pois)
  lambda_samples <- samples_mat_pois[, 1:K]
  w_samples <- samples_mat_pois[, (K+1):(2*K)]
  
  # ----------------
  # Compute Density
  # ----------------
  density <- matrix(0, nrow = nrow(w_samples), ncol = length(grid))
  
  for (g in 1:length(grid)) {
    for (k in 1:K) {
      density[, g] <- density[, g] + w_samples[, k] * dpois((grid[g]), lambda = lambda_samples[,k])
    }
  }
  return(list(samples = samples_pois, density = density, grid = grid))
}

#----------------------------------------------------------------------------------------#
#' Poisson Mixture Model with MCMC via Nimble
#----------------------------------------------------------------------------------------#
#' @param y Observed count data (a vector of integers)
#' @param grid Grid of values for evaluating densities (for plotting or other purposes)
#' @param K Number of mixture components
#' @param alpha Prior shape parameter for the Gamma prior on the Poisson rate (lambda)
#' @param beta Prior rate parameter for the Gamma prior on the Poisson rate (lambda)
#' @param alpha_dirichlet Prior concentration parameters for the Dirichlet distribution of the weights
#' @param nsim Number of MCMC iterations to run
#' @param nburn Number of burn-in iterations
#' @return A list containing:
#' \item{samples}{Posterior samples of the parameters from the MCMC}
#' \item{density}{Matrix of densities for each component at each grid point for all MCMC iterations}
#' \item{grid}{Grid points used for evaluating the densities}
pmm_nimble <- function(y, grid, K, alpha, beta, alpha_dir, nsim = 5000, nburn = 1000){
  
  # Poisson mixture model code for Nimble
  pois_mix_mod <- nimbleCode({
    for (i in 1:n){      
      z[i] ~ dcat(w[1:K])      
      y[i] ~ dpois(lambda[z[i]])
    }
    # Priors
    for (k in 1:K) {
      lambda[k] ~ dgamma(alpha, beta) 
    }
    
    # weights
    w[1:K] ~ ddirch(alpha_dirichlet[1:K]) 
    
  })
  # --------------
  data_list <- list(
    y = y,
    n = length(y),
    alpha = alpha,
    beta = beta
  )

  consts <- list(K = K, alpha_dirichlet = alpha_dir)

  inits <- function(){
    list(
      z = sample(1:K, size = length(y), replace = TRUE),
      lambda = rgamma(K, alpha, beta),
      w = rep(1/K, K)
    )
  }
  
  params_pois <- c("lambda", "w")
  
  # Config and compile NIMBLE
  n_model <- nimbleModel(code = pois_mix_mod, data = data_list, inits = inits(), constants = consts)
  mcmc_conf <- configureMCMC(n_model)
  mcmc <- buildMCMC(mcmc_conf)
  compiled_model <- compileNimble(n_model)
  compiled_mcmc <- compileNimble(mcmc, project = n_model)
  
  # sampling
  samples <- runMCMC(compiled_mcmc, niter = nsim)
  samples_matrix <- as.matrix(samples)
  lambda_samples <- samples_matrix[, grep("lambda", colnames(samples_matrix))]
  w_samples <- samples_matrix[, grep("w", colnames(samples_matrix))]
  
  # Initialize matrix of density
  density <- matrix(0, nrow = nrow(w_samples), ncol = length(grid))
  
  # ----------------
  # Compute Density
  # ----------------
  for (g in 1:length(grid)) {
    for (k in 1:K) {
      density[, g] <- density[, g] + w_samples[, k] * dpois(grid[g], lambda = lambda_samples[, k])
    }
  }
  
  return(list(samples = samples, density = density, grid = grid))
}

#----------------------------------------------------------------------------------------#
#' Dirichlet Process Poisson Mixture Model with MCMC via JAGS
#'
#' This function fits a Dirichlet Process Poisson Mixture Model (DP-PMM) using MCMC in JAGS.
#' The model assumes Poisson-distributed data for each component, with parameters sampled
#' from prior distributions, and a Dirichlet Process prior on the mixture weights using the
#' stick-breaking construction.
#'
#' @param y Numeric vector of observed count data.
#' @param grid Sequence of grid points over which density estimates are calculated.
#' @param L Truncation level for the number of mixture components in the model.
#' @param alambda Shape parameter for the prior Gamma distribution of the component rates `lambda_l`.
#' @param blambda Rate parameter for the prior Gamma distribution of the component rates `lambda_l`.
#' @param aalpha Shape parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param balpha Rate parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param n_iter Number of MCMC iterations for Gibbs sampling.
#' @param nburn Number of burn-in iterations to discard before computing results.
#'
#' @return A list containing:
#' \item{samples}{MCMC samples of the parameters from the posterior distribution.}
#' \item{density}{A matrix of density estimates for each grid point.}
#' \item{grid}{Grid points over which the density was computed.}

dp_pmm <- function(y, grid, L, alambda, blambda, aalpha, balpha, n_iter, nburn) {
  library(rjags)
  
  # Hierarchical structure
  dp_pmm_model <- "
  model{
    for(i in 1:n){
      # Latent Variable
      z[i] ~ dcat(w[])
      
      # Likelihood
      y[i] ~ dpois(lambda[z[i]])
    }
    
    # Priors for component rates
    for (l in 1:L){
      lambda[l] ~ dgamma(alambda, blambda)
    }
    
    # Stick-breaking construction for mixture weights
    for (l in 1:(L-1)) {
      v[l] ~ dbeta(1, alpha)
    }
    v[L] <- 1
    
    w[1] <- v[1]
    for (l in 2:L) {
      w[l] <- v[l] * (1-v[l-1]) * w[l-1] / v[l-1]
    }
    
    # Hyperprior on concentration parameter alpha
    alpha ~ dgamma(aalpha, balpha)
  }
  "
  data_list <- list(
    y = y,
    n = length(y),
    L = L,
    alambda = alambda,
    blambda = blambda,
    aalpha = aalpha,
    balpha = balpha
  )
  params_pois <- c("lambda", "w", "alpha", "z")
  
  # Compile the model
  dp_pois_model <- jags.model(textConnection(dp_pmm_model), data = data_list, 
                              n.chains = 1, n.adapt = 1000)
  
  # Burn-in period
  update(dp_pois_model, nburn)
  
  # Sampling
  samples_pois <- coda.samples(dp_pois_model, variable.names = params_pois, 
                               n.iter = n_iter)
  samples_matrix_dp_pmm <- as.matrix(samples_pois)
  w_samples <- samples_matrix_dp_pmm[, grep("^w\\[", colnames(samples_matrix_dp_pmm))]
  lambda_samples <- samples_matrix_dp_pmm[, grep("^lambda\\[", colnames(samples_matrix_dp_pmm))]
  
  # store density
  density <- matrix(0, nrow = nrow(w_samples), ncol = length(grid))
  
  # ----------------
  # Compute Density
  # ----------------
  for (g in 1:length(grid)) {
    for (l in 1:L) {
      density[, g] <- density[, g] + w_samples[, l] * dpois(grid[g], lambda = lambda_samples[, l])
    }
  }
  
  return(list(samples = samples_pois, density = density, grid = grid))
}

#----------------------------------------------------------------------------------------#
#' Dirichlet Process Poisson Mixture Model with MCMC via NIMBLE
#'
#' This function fits a Dirichlet Process Poisson Mixture Model (DP-PMM) using MCMC in NIMBLE.
#' The model assumes Poisson-distributed data for each component, with parameters sampled
#' from prior distributions, and a Dirichlet Process prior on the mixture weights using the
#' stick-breaking construction.
#'
#' @param y Numeric vector of observed count data.
#' @param grid Sequence of grid points over which density estimates are calculated.
#' @param L Truncation level for the number of mixture components in the model.
#' @param alambda Shape parameter for the prior Gamma distribution of the component rates `lambda_l`.
#' @param blambda Rate parameter for the prior Gamma distribution of the component rates `lambda_l`.
#' @param aalpha Shape parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param balpha Rate parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param n_iter Number of MCMC iterations.
#' @param nburn Number of burn-in iterations to discard before computing results.
#'
#' @return A list containing:
#' \item{samples}{MCMC samples of the parameters from the posterior distribution.}
#' \item{density}{A matrix of density estimates for each grid point.}
#' \item{grid}{Grid points over which the density was computed.}

dp_pmm_nimble <- function(y, grid, L, alambda, blambda, aalpha, balpha, n_iter, nburn) {
  
  # Hierarchy
  dp_pmm_code <- nimbleCode({
    for (i in 1:N) {
      # Latent variable
      z[i] ~ dcat(w[1:L])
      
      # Likelihood
      y[i] ~ dpois(lambda[z[i]])
    }
    
    for (l in 1:L) {
      lambda[l] ~ dgamma(alambda, blambda)
    }
    
    # Stick-breaking construction for mixture weights
    for (l in 1:(L - 1)) {
      v[l] ~ dbeta(1, alpha)
    }
    v[L] <- 1  # Set the last v[L] to 1 to complete the stick-breaking process
    
    # Use the stick_breaking function to construct weights w[]
    w[1:L] <- stick_breaking(v[1:(L-1)])
    
    # Prior for the concentration parameter alpha
    alpha ~ dgamma(aalpha, balpha)
  })
  
  data_list <- list(y = y)
  constants <- list(N = length(y), L = L, alambda = alambda, blambda = blambda,
                    aalpha = aalpha, balpha = balpha)
  
  # Create the Nimble model
  dp_pmm_model <- nimbleModel(code = dp_pmm_code, data = data_list, 
                              constants = constants)
  
  # Configure the MCMC sampler
  mcmc_conf <- configureMCMC(dp_pmm_model)
  
  # Use a random-walk Metropolis-Hastings sampler for the concentration parameter alpha
  mcmc_conf$removeSampler("alpha")
  mcmc_conf$addSampler(target = "alpha", type = "RW")
  mcmc_conf$addMonitors(c("lambda", "w", "alpha", "z"))
  
  # Build the MCMC object
  mcmc <- buildMCMC(mcmc_conf)
  compiled_model <- compileNimble(dp_pmm_model)
  compiled_mcmc <- compileNimble(mcmc, project = dp_pmm_model)
  
  # Sampling
  samples <- runMCMC(compiled_mcmc, niter = n_iter, nburnin = nburn)
  samples_mcmc <- as.mcmc(samples)
  w_samples <- samples[, grep("^w\\[", colnames(samples))]
  lambda_samples <- samples[, grep("^lambda\\[", colnames(samples))]
  
  # Compute densities over grid points
  density <- matrix(0, nrow = nrow(w_samples), ncol = length(grid))
  
  # ----------------
  # Compute Density
  # ----------------
  for (g in 1:length(grid)) {
    for (l in 1:L) {
      density[, g] <- density[, g] + w_samples[, l] * dpois(grid[g], lambda = lambda_samples[, l])
    }
  }
  
  return(list(samples = samples_mcmc, density = density, grid = grid))
}

#--------------------------------------------------------------------------#
#' Fitting Criteria for Poisson Mixture Model
#--------------------------------------------------------------------------#
#' This function computes WAIC and LPML for model performance metric
#'
#' @param y Numeric vector of observed data.
#' @param P Weights sample
#' @param Lambda Rates sample
#' @param nburn burn-in period (if needed)
#' 
#' @return a list containg:
#' \item{CPO}: CPO
#' \item{LPML}: LPML
#' \item{WAIC}: WAIC
good_fit_criteria_pois <- function(y, P, Lambda, nburn) {
  nsim <- nrow(P)  # number of MCMC iterations
  K <- ncol(Lambda)  # number of mixture components
  n <- length(y)  # number of observations

  Densym <- matrix(0, nrow = (nsim - nburn), ncol = n)

  # Loop over MCMC samples post-burn-in
  for (i in (nburn + 1):nsim) {
    for (k in 1:K) {
      Densym[i - nburn, ] <- Densym[i - nburn, ] + P[i, k] * dpois(y, lambda = Lambda[i, k])
    }
  }

  # Compute CPO and LPML
  cpoinv <- apply(1 / Densym, 2, mean)
  cpo <- 1 / cpoinv
  lpml <- sum(log(cpo))

  # Compute WAIC
  lpd <- sum(log(apply(Densym, 2, mean)))
  p2 <- sum(apply(log(Densym), 2, var))
  waic <- -2 * (lpd - p2)

  return(list(CPO = cpo, LPML = lpml, WAIC = waic))
}