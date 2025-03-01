require(MCMCpack)
require(coda)
require(rjags)
require(nimble)
#--------------------------------------------------------------------------#
#' Student-T Mixture Model with MCMC via NIMBLE
#--------------------------------------------------------------------------#
#'
#' This function fits a Student-T Mixture Model (TMM) with `K` components using 
#' Markov Chain Monte Carlo (MCMC) in NIMBLE.
#'
#' @param y Numeric vector of observed data.
#' @param grid Sequence of grid points over which density estimates are calculated.
#' @param K Number of mixture components in the model.
#' @param amu Mean of the normal prior for the component means `mu_k`.
#' @param b2mu Precision of the normal prior for the component means `mu_k`.
#' @param asigma2 Shape parameter for the prior Gamma distribution of the precisions `tau_k`.
#' @param bsigma2 Rate parameter for the prior Gamma distribution of the precisions `tau_k`.
#' @param rate_nu Rate parameter for the prior Exponential distribution on the degrees of freedom `nu_k`.
#' @param alpha_dirichlet Vector of parameters for the Dirichlet prior on the component weights.
#' @param nsim Number of MCMC iterations for Gibbs sampling.
#' @param nburn Number of burn-in iterations to discard before computing results.
#'
#' @return A list containing:
#' \item{samples}{MCMC samples of the parameters from the posterior distribution.}
#' \item{density}{A matrix of density estimates for each grid point.}
#' \item{grid}{Grid points over which the density was computed.}
tmix_nimble <- function(y, grid, K, amu, b2mu, asigma2, bsigma2, rate_nu, alpha_dirichlet, nsim, nburn) {
  
  # Hierarchy structure
  t_mixture_code <- nimbleCode({
    for (i in 1:N) {
      # Latent variable 
      z[i] ~ dcat(w[1:K])
      
      # Scaling factor omega for heavy tails
      omega[i] ~ dgamma(nu[z[i]] / 2, nu[z[i]] / 2)
      
      # Likelihood
      y[i] ~ dnorm(mu[z[i]], tau[z[i]] * omega[i])
    }
    
    for (k in 1:K) {
      # Priors for the component means and precisions (tau = 1/sigma^2)
      mu[k] ~ dnorm(amu, b2mu)
      tau[k] ~ dgamma(asigma2, bsigma2)
      
      # Prior for the degrees of freedom (using Exponential distribution as a simple prior)
      nu[k] ~ dexp(rate_nu)
    }
    
    # Prior for the mixture weights (Dirichlet distribution)
    w[1:K] ~ ddirch(alpha_dirichlet[1:K])
  })
  # --------------
  data_t_mixture <- list(y = y)
  constants <- list(N = length(y), K = K, amu = amu, b2mu = b2mu, 
                    asigma2 = asigma2, bsigma2 = bsigma2, rate_nu = rate_nu, alpha_dirichlet = alpha_dirichlet)

  # Initial values for the MCMC chain
  initsFunction <- function() {
    list(
      mu = rnorm(K, amu, sqrt(1/b2mu)),
      tau = rgamma(K, asigma2, bsigma2),
      nu = rexp(K, rate = rate_nu),
      z = sample(1:K, length(y), replace = TRUE),
      w = rep(1/K, K),
      omega = rgamma(length(y), 1, 1)  # Initialize omega for each observation
    )
  }
  
  # Configure and compile model
  t_mixture_model <- nimbleModel(t_mixture_code, data = data_t_mixture, constants = constants, inits = initsFunction())
  mcmc_conf <- configureMCMC(t_mixture_model)
  
  # Replace samplers for nu with RW Metropolis-Hastings samplers
  for (k in 1:K) {
    mcmc_conf$removeSampler(paste0("nu[", k, "]"))
    mcmc_conf$addSampler(target = paste0("nu[", k, "]"), type = "RW")
  }
  
  mcmc_conf$addMonitors(c("mu", "tau", "nu", "w", "omega"))
  mcmc <- buildMCMC(mcmc_conf)
  
  # Compile the model and the MCMC object for faster computation
  compiled_model <- compileNimble(t_mixture_model)
  compiled_mcmc <- compileNimble(mcmc, project = t_mixture_model)
  
  # Sampling
  samples <- runMCMC(compiled_mcmc, niter = nsim, nburnin = nburn)
  samples_mcmc <- as.mcmc(samples)
  w_samples <- samples[, grep("w", colnames(samples))]
  mu_samples <- samples[, grep("mu", colnames(samples))]
  tau_samples <- samples[, grep("tau", colnames(samples))]
  nu_samples <- samples[, grep("nu", colnames(samples))]
  omega_samples <- samples[, grep("omega", colnames(samples))]
  
  # ------------------
  # Compute Density
  # ------------------
  density <- matrix(0, nrow = nrow(mu_samples), ncol = length(grid))
  for (g in 1:length(grid)) {
    for (k in 1:K) {
      density[, g] <- density[, g] + w_samples[, k] * dt((grid[g] - mu_samples[, k]) / sqrt(1/tau_samples[, k]), df = nu_samples[, k]) / sqrt(1/tau_samples[, k])
    }
  }
  
  return(list(samples = samples_mcmc, density = density, grid = grid))
}

#----------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------#
#' DP Student-T Mixture Model with MCMC via NIMBLE
#--------------------------------------------------------------------------#
#'
#' This function fits a Student-T Mixture Model (TMM) with `K` components using 
#' Markov Chain Monte Carlo (MCMC) in NIMBLE.
#'
#' @param y Numeric vector of observed data.
#' @param grid Sequence of grid points over which density estimates are calculated.
#' @param L Truncation level for the number of mixture components in the model.
#' @param amu Mean of the normal prior for the component means `mu_k`.
#' @param b2mu Precision of the normal prior for the component means `mu_k`.
#' @param asigma2 Shape parameter for the prior Gamma distribution of the precisions `tau_k`.
#' @param bsigma2 Rate parameter for the prior Gamma distribution of the precisions `tau_k`.
#' @param rate_nu Rate parameter for the prior Exponential distribution on the degrees of freedom `nu_k`.
#' @param aalpha Shape parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param balpha Rate parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param n_iter Number of MCMC iterations for Gibbs sampling.
#' @param nburn Number of burn-in iterations to discard before computing results.
#'
#' @return A list containing:
#' \item{samples}{MCMC samples of the parameters from the posterior distribution.}
#' \item{density}{A matrix of density estimates for each grid point.}
#' \item{grid}{Grid points over which the density was computed.}
DP_tmix_nimble <- function(y, grid, L, amu, b2mu, asigma2, bsigma2, rate_nu, aalpha, balpha, n_iter, nburn) {
  require(nimble, warn.conflicts = FALSE)
  
  # Hierarchy structure
  dp_t_mixture_code <- nimbleCode({
    for (i in 1:N) {
      # Latent variable
      z[i] ~ dcat(w[1:L])
      
      # Scaling factor omega for heavy tails
      omega[i] ~ dgamma(nu[z[i]] / 2, nu[z[i]] / 2)
      
      # Likelihood
      y[i] ~ dnorm(mu[z[i]], tau[z[i]] * omega[i])
    }
    
    for (l in 1:L) {
      mu[l] ~ dnorm(amu, b2mu)
      tau[l] ~ dgamma(asigma2, bsigma2)
      nu[l] ~ dexp(rate_nu)
    }
    
    # Stick-breaking construction for the mixture weights
    for (l in 1:(L - 1)) {
      v[l] ~ dbeta(1, alpha)
    }
    v[L] <- 1  # Ensure the last v[L] is set to 1
    
    # Construct the mixture weights w[] using stick-breaking
    w[1:L] <- stick_breaking(v[1:(L-1)])
    
    # Prior for the concentration parameter alpha
    alpha ~ dgamma(aalpha, balpha)
  })
  # --------------
  data_list <- list(y = y)
  constants <- list(N = length(y), L = L, amu = amu, b2mu = b2mu,
                    asigma2 = asigma2, bsigma2 = bsigma2,
                    rate_nu = rate_nu, aalpha = aalpha, balpha = balpha)
  
  # Configure and compile the model
  dp_t_mixture_model <- nimbleModel(code = dp_t_mixture_code, data = data_list, constants = constants)
  mcmc_conf <- configureMCMC(dp_t_mixture_model)
  
  # Remove the existing slice sampler for `nu` and `v`
  for (l in 1:L) {
    mcmc_conf$removeSampler(paste0("nu[", l, "]"))
  }
  for (l in 1:(L - 1)) {
    mcmc_conf$removeSampler(paste0("v[", l, "]"))
  }
  mcmc_conf$removeSampler("alpha")
  
  # Add Metropolis-Hastings samplers for `nu`, `v`, and `alpha`
  for (l in 1:L) {
    mcmc_conf$addSampler(target = paste0("nu[", l, "]"), type = "RW", control = list(scale = 0.5))  # Adjust scale as needed
  }
  for (l in 1:(L - 1)) {
    mcmc_conf$addSampler(target = paste0("v[", l, "]"), type = "RW", control = list(scale = 0.5))  # Adjust scale as needed
  }
  mcmc_conf$addSampler(target = "alpha", type = "RW", control = list(scale = 0.5))  # Adjust scale as needed
  mcmc_conf$addMonitors(c("mu", "tau", "nu", "w", "alpha", "z", "omega"))

  mcmc <- buildMCMC(mcmc_conf, useConjugacy = FALSE)
  compiled_model <- compileNimble(dp_t_mixture_model)
  compiled_mcmc <- compileNimble(mcmc, project = dp_t_mixture_model)
  
  # Sampling
  samples <- runMCMC(compiled_mcmc, niter = n_iter, nburnin = nburn)
  samples_mcmc <- as.mcmc(samples)
  w_samples <- samples[, grep("^w\\[", colnames(samples))]
  mu_samples <- samples[, grep("^mu\\[", colnames(samples))]
  tau_samples <- samples[, grep("^tau\\[", colnames(samples))]
  nu_samples <- samples[, grep("^nu\\[", colnames(samples))]
  omega_samples <- samples[, grep("^omega\\[", colnames(samples))]
  
  # ------------------
  # Compute Density
  # ------------------
  density <- matrix(0, nrow = nrow(mu_samples), ncol = length(grid))
  for (g in 1:length(grid)) {
    for (l in 1:L) {
      standardized_x <- (grid[g] - mu_samples[, l]) * sqrt(tau_samples[, l])
      density[, g] <- density[, g] + w_samples[, l] * dt(standardized_x, df = nu_samples[, l]) * sqrt(tau_samples[, l])
    }
  }
  
  return(list(samples = samples_mcmc, density = density, grid = grid))
}

#----------------------------------------------------------------------------------------#
#' Fitting Criteria for Student-t Mixture Model
#--------------------------------------------------------------------------#
#' This function computes WAIC and LPML for model performance metric
#'
#' @param y Numeric vector of observed data.
#' @param P Weights sample
#' @param Mu Mean samples
#' @param Sigma2 Variance samples
#' @param df_samples degree of freedom samples
#' @param nburn burn-in period (if needed)
#' 
#' @return a list containg:
#' \item{CPO}: CPO
#' \item{LPML}: LPML
#' \item{WAIC}: WAIC

good_fit_criteria_t <- function(y, P, Mu, Sigma2, df_samples, nburn) {
  nsim <- nrow(P)   # number of MCMC iterations
  K <- ncol(Mu)     # number of mixture components
  n <- length(y)    # number of observations
  
  Densy <- array(0, c((nsim - nburn), n, K))            # density of each component for each observation
  Densym <- matrix(0, nrow = (nsim - nburn), ncol = n)  # total mixture density for each observation
  
  # Loop over MCMC samples post-burn-in
  for (i in (nburn + 1):nsim) {
    for (k in 1:K) {
      Densy[i - nburn, , k] <- P[i, k] * dt((y - Mu[i, k]) / sqrt(Sigma2[i, k]), df = df_samples[i, k]) / sqrt(Sigma2[i, k])
    }
    for (j in 1:n) {
      Densym[i - nburn, j] <- sum(Densy[i - nburn, j, ])
    }
  }
  
  # Compute LPML
  cpoinv <- apply(1 / Densym, 2, mean)
  cpo <- 1 / cpoinv
  lpml <- sum(log(cpo))
  
  # Compute WAIC
  lpd <- sum(log(apply(Densym, 2, mean)))
  p2 <- sum(apply(log(Densym), 2, var))
  waic <- -2 * (lpd - p2)
  
  return(list(CPO = cpo, LPML = lpml, WAIC = waic))
}