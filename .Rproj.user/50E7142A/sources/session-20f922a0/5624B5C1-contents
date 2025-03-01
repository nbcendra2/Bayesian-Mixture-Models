require(MCMCpack)
require(coda)
require(rjags)
require(nimble)
#--------------------------------------------------------------------------#
#' Gaussian Mixture Model with MCMC via JAGS
#--------------------------------------------------------------------------#
#'
#' This function fits a Gaussian Mixture Model (GMM) with `K` components using MCMC in JAGS.
#'
#' @param y Numeric vector of observed data.
#' @param grid Sequence of grid points over which density estimates are calculated.
#' @param K Number of mixture components in the model.
#' @param amu Mean of the normal prior for the component means `mu_k`.
#' @param b2mu Precision of the normal prior for the component means `mu_k`.
#' @param asigma2 Shape parameter for the prior Gamma distribution of the variances `sigma2_k`.
#' @param bsigma2 Rate parameter for the prior Gamma distribution of the variances `sigma2_k`.
#' @param alpha_dirichlet Vector of parameters for the Dirichlet prior on the component weights.
#' @param nsim Number of MCMC iterations for Gibbs sampling.
#' @param nburn Number of burn-in iterations to discard before computing results.
#' 
#' @return A list containing:
#' \item{samples}{MCMC samples of the parameters from the posterior distribution.}
#' \item{density}{A matrix of density estimates for each grid point.}
#' \item{grid}{Grid points over which the density was computed.}
gmm <- function(y, grid, K, amu, b2mu, asigma2, bsigma2, alpha_dirichlet, nsim, nburn) {
  
  # Hierarchical model structure
  gmm_model <- "
  model {
    for (i in 1:n) {
      # Latent variable
      z[i] ~ dcat(w[]) 
      
      # Likelihood
      y[i] ~ dnorm(mu[z[i]], tau[z[i]])  
    }
  
    for (k in 1:K) {
      # Priors for the means and precisions of the components
      mu[k] ~ dnorm(amu, b2mu)
      tau[k]  ~ dgamma(asigma2, bsigma2)
    }
    
    # Components weights
    w[1:K] ~ ddirch(alpha_dirichlet[])
  }
  "
  # --------------
  data_gmm <- list(y = y, n = length(y), K = K, amu = amu, b2mu = b2mu, 
                   asigma2 = asigma2, bsigma2 = bsigma2, alpha_dirichlet = alpha_dirichlet)
  
  inits_gmm <- function() {
    list(
      mu = seq(-5,5, length.out = K),
      tau = rep(1, K),
      z = sample(1:K, length(y), replace = TRUE)
    )
  }
  
  # Parameters to monitor (mu, tau, and w)
  params_gmm <- c("mu", "tau", "w")
  
  jagsModel_gmm <- jags.model(textConnection(gmm_model), data = data_gmm, inits = inits_gmm, 
                              n.chains = 1, n.adapt = 1000)
  
  # Burn-in and sampling
  update(jagsModel_gmm, nburn)
  samples_gmm <- coda.samples(jagsModel_gmm, variable.names = params_gmm, n.iter = nsim)
  samples_matrix_gmm <- as.matrix(samples_gmm)
  
  # Extract posterior samples for weights, means, and variances
  w_samples_gmm <- samples_matrix_gmm[, grep("w", colnames(samples_matrix_gmm))]
  mu_samples_gmm <- samples_matrix_gmm[, grep("mu", colnames(samples_matrix_gmm))]
  tau_samples_gmm <- samples_matrix_gmm[, grep("tau", colnames(samples_matrix_gmm))]
  
  # ------------
  # Evaluate Density
  # ------------
  density <- matrix(0, nrow = nrow(w_samples_gmm), ncol = length(grid))
  for (g in 1:length(grid)) {
    for (k in 1:K) {
      density[, g] <- density[, g] + w_samples_gmm[, k] * dnorm(grid[g], 
                                                                mean = mu_samples_gmm[, k], 
                                                                sd = sqrt(1/tau_samples_gmm[, k]))
    }
  }
  return(list(samples = samples_gmm, density = density, grid = grid))
}

#--------------------------------------------------------------------------#
#' Gaussian Mixture Model with MCMC via NIMBLE
#--------------------------------------------------------------------------#
#' 
#' @param y Observed data, a vector of length n
#' @param grid A sequence of points over which to evaluate the density
#' @param K Number of components in the mixture model
#' @param amu Prior mean for the normal distribution of the component means (mu)
#' @param b2mu Precision (1/variance) for the prior on the component means
#' @param asigma2 Shape parameter for the gamma prior on the component variances (sigma^2)
#' @param bsigma2 Rate parameter for the gamma prior on the component variances (sigma^2)
#' @param alpha_dirichlet Parameters of the Dirichlet prior on the mixture weights (w)
#' @param nsim Number of MCMC iterations to run
#' @param nburn Number of MCMC burn-in iterations to discard
#' @return A list containing:
#' \item{samples}{MCMC samples of the parameters from the posterior distribution.}
#' \item{density}{A matrix of density estimates for each grid point.}
#' \item{grid}{Grid points over which the density was computed.}
gmm_nimble <- function(y, grid, K, amu, b2mu, asigma2, bsigma2, alpha_dirichlet, nsim, nburn){
  # Hiearchical model structure
  norm_mix <- nimbleCode({
    for(i in 1:n){
      # Latent variable
      z[i] ~ dcat(w[1:K])
      
      # Likelihood
      y[i] ~ dnorm(mu[z[i]], tau[z[i]])
    }
    
    # Loop over each component
    for(k in 1:K){ 
      # Prior for the component means mu and tau
      mu[k] ~ dnorm(amu, b2mu)
      tau[k] ~ dgamma(asigma2, bsigma2)
    }
    w[1:K] ~ ddirch(alpha_dirichlet[1:K])
  })  
  
  data_mix <- list(y = y)
  constants <- list(n = length(y), K = K, amu = amu, b2mu = b2mu, 
                        asigma2 = asigma2, bsigma2 = bsigma2, 
                        alpha_dirichlet = alpha_dirichlet)
  
  # Initial values for the MCMC chain
  initsFunction <- function() {
    list(
      mu = seq(-5, 5, length.out = K),
      tau = rep(1, K),
      z = sample(1:K, length(y), replace = TRUE)
    )
  }
  
  # -------------
  # Compile and run NIMBLE
  # -------------
  n_model <- nimbleModel(norm_mix, data = data_mix, constants = constants, inits = initsFunction())
  mcmc_conf <- configureMCMC(n_model)
  mcmc_conf$addMonitors("tau")
  mcmc <- buildMCMC(mcmc_conf)
  compiled_model <- compileNimble(n_model)
  compiled_mcmc <- compileNimble(mcmc, project = n_model)
  
  # Sample from model
  samples <- runMCMC(compiled_mcmc, niter = nsim)
  samples_mcmc <- as.mcmc(samples)
  
  # Extract the posterior samples for specific components
  w_samples <- samples[, grep("w", colnames(samples))]      # Mixture weights
  mu_samples <- samples[, grep("mu", colnames(samples))]    # Component means
  tau_samples <- samples[, grep("tau", colnames(samples))]  # Component precision

  n_samples <- nrow(mu_samples)  # Number of MCMC samples
  comp <- ncol(mu_samples)  # Number of mixture components
  density <- matrix(0, nrow = n_samples, ncol = length(grid))
  
  # ------------
  # Evaluate Density
  # ------------
  for (g in 1:length(grid)) {
    for (k in 1:comp) {
      # Add the weighted density for each component at each grid point
      density[, g] <- density[, g] + w_samples[, k] * dnorm(grid[g], mean = mu_samples[, k], sd = sqrt(1/tau_samples[, k]))
    }
  }

  
  return(list(samples = samples_mcmc, density = density, grid = grid))
}


#--------------------------------------------------------------------------#
#'Dirichlet Process Gaussian Mixture Model with MCMC via JAGS
#--------------------------------------------------------------------------#
#' This function fits a Dirichlet Process Gaussian Mixture Model (DP-GMM) using MCMC in JAGS.
#' Blocked Gibbs sampler is used
#'
#' @param y Numeric vector of observed data.
#' @param grid Sequence of grid points over which density estimates are calculated.
#' @param L Truncation level for the number of mixture components in the model.
#' @param amu Mean of the normal prior for the component means `mu_l`.
#' @param b2mu Precision of the normal prior for the component means `mu_l`.
#' @param asigma2 Shape parameter for the prior Gamma distribution of the precisions `tau_l`.
#' @param bsigma2 Rate parameter for the prior Gamma distribution of the precisions `tau_l`.
#' @param aalpha Shape parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param balpha Rate parameter for the prior Gamma distribution of the concentration parameter `alpha`.
#' @param n_iter Number of MCMC iterations for Gibbs sampling.
#' @param nburn Number of burn-in iterations to discard before computing results.
#' 
#' @return A list containing:
#' \item{samples}{MCMC samples of the parameters from the posterior distribution.}
#' \item{density}{A matrix of density estimates for each grid point.}
#' \item{grid}{Grid points over which the density was computed.}
dp_gmm <- function(y, grid, L, amu, b2mu, asigma2, bsigma2, aalpha, balpha, n_iter, nburn) {
  # Hierarchical structure
  dp_gmm_model <- "
  model {
    for (i in 1:n) {
      # Latent variable 
      z[i] ~ dcat(w[])
      # Likelihood
      y[i] ~ dnorm(mu[z[i]], tau[z[i]])
    }

    for (l in 1:L) {
      # Priors for the means and precisions of the components
      mu[l] ~ dnorm(amu, b2mu)
      tau[l] ~ dgamma(asigma2, bsigma2)
    }

    # Stick-breaking construction for mixture weights
    for (l in 1:(L-1)) {
      v[l] ~ dbeta(1, alpha)
    }
    
    v[L] <- 1  # Ensure the last v[L] is set to 1

    # Construct the mixture weights w[]
    w[1] <- v[1]
    for (l in 2:L) {
      w[l] <- v[l] * (1 - v[l-1]) * w[l-1] / v[l-1]
    }

    # Prior for the concentration parameter alpha
    alpha ~ dgamma(aalpha, balpha)
  }
  "
  # --------------
  data_dp_gmm <- list(y = y, n = length(y), L = L, amu = amu, b2mu = b2mu, 
                      asigma2 = asigma2, bsigma2 = bsigma2, aalpha = aalpha, balpha = balpha)
  
  # Parameters to monitor (mu, tau, w, alpha)
  params_dp_gmm <- c("mu", "tau", "w", "alpha", "z")
  
  
  # Run model
  jagsModel_dp_gmm <- jags.model(textConnection(dp_gmm_model), data = data_dp_gmm 
                                 , n.chains = 1, n.adapt = 1000)
  # Burn-in and sampling
  update(jagsModel_dp_gmm, nburn)
  samples_dp_gmm <- coda.samples(jagsModel_dp_gmm, variable.names = params_dp_gmm, 
                                 n.iter = n_iter)
  samples_matrix_dp_gmm <- as.matrix(samples_dp_gmm)
  
  # Extract posterior samples for weights, means, precisions, and alpha
  w_samples_dp_gmm <- samples_matrix_dp_gmm[, grep("^w\\[", colnames(samples_matrix_dp_gmm))]
  mu_samples_dp_gmm <- samples_matrix_dp_gmm[, grep("^mu\\[", colnames(samples_matrix_dp_gmm))]
  tau_samples_dp_gmm <- samples_matrix_dp_gmm[, grep("^tau\\[", colnames(samples_matrix_dp_gmm))]
  alpha_samples_dp_gmm <- samples_matrix_dp_gmm[, "alpha"]
  
  # ----------------
  # Compute Density
  # ----------------
  density <- matrix(0, nrow = nrow(w_samples_dp_gmm), ncol = length(grid))
  # Loop over grid points to compute density
  for (g in 1:length(grid)) {
    for (l in 1:L) {
      density[, g] <- density[, g] + w_samples_dp_gmm[, l] * dnorm(grid[g], 
                                                             mean = mu_samples_dp_gmm[, l], 
                                                             sd = sqrt(1 / tau_samples_dp_gmm[, l]))
    }
  }
  
  return(list(samples = samples_dp_gmm, density = density, grid = grid))
}

#--------------------------------------------------------------------------#
#' Fitting Criteria for Gaussian Mixture Model
#--------------------------------------------------------------------------#
#' This function computes WAIC and LPML for model performance metric
#'
#' @param y Numeric vector of observed data.
#' @param P Weights sample
#' @param Mu Mean samples
#' @param Sigma2 Variance samples
#' @param nburn burn-in period (if needed)
#' 
#' @return a list containg:
#' \item{CPO}: CPO
#' \item{LPML}: LPML
#' \item{WAIC}: WAIC
good_fit_criteria <- function(y, P, Mu, Sigma2, nburn) {
  nsim <- nrow(P)
  K <- ncol(Mu)
  n <- length(y)
  
  # store the density of each component k for each observation in y at each iteration after nburn
  Densy <- array(0, c((nsim - nburn), n, K))
  
  # stores total mixture density for each observation at each mcmc iteration
  Densym <- matrix(0, nrow = (nsim - nburn), ncol = n)
  
  for (i in (nburn + 1):nsim) {
    for (k in 1:K) {
      Densy[i - nburn, , k] <- P[i, k] * dnorm(y, Mu[i,
                                                     k], sqrt(Sigma2[i, k]))
    }
    for (j in 1:n) {
      Densym[i - nburn, j] <- sum(Densy[i - nburn, j, ])
    }
  }
  
  cpoinv <- apply(1/Densym, 2, mean)
  cpo <- 1/cpoinv
  lpml <- sum(log(cpo))
  
  lpd <- sum(log(apply(exp(log(Densym)), 2, mean)))
  p2 <- sum(apply(log(Densym), 2, var))
  waic <- -2 * (lpd - p2)
  
  return(list(CPO = cpo, LPML = lpml, WAIC = waic))
}
