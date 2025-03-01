#--------------------------------------------------------------------------#
#' Covariate-Dependent Variance for the Linear Dirichlet Dependent Mixture Model
#--------------------------------------------------------------------------#
#'
#' This function fits a Linear Dirichlet Dependent Process (LDDP) model 
#' with covariate-dependent mean and variance using Markov Chain Monte Carlo (MCMC) in JAGS.
#'
#' @param y Numeric vector of observed response values.
#' @param X Numeric matrix of predictor variables (design matrix).
#' @param prior A list specifying prior hyperparameters:
#'   \itemize{
#'     \item \code{L}: Maximum number of mixture components.
#'     \item \code{m0}: Mean vector for the normal prior on regression coefficients \code{beta}.
#'     \item \code{S0}: Covariance matrix for the normal prior on \code{beta}.
#'     \item \code{a}, \code{b}: Mean vector (\code{a}) and inverse covariance matrix (\code{b}) 
#'           for the normal prior on variance regression coefficients \code{gamma}.
#'     \item \code{Psi}: Scale matrix for the Wishart prior on covariance.
#'     \item \code{nu}: Degrees of freedom for the Wishart prior.
#'     \item \code{aalpha}, \code{balpha}: Shape and rate parameters for the Gamma prior on concentration parameter \code{alpha}.
#'   }
#' @param nburn Number of burn-in iterations to discard before computing results.
#' @param niter Number of MCMC iterations.
#' @param thin Thinning factor for MCMC samples.
#' @param n_chains Number of MCMC chains.
#' @param adapt_steps Number of adaptation steps for JAGS.
#' @param seed Random seed for reproducibility.
#'
#' @return A list containing:
#' \item{mcmc_samples}{MCMC samples from the posterior distribution as a `coda` object.}
#' \item{w_samples}{Matrix of posterior samples for the mixture weights.}
#' \item{beta_samples}{3D array of posterior samples for regression coefficients \code{beta} (iterations x L x P).}
#' \item{gamma_samples}{3D array of posterior samples for variance coefficients \code{gamma} (iterations x L x P).}
#' \item{alpha_samples}{Matrix of posterior samples for the concentration parameter \code{alpha}.}
#' \item{z_samples}{Matrix of posterior samples for latent component assignments \code{z}.}
lddp_var <- function(y, X, prior, nburn = 1000, niter = 5000, thin = 1, n_chains = 1, 
                 adapt_steps = 1000, seed = 123) {
  # Library
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("The 'rjags' package is required but not installed.")
  }
  require(rjags)
  # model
  lddp_model <- "
    model {
      for (i in 1:n){
        z[i] ~ dcat(w[])
        y[i] ~ dnorm(inprod(x[i, ], beta[z[i], ]), exp(-inprod(x[i, ], gamma[z[i], ])))
      }
      for (l in 1:L){
        beta[l, 1:P] ~ dmnorm(abeta[], Sigmainv[, ])
        gamma[l, 1:P] ~ dmnorm(agamma[], bgamma[,])
      }
      abeta[1:P] ~ dmnorm(amu[], b2mu[,])
      # for (i in 1:P) {
      #     psi_adj[i, i] <- psi[i, i] + 1e-6
      # }
      # # ensures positive definite matrix, adding a regularization term.
      # for (i in 1:P) {
      #     for (j in (i+1):P) {
      #         psi_adj[i, j] <- psi[i, j]
      #         psi_adj[j, i] <- psi[j, i]
      #     }
      # }
      # Sigmainv[1:P, 1:P] ~ dwish(nu * psi_adj[, ], nu)
      Sigmainv[1:P, 1:P] ~ dwish(nu * psi[, ], nu)
      
      # Stick breaking construction
      # T(0.001,0.999) to avoid computational error
      # jags slice sampler doesnt work when prob dist of sampled var is infinite at some point
      # This happend when beta dist at 0 or 1
      for (l in 1:(L-1)) {
        v[l] ~ dbeta(1, alpha) T(0.001,0.999)
      }
      v[L] <- 1
      w[1] <- v[1]
      for (l in 2:L) {
        w[l] <- v[l] * (1 - v[l-1]) * w[l-1] / v[l-1]
      }
      alpha ~ dgamma(aalpha, balpha)
    }
  "
  
  # Prepare data for JAGS
  data_list <- list(
    n = length(y),           # Obs
    L = prior$L,             # Num of upper bound comp
    P = ncol(X),             # Num of predictors
    y = y,                   # Response variable
    x = X,                   # Desgn matrix
    agamma = prior$a,       
    bgamma = prior$b,       
    amu = prior$m0,          # Mean vector for beta coefficients
    b2mu = solve(prior$S0),         # Precision matrix for beta coefficients -> perlu di inverse? 
    psi = prior$Psi,         # Scale matrix for Wishart prior
    nu = prior$nu,           # Degrees of freedom for Wishart prior
    aalpha = prior$aalpha,   # Shape parameter for Gamma prior on alpha
    balpha = prior$balpha    # Rate parameter for Gamma prior on alpha
  )
  
  # MCMC configuration
  mcmc_config <- list(
    nburn = nburn,
    niter = niter,
    thin = thin
  )
  L <- prior$L
  P <- ncol(X)
  inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed,
               beta = array(matrix(prior$m0, ncol = ncol(X), nrow = prior$L, byrow = TRUE),
                   dim = c(prior$L, ncol(X))),
               gamma = array(matrix(prior$a, ncol = ncol(X), nrow = prior$L, byrow = TRUE),
                   dim = c(prior$L, ncol(X))),
               z = sample(1:L, n, replace = TRUE),
               Sigmainv = diag(0.1,P),
               alpha = 1
               )
  # inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed)
  # Define monitored parameters
  params_lddp_gmm <- c("beta", "gamma", "w", "alpha", "z")
  
  # Initialize JAGS model
  jags_model <- jags.model(textConnection(lddp_model), 
                          data = data_list, 
                          n.chains = n_chains, 
                          n.adapt = adapt_steps,
                          inits = inits)
  # Burn-in phase
  update(jags_model, mcmc_config$nburn)
  
  # Sampling
  samples <- coda.samples(jags_model, 
                          variable.names = params_lddp_gmm,
                          n.iter = mcmc_config$niter,
                          thin = mcmc_config$thin)
  samples_matrix <- as.matrix(samples)
  w_samples <- samples_matrix[, grep("^w\\[", colnames(samples_matrix))]
  beta_samples <- samples_matrix[, grep("^beta\\[", colnames(samples_matrix))]
  gamma_samples <- samples_matrix[, grep("^gamma\\[", colnames(samples_matrix))]
  alpha_samples <- samples_matrix[, grep("^alpha", colnames(samples_matrix))]
  z_samples <- samples_matrix[, grep("^z\\[", colnames(samples_matrix))]
  # Reshape beta into 3D array (iterations x L x P)
  L <- prior$L
  P <- ncol(X)
  beta_3d <- array(beta_samples, dim = c(nrow(samples_matrix), L, P))
  gamma_3d <- array(gamma_samples, dim = c(nrow(samples_matrix),L, P))
  
  result <- list(
    mcmc_samples = samples,           # Full coda mcmc object
    w_samples = w_samples,            # Matrix of w samples
    beta_samples = beta_3d,           # 3D array of beta samples
    gamma_samples = gamma_3d,         # Matrix of tau samples
    alpha_samples = alpha_samples,    # Matrix of alpha samples
    z_samples = z_samples             # Matrix of z samples
  )
  
  return(result)
}
#--------------------------------------------------------------------------#
#' Covariate-Dependent Variance for the Linear Dirichlet Dependent Mixture Model
#--------------------------------------------------------------------------#
#'
#' This function fits a Linear Dirichlet Dependent Process (LDDP) model 
#' with covariate-dependent mean and variance using Markov Chain Monte Carlo (MCMC) in NIMBLE.
#'
#' @param y Numeric vector of observed response values.
#' @param X Numeric matrix of predictor variables (design matrix).
#' @param prior A list specifying prior hyperparameters:
#'   \itemize{
#'     \item \code{L}: Maximum number of mixture components.
#'     \item \code{m0}: Mean vector for the normal prior on regression coefficients \code{beta}.
#'     \item \code{S0}: Covariance matrix for the normal prior on \code{beta}.
#'     \item \code{a}, \code{b}: Mean vector (\code{a}) and inverse covariance matrix (\code{b}) 
#'           for the normal prior on variance regression coefficients \code{gamma}.
#'     \item \code{Psi}: Scale matrix for the Wishart prior on covariance.
#'     \item \code{nu}: Degrees of freedom for the Wishart prior.
#'     \item \code{aalpha}, \code{balpha}: Shape and rate parameters for the Gamma prior on concentration parameter \code{alpha}.
#'   }
#' @param nburn Number of burn-in iterations to discard before computing results.
#' @param niter Number of MCMC iterations.
#' @param thin Thinning factor for MCMC samples.
#' @param n_chains Number of MCMC chains.
#'
#' @return A list containing:
#' \item{mcmc_samples}{MCMC samples from the posterior distribution as a `coda` object.}
#' \item{w_samples}{Matrix of posterior samples for the mixture weights.}
#' \item{beta_samples}{3D array of posterior samples for regression coefficients \code{beta} (iterations x L x P).}
#' \item{gamma_samples}{3D array of posterior samples for variance coefficients \code{gamma} (iterations x L x P).}
#' \item{alpha_samples}{Matrix of posterior samples for the concentration parameter \code{alpha}.}
#' \item{z_samples}{Matrix of posterior samples for latent component assignments \code{z}.}
lddp_nimble_var <- function(y, X, prior, nburn = 1000, 
                        niter = 5000, thin = 1, 
                        n_chains = 1) {
  # Hierarchy model
  lddp_model <- nimbleCode({
    for (i in 1:n){
        z[i] ~ dcat(w[1:L])
        y[i] ~ dnorm(inprod(x[i, 1:P], beta[z[i], 1:P]), exp(-inprod(x[i,1:P], gamma[z[i], 1:P])))
      }
      for (l in 1:L){
        beta[l, 1:P] ~ dmnorm(abeta[1:P], Sigmainv[1:P, 1:P])
        gamma[l, 1:P] ~ dmnorm(agamma[1:P], bgamma[1:P,1:P])
      }
      abeta[1:P] ~ dmnorm(amu[1:P], b2mu[1:P, 1:P])
      Psi_d[1:P, 1:P] <- nu*psi[1:P, 1:P]
      # dwish in nimble uses precision/scale matrix instead of cov
      Sigmainv[1:P, 1:P] ~ dwish(Psi_d[1:P, 1:P], nu)
      
      # Stick break
      for (l in 1:(L-1)) {
        v[l] ~ dbeta(1, alpha)
      }
      v[L] <- 1
      w[1:L] <- stick_breaking(v[1:(L-1)])
      alpha ~ dgamma(aalpha, balpha)
  })
  data_list <- list(
    y = y,                   # Response variable
    x = X                    # Design matrix
  )
  
constants <- list(
    n = length(y),
    L = prior$L,
    P = ncol(X),
    agamma = prior$a,
    bgamma = prior$b,
    amu = prior$m0,
    b2mu = solve(prior$S0),
    psi = prior$Psi,
    nu = prior$nu,
    aalpha = prior$aalpha,
    balpha = prior$balpha
          )
  L <- prior$L
  P = ncol(X)
  
  
  inits <- list(
  beta = array(0, dim = c(L, P)),
  gamma = array(0, dim = c(L, P)),
  z = sample(1:L, n, replace = TRUE)  # Random initialization for z
  )

  # MCMC confg
  mcmc_config <- list(
    nburn = nburn,
    niter = niter,
    thin = thin
  )
  
  nimble_model <- nimbleModel(code = lddp_model, data = data_list, 
                              constants = constants,
                              inits = inits)

  
  mcmc_conf <- configureMCMC(nimble_model, monitors =c("beta", "gamma", "w", "alpha", "z") )
  
  mcmc_conf$removeSamplers("gamma")
  
  mcmc_conf$addSampler(target = "gamma",
                       type = "RW_block",
                       control = list(scale = 0.1))
  
  # Build the MCMC object
  mcmc <- buildMCMC(mcmc_conf)
  
  # Compile the model and the MCMC object
  compiled_model <- compileNimble(nimble_model)
  compiled_mcmc <- compileNimble(mcmc, project = nimble_model)
  
  # Run the MCMC to obtain samples (after nburn iterations are discarded)
  # note sample size = niter - nburn in nimble, different from jags!
  # Sampling
  samples <- runMCMC(compiled_mcmc, niter = niter, nburnin = nburn, thin = thin)
  # Convert the samples to MCMC format (from matrix to MCMC object)
  samples_mcmc <- as.mcmc(samples)
  w_samples <- samples[, grep("^w\\[", colnames(samples))]
  beta_samples <- samples[, grep("^beta\\[", colnames(samples))]
  gamma_samples <- samples[, grep("^gamma\\[", colnames(samples))]
  alpha_samples <- samples[, grep("^alpha", colnames(samples))]
  z_samples <- samples[, grep("^z\\[", colnames(samples))]
  
  # Reshape beta into 3D array (iterations x L x P)
  L <- prior$L
  P <- ncol(X)
  beta_3d <- array(beta_samples, dim = c(nrow(samples), L, P))
  
  gamma_3d <- array(gamma_samples, dim = c(nrow(samples), L, P))
  
  # return list
  result <- list(
    mcmc_samples = samples_mcmc,    # Full coda mcmc object
    w_samples = w_samples,          # Matrix of w samples
    beta_samples = beta_3d,         # 3D array of beta samples
    gamma_samples = gamma_3d,       # Matrix of tau samples
    alpha_samples = alpha_samples,  # matrix of alpha samples
    z_samples = z_samples           # Matrix of z samples
  )  
  
  return(result)
}

#--------------------------------------------------------------------------#
#' Posterior Prediction for Covariate-Dependent Variance in LDDP Model
#--------------------------------------------------------------------------#
#'
#' This function generates posterior predictions for a Linear Dirichlet Dependent 
#' Process (LDDP) model with covariate-dependent variance using MCMC samples.
#'
#' @param samples MCMC samples from the LDDP_var model.
#' @param X_predict Numeric matrix of predictor variables for new observations.
#' @param L Maximum number of mixture components.
#' @param P Number of covariates (predictors).
#' @param ci_prob Numeric vector specifying the lower and upper quantiles 
#'        for credible intervals (default: `c(0.025, 0.975)`).
#' @param m_true Optional numeric vector of true response values for evaluation.
#'
#' @return A list containing:
#' \item{predicted_mean}{Vector of posterior predictive means.}
#' \item{credible_intervals}{List with lower and upper bounds of credible intervals.}
#' \item{predicted_var}{Vector of posterior predictive variance estimates.}
#' \item{credible_intervals_var}{List with lower and upper bounds of variance credible intervals.}
#' \item{evaluation}{(If `m_true` is provided) A list with: 
#'    \itemize{
#'      \item \code{l2_error}: L2 norm of the prediction error.
#'      \item \code{empirical_coverage}: Proportion of true values covered by the credible interval.
#'      \item \code{mean_ci_length}: Average length of the credible intervals.
#'    }}
prediction_lddp_var_dep <- function(samples, X_predict, L, P, ci_prob = c(0.025, 0.975), m_true = NULL) {
  # m_true, true response values for new data points
  # L: mix comp upper bound
  # P: num of covariates
  # sampels: mcmc samples

  # Convert samples to a matrix for easier manipulation
  samples_matrix <- as.matrix(samples)
  
  # Extract relevant parameters
  w_samples <- samples_matrix[, grep("^w\\[", colnames(samples_matrix))]
  beta_samples <- samples_matrix[, grep("^beta\\[", colnames(samples_matrix))]
  # tau_samples <- samples_matrix[, grep("^tau\\[", colnames(samples_matrix))]
  gamma_samples <- samples_matrix[, grep("^gamma\\[", colnames(samples_matrix))]

  # Reshape beta into 3D array: (iterations x L x P)
  beta_3d <- array(beta_samples, dim = c(nrow(samples_matrix), L, P))

  gamma_3d <- array(gamma_samples, dim = c(nrow(samples_matrix), L, P))
  
  # Number of new predictions
  n_new <- nrow(X_predict)
  
  # Initialize matrix to store predictions
  m_pred <- matrix(0, nrow = n_new, ncol = nrow(samples_matrix))
  
  # Compute predictions for each posterior sample
  for (iter in 1:nrow(samples_matrix)) {
    for (j in 1:n_new) {
      # Weighted sum of predictions across components
      m_pred[j, iter] <- sum(w_samples[iter, ] * X_predict[j, ] %*% t(beta_3d[iter, , ]))
    }
  }
  
  m_var <- matrix(0, nrow = n_new, ncol = nrow(samples_matrix))
  # Variance of mix regression
  for(iter in 1:nrow(samples_matrix)){
    for(j in 1:n_new){
      sigma2_samples <- exp(X_predict[j,] %*% t(gamma_3d[iter, , ]))
      term1 <- sum(w_samples[iter,] * (sigma2_samples + (X_predict[j,] %*% t(beta_3d[iter,,]))^2
                                         ))
      term2 <- (sum(w_samples[iter,] * X_predict[j,] %*% t(beta_3d[iter,,])))^2
      
      m_var[j, iter] <- term1 - term2
    }
  }
  var_mean <- apply(m_var, 1, mean)
  var_l <- apply(m_var, 1, quantile, prob = 0.025)
  var_h <- apply(m_var, 1, quantile, prob = 0.975)
  
  
  # Compute posterior predictive mean and credible intervals
  pred_mean <- apply(m_pred, 1, mean)
  pred_l <- apply(m_pred, 1, quantile, prob = ci_prob[1])
  pred_h <- apply(m_pred, 1, quantile, prob = ci_prob[2])
  
  
  # res
  results <- list(
    predicted_mean = pred_mean,
    credible_intervals = list(lower = pred_l, upper = pred_h),
    predicted_var = var_mean,
    credible_intervals_var = list(lower = var_l, upper=var_h)
  )
  
  # If m_true not null
  if (!is.null(m_true)) {
    # Regression error
    l2_error <- sum(((m_true - pred_mean)^2)/n_new)^.5
    empirical_coverage <- sum((m_true >= pred_l)&(m_true <= pred_h))/n_new
    
    # Evaluation metrics
    results$evaluation <- list(
      l2_error = l2_error,
      empirical_coverage = empirical_coverage,
      mean_ci_length = mean(pred_h - pred_l)  # Average CI length
    )
  }
  
  return(results)
}



inf_criteria_v <-
  function(y, X, res){
    n <- length(y)
    
    if(is.null(res$P)){
      L <- 1
    } else{
      L <- ncol(res$Beta)
    }
    
    if(L > 1){
      p <- res$P
    }
    if(L == 1){
      p <- NULL
    }
    
    beta <- res$Beta
    gamma <- res$gamma
    niter <- nrow(beta)
    
    if(L == 1){
      term <- matrix(0, nrow = niter, ncol = n)
      for(k in 1:niter) {
        sigma2_k <- exp(X%*% gamma[k,])
        term[k,] <- dnorm(y, mean = X%*%beta[k,], sd = sqrt(sigma2_k))
      }        
    }
    
    if(L > 1){
      term_1 <- array(0, c(niter, L, n))
      term <- matrix(0, nrow = niter, ncol = n)
      
      for(i in 1:n) {
        for(l in 1:L) {
          # extract all comps l for ith obs
          sigma2_l <- exp(X[i,] %*% t(gamma[, l,]))
          term_1[,l,i] <- p[,l]*dnorm(y[i], mean = c(X[i,]%*%t(beta[,l,])), sd = sqrt(sigma2_l))
        }
        # adds up over the comp for each iteration
        term[,i] <- apply(term_1[,,i], 1, function(x) sum(x))
      }
    }
    
    term
  }

#--------------------------------------------------------------------------#
#' Compute WAIC for LDDP_var output
#--------------------------------------------------------------------------#
#'
#' This function calculates the Watanabe Akaike Information Criterion (WAIC) 
#' for Covariate-dependent mean and variance LDDP_var 
#' based on MCMC posterior samples.
#'
#' @param y Numeric vector of observed response values.
#' @param X Numeric matrix of predictor variables.
#' @param res A list containing MCMC posterior samples, including:
#'   \itemize{
#'     \item \code{P}: Number of predictors.
#'     \item \code{Beta}: Matrix of posterior samples for regression coefficients.
#'     \item \code{Gamma}: Matrix of posterior samples for regression coefficients of variance.
#'   }
#' @param L Number of mixture components.
#' @param termsum Optional precomputed matrix of likelihood terms 
#'        (default: computed using `inf_criteria_v` function).
#'
#' @return A list containing:
#' \item{pW}{Effective number of parameters in the model.}
#' \item{WAIC}{Computed WAIC value.}
waicnp_v <-
  function(y, X, res, L, termsum = NULL) {
    n <- length(y)  
    p <- res$P
    beta <- res$Beta
    niter <- nrow(p)
    
    if(is.null(termsum)) {
      termsum <- inf_criteria_v(y, X, res)
    }
    logtermsum <- log(termsum)
    
    lpd <- sum(log(apply(exp(logtermsum),2,mean)))
    p2 <- sum(apply(logtermsum,2,var))
    waic <- -2*(lpd-p2)
    
    res <- list()
    res$pW <- p2
    res$WAIC <- waic
    res
  }

#--------------------------------------------------------------------------#
#' Compute Log Pseudo-Marginal Likelihood (LPML) for LDDP_var
#--------------------------------------------------------------------------#
#'
#' This function calculates the Log Pseudo-Marginal Likelihood (LPML) 
#' for Covariate-dependent mean and variance LDDP_var.
#'
#' @param y Numeric vector of observed response values.
#' @param X Numeric matrix of predictor variables.
#' @param res A list containing MCMC posterior samples, including:
#'   \itemize{
#'     \item \code{P}: Number of predictors.
#'     \item \code{Beta}: Matrix of posterior samples for regression coefficients.
#'     \item \code{Gamma}: Matrix of posterior samples for regression coefficients of variance.
#'   }
#' @param L Number of mixture components.
#' @param termsum Optional precomputed matrix of likelihood terms 
#'        (default: computed using `inf_criteria_v` function).
#'
#' @return A list containing:
#' \item{cpo}{Conditional predictive ordinate (CPO) values for each observation.}
#' \item{lpml}{Computed LPML value.}
lpml_v <- function(y, X, res, L, termsum = NULL) {
    n <- length(y)
    p <- res$P
    beta <- res$Beta
    niter <- nrow(beta)
    
    if(is.null(termsum)) {
      termsum <- inf_criteria_v(y, X, res)
    }
    
    aux <- 1/termsum
    omegabari <- apply(aux, 2, mean)
    omegabari_1 <- sqrt(niter) * omegabari
    omegatilde <- matrix(0, nrow = niter, ncol = n)
    
    for(i in 1:n) {
      omegatilde[,i] <- pmin(aux[,i], omegabari_1[i])  
    }
    
    sum_omegatilde <- apply(omegatilde,2,sum)
    sum_term_omegatilde <- apply(termsum*omegatilde, 2, sum)
    cpo <- sum_term_omegatilde/sum_omegatilde
    
    lpml <- sum(log(cpo))
    
    res <- list()
    res$cpo <- cpo
    res$lpml <- lpml
    res 
}
