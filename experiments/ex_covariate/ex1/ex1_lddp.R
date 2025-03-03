source("code/Bayes_lddp.R")
require(ggplot2)
require(mcclust.ext)
require(nor1mix)
require(splines)

set.seed(1010101)
n <- 200
p <- 2
x <-  matrix(0, n, p)

x[,1] <- runif(n, -1, 8)
x[,2] <- (x[, 1] - 3.5)^2/3 - 1 + rnorm(n, 0, 0.05)

y <- 5 - log(x[, 1] + 2) + rnorm(n, 0, 0.05)
# data free noise, y without epsilon noise
m_true <- 5 - log(x[, 1] + 2)

# Prediction dataset
n_new <- 800
x_new <- matrix(0, n_new, p)

x_new[, 1] <- runif(n_new, -1, 8)
x_new[, 2] <- (x_new[, 1] - 3.5)^2/3 - 1 + rnorm(n_new, 0, 0.05)

m_true_new <-  5-log(x_new[, 1] + 2)


# True density
y_grid <- seq(2.4, 5.3, .02)
m2 <- length(y_grid)
f_true_new <- dnorm(matrix(y_grid, nrow = m2, ncol = n_new), t(matrix(m_true_new, nrow = n_new, ncol = m2)), 0.05)

### LDDP model
# Design and prediction matrices
# Design matrix with covariate xi1 and xi2 as made above, with intercept 1
X <- cbind(rep(1, n),
           x[, 1],
           x[, 2]
           )
# Num of covariates, 3
k <- ncol(X)
X_predict <-cbind(rep(1, n_new),
                  x_new[, 1],
                  x_new[, 2]
                  )

# -------------------
# Non B-spline
# -------------------
# LDDP model

# Data-driven prior
fit_lm <- lm(y ~ X - 1)
prior <- list(m0 = fit_lm$coefficients,
              S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2), 
              nu = k + 2, 
              Psi = 30*(solve(t(X)%*%X)*(sigma(fit_lm)^2)),
              a = 2,
              b = (sigma(fit_lm)^2)/2,
              aalpha = 1,
              balpha = 1,
              L = 20)

# fitting lddp
fit_lddp <- lddp(y, X, prior, nburn = 5000, niter = 5000, seed = 123)

# calculate E[Y|x] and Var(Y|x)
pred <- prediction_lddp(fit_lddp$mcmc_samples, X_predict, L = 20, P = k,
                        ci_prob = c(0.025, 0.975))


df_lddp_ex1 <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- pred$predicted_mean,
  lower_c <- pred$credible_intervals$lower,
  upper_c <- pred$credible_intervals$upper
)

ggplot(df_lddp_ex1, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 5)

#waic and lpml
samples_res_ex1 <- list(z = fit_lddp$z_samples,
                    P = fit_lddp$w_samples,
                    Beta = fit_lddp$beta_samples,
                    Sigma2 = 1/fit_lddp$tau_samples)

waic_lddp_ex1 <- waicnp(y, X, res = samples_res_ex1, L = 20, termsum = NULL)
waic_lddp_ex1 

lpml_lddp_ex1 <- lpml(y,X,res = samples_res_ex1, L = 20, termsum = NULL)
lpml_lddp_ex1 

# plotting variance
true_var <- rep(0.0025, n_new)
# plotting variance
var_mean <- pred$predicted_var
var_l <- pred$credible_intervals_var$lower
var_h <- pred$credible_intervals_var$upper

df_lddp_ex1_var <- data.frame(
  x = x_new[,1],            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_mean,  # Posterior mean variance
  lower_var = var_l,    # Lower credible interval
  upper_var = var_h     # Upper credible interval
)

ggplot(df_lddp_ex1_var, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 0.05)



# -------------------
# knot 0
# -------------------
knots0 <- c()
X0 <- cbind(rep(1, n),
            bs(x[ , 1], degree = 3, knots = knots0, intercept = FALSE),
            bs(x[ , 2], degree = 3, knots = knots0, intercept = FALSE)
)

X0_predict <-cbind(rep(1, n_new),
                   predict(bs(x[, 1], degree = 3, knots = knots0, intercept = FALSE),
                           x_new[, 1]),
                   predict(bs(x[, 2], degree = 3, knots = knots0, intercept = FALSE), 
                           x_new[, 2])
)
k0 <- ncol(X0)

fit_lm0 <- lm(y ~ X0 - 1)
prior0 <- list(m0 = fit_lm0$coefficients,
               S0 = solve(t(X0)%*%X0)*(sigma(fit_lm0)^2), 
               nu = k0 + 2, 
               Psi = 30*(solve(t(X0)%*%X0)*(sigma(fit_lm0)^2)),
               a = 2,
               b = (sigma(fit_lm0)^2)/2,
               aalpha = 1,
               balpha = 1,
               L = 20)

fit_lddp0 <- lddp(y, X0, prior0, nburn = 5000, niter = 5000, adapt_steps = 1000, seed = 123)
pred0 <- prediction_lddp(fit_lddp0$mcmc_samples, X0_predict, L = 20, P = ncol(X0),
                        ci_prob = c(0.025, 0.975))

df_lddp_ex10 <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- pred0$predicted_mean,
  lower_c <- pred0$credible_intervals$lower,
  upper_c <- pred0$credible_intervals$upper
)

# plot
ggplot(df_lddp_ex10, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))   + ylim(NA, 5)


#waic and lpml
samples_res_ex10 <- list(z = fit_lddp0$z_samples,
                    P = fit_lddp0$w_samples,
                    Beta = fit_lddp0$beta_samples,
                    Sigma2 = 1/fit_lddp0$tau_samples)

waic_lddp_ex10 <- waicnp(y, X0, res = samples_res_ex10, L = 20, termsum = NULL)
waic_lddp_ex10 

lpml_lddp_ex10 <- lpml(y,X0,res = samples_res_ex10, L = 20, termsum = NULL)
lpml_lddp_ex10 

# plotting variance
var_mean0 <- pred0$predicted_var
var_l0 <- pred0$credible_intervals_var$lower
var_h0 <- pred0$credible_intervals_var$upper

df_lddp_ex10_var <- data.frame(
  x = x_new[,1],            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_mean0,  # Posterior mean variance
  lower_var = var_l0,    # Lower credible interval
  upper_var = var_h0     # Upper credible interval
)

ggplot(df_lddp_ex10_var, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 0.05)


# -------------------
# Bayesian Linear Regression
# -------------------
lin_reg_model <- "model{
  for (i in 1:n){
    y[i] ~ dnorm(inprod(x[i,], beta[1:P]), tau)
  }
  # Mulvar
  beta[1:P] ~ dmnorm(abeta, Sigmainv[,])
  
  tau ~ dgamma(asigma2, bsigma2)
  sigma2 <- 1/tau
}"
fit_norm2 <- lm(y ~ X -1)
prior_linreg <- list(m0 = fit_norm2$coefficients,
                     S0 = solve(solve(t(X)%*%X)*(sigma(fit_norm2)^2)),
                     a = 2,
                     b = (sigma(fit_norm2)^2)/2)
# set.seed
inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 123)

data_linreg <- list(y = y, 
                    x = X, 
                    n = n, 
                    abeta = prior_linreg$m0,
                    Sigmainv = solve(prior_linreg$S0),
                    asigma2 = prior_linreg$a,
                    bsigma2 = prior_linreg$b,
                    P = k)
model_linreg = jags.model(textConnection(lin_reg_model), n.chains = 1, 
                   data = data_linreg, inits = inits)
update(model_linreg, 5000)
samples_linreg <- coda.samples(model_linreg, variable.names = c("beta","sigma2"), n.iter = 5000)
samples_linreg_mat <- as.matrix(samples_linreg)
# grab samples
beta_linreg <- samples_linreg_mat[, grep("^beta\\[", colnames(samples_linreg_mat))]
sigma2_linreg <- samples_linreg_mat[,"sigma2"]
# matrix for fitted values/exp
m_pred_linreg <- matrix(0, nrow = 800, ncol = 5000)
for(i in 1: 5000){
  m_pred_linreg[,i] <- beta_linreg[i,1] + beta_linreg[i,2] * X_predict[,2] + beta_linreg[i,3] * X_predict[,3]
}

m_pred_linreg_m <- apply(m_pred_linreg, 1, mean)
m_pred_linreg_l <- apply(m_pred_linreg, 1, quantile, prob = 0.025)
m_pred_linreg_h <- apply(m_pred_linreg, 1, quantile, prob = 0.975)

df_linreg_ex1 <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- m_pred_linreg_m,
  lower_c <- m_pred_linreg_l,
  upper_c <- m_pred_linreg_m
)

ggplot(df_linreg_ex1, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 5)


# compute waic and lpml for normal linear regression
samples_res_linreg <- list(z = NULL,
                    P = NULL,
                    Beta = beta_linreg,
                    Sigma2 = sigma2_linreg)

waic_linreg_ex1 <- waicnp(y, X, samples_res_linreg, L = 1)
waic_linreg_ex1 #-492.8921

lpml_linreg_ex1 <- lpml(y, X, samples_res_linreg, L = 1)
lpml_linreg_ex1 #246.4346


# plotting conditional variance
var_meang <- mean(sigma2_linreg) 
var_lg <- quantile(sigma2_linreg, probs = 0.025)  
var_hg <- quantile(sigma2_linreg, probs = 0.975) 

df_linreg_ex1_var <- data.frame(
  x = x_new[,1],            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_meang,  # Posterior mean variance
  lower_var = var_lg,    # Lower credible interval
  upper_var = var_hg     # Upper credible interval
)

ggplot(df_linreg_ex1_var, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 0.05)

