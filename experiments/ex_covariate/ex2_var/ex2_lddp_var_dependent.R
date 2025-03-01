source("code/Bayes_lddp_var.R")
require(ggplot2)
require(mcclust.ext)
require(nor1mix)
require(splines)


set.seed(1010101)
n <- 400
p <- 1
x <- matrix(0, n, p)
y <- rep(0, n)

x[, 1] <- runif(n, -2, 10)

m  <-  rep(0, n)
m[x <= 2] <- 0
m[x > 2 & x <= 5] <- 2*(x[x > 2 & x <=5 ] -2)
m[x > 5] <- 6
eps <- rep(0, n)
eps[x <= 2] <- rnorm(sum(x <= 2), 0, 0.2)
eps[x > 2 & x <= 5] <- rnorm(sum(x > 2 & x <= 5), 0, 0.05)
eps[x > 5] <- rnorm(sum(x > 5), 0, (x[x > 5] - 5)^2/15 + 0.01)
y <- m + eps

# Prediction
n_new <- 800
x_new <- matrix(0, n_new, p)

x_new[, 1] <- seq(-2, 10, 12/(n_new - 1))

m_true_new <- rep(0, n_new)
m_true_new[x_new > 2 & x_new <= 5] <- 2*(x_new[x_new > 2 & x_new <= 5] - 2)
m_true_new[x_new > 5] <- 6
# plot(x,y)
y_grid <- seq(-1, 10, .02)
m2 <- length(y_grid)
f_true_new <- matrix(0, m2, n_new)
f_true_new[, x_new <= 2] <- dnorm(matrix(y_grid, nrow = m2, ncol = sum(x_new <= 2)),
                                  t(matrix(m_true_new[x_new <= 2],
                                           nrow = sum(x_new <= 2), ncol = m2)), 0.2)
f_true_new[, x_new > 2 & x_new <= 5] <- dnorm(matrix(y_grid, nrow = m2, 
                                                     ncol = sum(x_new > 2 & x_new <= 5)),
                                              t(matrix(m_true_new[x_new > 2 & x_new <= 5],
                                                       nrow = sum(x_new > 2 & x_new <= 5),
                                                       ncol = m2)), 0.05)
f_true_new[, x_new > 5] <- dnorm(matrix(y_grid, nrow = m2, ncol = sum(x_new > 5)),
                                 t(matrix(m_true_new[x_new > 5],
                                          nrow = sum(x_new > 5), ncol = m2)),
                                 t(matrix((x_new[x_new > 5] - 5)^2/15 + 0.01,
                                          nrow = sum(x_new > 5), ncol = m2)))
# -------------------
# Non B-spline
# -------------------
# LDDP model
# Design and Prediction matrices
X <- cbind(rep(1, n), x[, 1])
k <- ncol(X)
X_predict <-cbind(rep(1, n_new), x_new[, 1])

# Data-driven prior
fit_lm <- lm(y ~ X - 1)
prior_v <- list(m0 = fit_lm$coefficients,
              S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2), 
              nu = k + 2, 
              Psi = 30*(solve(t(X)%*%X)*(sigma(fit_lm)^2)),
              a = rep(0,ncol(X)),
              b = diag(1,ncol(X)),
              aalpha = 1,
              balpha = 1,
              L = 20)
res_lm_square <- residuals(fit_lm)^2
fit_gamma <- lm(log(res_lm_square) ~ X -1 )
prior_v <- list(m0 = fit_lm$coefficients,
              S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2), 
              nu = k + 2, 
              Psi = 30*(solve(t(X)%*%X)*(sigma(fit_lm)^2)),
              a = fit_gamma$coefficients,
              b = diag(1,ncol(X)),
              aalpha = 1,
              balpha = 1,
              L = 20)


# fitting lddp
# seed = 345
seed = 123
fit_lddp_v <- lddp_var(y, X, prior_v, nburn = 5000, niter = 5000, adapt_steps = 1000, seed = seed)

pred_v <- prediction_lddp_var_dep(fit_lddp_v$mcmc_samples, X_predict, L = 20, P = k,
                        ci_prob = c(0.025, 0.975))

df_lddp_ex2_v <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- pred_v$predicted_mean,
  lower_c <- pred_v$credible_intervals$lower,
  upper_c <- pred_v$credible_intervals$upper
)

ggplot(df_lddp_ex2_v, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15)+ theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) 

#waic and lpml
samples_res_ex2_v <- list(z = fit_lddp_v$z_samples,
                    P = fit_lddp_v$w_samples,
                    Beta = fit_lddp_v$beta_samples,
                    gamma = fit_lddp_v$gamma_samples)


waic_lddp_ex2_v <- waicnp_v(y, X, res = samples_res_ex2_v, L = 20, termsum = NULL)
waic_lddp_ex2_v #552.0286

lpml_lddp_ex2_v <- lpml_v(y, X,res = samples_res_ex2_v, L = 20, termsum = NULL)
lpml_lddp_ex2_v #-283.6174


# plotting variance
true_var <- rep(0, length(x_new))

# Variance for x <= 2
true_var[x_new <= 2] <- 0.2^2

# Variance for 2 < x <= 5
true_var[x_new > 2 & x_new <= 5] <- 0.05^2

# Variance for x > 5
true_var[x_new > 5] <- ((x_new[x_new > 5] - 5)^2 / 15 + 0.01)^2

# plotting variance
var_mean_v <- pred_v$predicted_var
var_l_v <- pred_v$credible_intervals_var$lower
var_h_v <- pred_v$credible_intervals_var$upper

df_lddp_ex2_v <- data.frame(
  x = x_new,            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_mean_v,  # Posterior mean variance
  lower_var = var_l_v,    # Lower credible interval
  upper_var = var_h_v     # Upper credible interval
)

ggplot(df_lddp_ex2_v, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))


# -------------------
# knot 0
# -------------------
knots0 <- c()
X0 <- cbind(rep(1, n),
           bs(x[, 1], degree = 3, knots = knots0, intercept = FALSE)
)
k0 <- ncol(X0)

X0_predict <-cbind(rep(1, n_new),
                  predict(bs(x[, 1], degree = 3, knots = knots0, intercept = FALSE),
                          x_new[, 1])
)
fit_lm0 <- lm(y ~ X0 - 1)
res_lm_square0 <- residuals(fit_lm0)^2
fit_gamma0 <- lm(log(res_lm_square0) ~ X0 -1 )
prior0_v <- list(m0 = fit_lm0$coefficients,
               S0 = solve(t(X0)%*%X0)*(sigma(fit_lm0)^2), 
               nu = k0 + 2, 
               Psi = 30*(solve(t(X0)%*%X0)*(sigma(fit_lm0)^2)),
               a = fit_gamma0$coefficients,
               b = diag(1,ncol(X0)),
               aalpha = 1,
               balpha = 1, 
               L = 20)

fit_lddpbs0_v <- lddp_var(y, X0, prior0_v, nburn = 5000, niter = 10000, adapt_steps = 5000, seed = 123)

predbs0_v <- prediction_lddp_var_dep(fit_lddpbs0_v$mcmc_samples, X0_predict, L = 20, P = k0,
                        ci_prob = c(0.025, 0.975))
df_lddpbs0_ex2_v <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- predbs0_v$predicted_mean,
  lower_c <- predbs0_v$credible_intervals$lower,
  upper_c <- predbs0_v$credible_intervals$upper
)

ggplot(df_lddpbs0_ex2_v, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15)+ theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))


#waic and lpml
samples_res_ex20_v <- list(z = fit_lddpbs0_v$z_samples,
                    P = fit_lddpbs0_v$w_samples,
                    Beta = fit_lddpbs0_v$beta_samples,
                    gamma = fit_lddpbs0_v$gamma_samples)

waic_lddpbs0_ex2_v <- waicnp_v(y, X0, res = samples_res_ex20_v, L = 20, termsum = NULL)
waic_lddpbs0_ex2_v #1062.952

lpml_lddp_ex2_v <- lpml_v(y,X0,res = samples_res_ex20_v, L = 20, termsum = NULL)
lpml_lddp_ex2_v #-530.8388


# plotting variance
true_var <- rep(0, length(x_new))

# Variance for x <= 2
true_var[x_new <= 2] <- 0.2^2

# Variance for 2 < x <= 5
true_var[x_new > 2 & x_new <= 5] <- 0.05^2

# Variance for x > 5
true_var[x_new > 5] <- ((x_new[x_new > 5] - 5)^2 / 15 + 0.01)^2

# plotting variance
var_mean0_v <- predbs0_v$predicted_var
var_l0_v <- predbs0_v$credible_intervals_var$lower
var_h0_v <- predbs0_v$credible_intervals_var$upper

df_lddp_ex20_v <- data.frame(
  x = x_new,            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_mean0_v,  # Posterior mean variance
  lower_var = var_l0_v,    # Lower credible interval
  upper_var = var_h0_v     # Upper credible interval
)

ggplot(df_lddp_ex20_v, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15)+ theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))


# -------------------
# knot 5
# -------------------
knots5 <- quantile(x[, 1], c(1/6,2/6,3/6,4/6,5/6))
X5 <- cbind(rep(1, n),
            bs(x[, 1], degree = 3, knots = knots5, intercept = FALSE)
)
k5 <- ncol(X5)

X5_predict <-cbind(rep(1, n_new),
                   predict(bs(x[, 1], degree = 3, knots = knots5, intercept = FALSE),
                           x_new[, 1])
)
fit_lm5 <- lm(y ~ X5 - 1)
res_lm_square5 <- residuals(fit_lm5)^2
fit_gamma5 <- lm(log(res_lm_square5) ~ X5 -1 )
prior5_v <- list(m0 = fit_lm5$coefficients,
               S0 = solve(t(X5)%*%X5)*(sigma(fit_lm5)^2), 
               nu = k5 + 2, 
               Psi = 30*(solve(t(X5)%*%X5)*(sigma(fit_lm5)^2)),
               a = fit_gamma5$coefficients,
               b = diag(1,ncol(X5)),
               aalpha = 1,
               balpha = 1, 
               L = 20)
# set.seed(123)
fit_lddp5 <- lddp_var(y, X5, prior5_v, nburn = 5000, niter = 10000, adapt_steps = 5000, seed = 123)
pred5 <- prediction_lddp_var_dep(fit_lddp5$mcmc_samples, X5_predict, L = 20, P = ncol(X5),
                        ci_prob = c(0.025, 0.975))

df_lddp_ex25 <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- pred5$predicted_mean,
  lower_c <- pred5$credible_intervals$lower,
  upper_c <- pred5$credible_intervals$upper
)

require(ggplot2)
ggplot(df_lddp_ex25, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15)+ theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))

# waic and lpml for cov dependent var
samples_res_v <- list(z = fit_lddp5$z_samples,
                    P = fit_lddp5$w_samples,
                    Beta = fit_lddp5$beta_samples,
                    gamma = fit_lddp5$gamma_samples)

dim(samples_res_v$gamma)
waic_v <- waicnp_v(y, X5, samples_res_v, L = 20)
waic_v 
lpml_v <- lpml_v(y, X5, samples_res_v, L = 20)
lpml_v

# plotting variance
true_var <- rep(0, length(x_new))

# Variance for x <= 2
true_var[x_new <= 2] <- 0.2^2

# Variance for 2 < x <= 5
true_var[x_new > 2 & x_new <= 5] <- 0.05^2

# Variance for x > 5
true_var[x_new > 5] <- ((x_new[x_new > 5] - 5)^2 / 15 + 0.01)^2


# plotting variance
var_mean5 <- pred5$predicted_var
var_l5 <- pred5$credible_intervals_var$lower
var_h5 <- pred5$credible_intervals_var$upper

df_lddp_ex25_var <- data.frame(
  x = x_new,            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_mean5,  # Posterior mean variance
  lower_var = var_l5,    # Lower credible interval
  upper_var = var_h5     # Upper credible interval
)

ggplot(df_lddp_ex25_var, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15)+ theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))
