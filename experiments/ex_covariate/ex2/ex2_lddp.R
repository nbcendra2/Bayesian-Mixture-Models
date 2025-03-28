source("code/Bayes_lddp.R")
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
# S0 cov matrix, di inverse di dalem function
prior <- list(m0 = fit_lm$coefficients,
              S0 = solve(t(X)%*%X)*(sigma(fit_lm)^2), 
              nu = k + 2, 
              Psi = 30*(solve(t(X)%*%X)*(sigma(fit_lm)^2)),
              a = 2,
              b = (sigma(fit_lm)^2)/2,
              aalpha = 1,
              balpha = 1,
              L = 20)
dim(prior$Psi)
# fitting lddp
fit_lddp <- lddp(y, X, prior, nburn = 5000, niter = 5000, seed = 123)

# calculate E[Y|x] and Var(Y|x)
pred <- prediction_lddp(fit_lddp$mcmc_samples, X_predict, L = 20, P = k,
                        ci_prob = c(0.025, 0.975))

df_lddp_ex2 <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- pred$predicted_mean,
  lower_c <- pred$credible_intervals$lower,
  upper_c <- pred$credible_intervals$upper
)

ggplot(df_lddp_ex2, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 8)

#waic and lpml
samples_res_ex2 <- list(z = fit_lddp$z_samples,
                    P = fit_lddp$w_samples,
                    Beta = fit_lddp$beta_samples,
                    Sigma2 = 1/fit_lddp$tau_samples)

waic_lddp_ex2 <- waicnp(y, X, res = samples_res_ex2, L = 20, termsum = NULL)
waic_lddp_ex2 

lpml_lddp_ex2 <- lpml(y,X,res = samples_res_ex2, L = 20, termsum = NULL)
lpml_lddp_ex2 


# plotting variance
true_var <- rep(0, length(x_new))

# Variance for x <= 2
true_var[x_new <= 2] <- 0.2^2

# Variance for 2 < x <= 5
true_var[x_new > 2 & x_new <= 5] <- 0.05^2

# Variance for x > 5
true_var[x_new > 5] <- ((x_new[x_new > 5] - 5)^2 / 15 + 0.01)^2

# plotting variance
var_mean <- pred$predicted_var
var_l <- pred$credible_intervals_var$lower
var_h <- pred$credible_intervals_var$upper

df_lddp_ex2_var <- data.frame(
  x = x_new,            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_mean,  # Posterior mean variance
  lower_var = var_l,    # Lower credible interval
  upper_var = var_h     # Upper credible interval
)

ggplot(df_lddp_ex2_var, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 8)


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

df_lddp_ex20 <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- pred0$predicted_mean,
  lower_c <- pred0$credible_intervals$lower,
  upper_c <- pred0$credible_intervals$upper
)

# plot
ggplot(df_lddp_ex20, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15))   + ylim(NA, 8)

#waic and lpml
samples_res_ex20 <- list(z = fit_lddp0$z_samples,
                    P = fit_lddp0$w_samples,
                    Beta = fit_lddp0$beta_samples,
                    Sigma2 = 1/fit_lddp0$tau_samples)

waic_lddp_ex20 <- waicnp(y, X0, res = samples_res_ex20, L = 20, termsum = NULL)
waic_lddp_ex20 

lpml_lddp_ex20 <- lpml(y,X0,res = samples_res_ex20, L = 20, termsum = NULL)
lpml_lddp_ex20 

# plotting variance
var_mean0 <- pred0$predicted_var
var_l0 <- pred0$credible_intervals_var$lower
var_h0 <- pred0$credible_intervals_var$upper

df_lddp_ex20_var <- data.frame(
  x = x_new,            # x values (e.g., grid)
  true_var = true_var,  # True variance
  mean_var = var_mean0,  # Posterior mean variance
  lower_var = var_l0,    # Lower credible interval
  upper_var = var_h0     # Upper credible interval
)

ggplot(df_lddp_ex20_var, aes(x = x))+ geom_line(aes(y = true_var), color = "black") +
  geom_line(aes(y = mean_var), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_var, ymax = upper_var), fill = "dodgerblue", alpha = 0.3) + labs(
     x = "x1", 
     y = "Var(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 8)


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
prior5 <- list(m0 = fit_lm5$coefficients,
               S0 = solve(t(X5)%*%X5)*(sigma(fit_lm5)^2), 
               nu = k5 + 2, 
               Psi = 30*(solve(t(X5)%*%X5)*(sigma(fit_lm5)^2)),
               a = 2,
               b = (sigma(fit_lm5)^2)/2,
               aalpha = 1,
               balpha = 1, 
               L = 20)

fit_lddp5 <- lddp(y, X5, prior5, nburn = 5000, niter = 5000, adapt_steps = 1000, seed = 123)
pred5 <- prediction_lddp(fit_lddp5$mcmc_samples, X5_predict, L = 20, P = ncol(X5),
                        ci_prob = c(0.025, 0.975))
df_lddp_ex25 <- data.frame(
  x = x_new[,1],
  true_y <- m_true_new,
  mean_y <- pred5$predicted_mean,
  lower_c <- pred5$credible_intervals$lower,
  upper_c <- pred5$credible_intervals$upper
)

ggplot(df_lddp_ex25, aes(x = x))+ geom_line(aes(y = true_y), color = "black") +
  geom_line(aes(y = mean_y), color = "blue", linetype = "dashed") + 
  geom_ribbon(aes(ymin = lower_c, ymax = upper_c), fill = "dodgerblue", alpha = 0.3) + 
  labs(
     x = "x1", 
     y = "E(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 8)

#waic and lpml
samples_res_ex25 <- list(z = fit_lddp5$z_samples,
                    P = fit_lddp5$w_samples,
                    Beta = fit_lddp5$beta_samples,
                    Sigma2 = 1/fit_lddp5$tau_samples)

waic_lddp_ex25 <- xwaicnp(y, X5, res = samples_res_ex25, L = 20, termsum = NULL)
waic_lddp_ex25 #52.478

lpml_lddp_ex25 <- lpml(y,X5,res = samples_res_ex25, L = 20, termsum = NULL)
lpml_lddp_ex25 #-26.78229

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
     y = "Var(Y | x)") + theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5, face="bold"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), text = element_text(size = 15)) + ylim(NA, 8)

                                                                                     