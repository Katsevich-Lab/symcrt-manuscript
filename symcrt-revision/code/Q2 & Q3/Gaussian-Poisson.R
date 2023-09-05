# computation of nu in different settings
library(symcrt)
library(simulatr)
library(katlabutils)
library(glmnet)
library(ggplot2)

## X|Z binomial; Y|Z Gaussian; varying dimension
parameter_grid <- data.frame(
  d = seq(11, 15, by = 1),      # dimension of Z
  n = 100                      # sample size
)

parameter_grid

get_ground_truth <- function(d){
  gamma <- rep(1, d)
  beta <- gamma
  return(list(gamma = gamma,
              beta = beta))
}

parameter_grid <- parameter_grid |> add_ground_truth(get_ground_truth)

parameter_grid


fixed_parameters <- list(
  B = 1e4,                      # number of data realizations
  seed = 4                      # seed to set prior to generating data and running methods
)


# define data-generating model based on the Gaussian linear model
generate_data_f <- function(n, d, ground_truth){
  Z <- katlabutils::fast_generate_mvn(rep(0, d), covariance =  diag(d), num_samples = n)
  X <- rnorm(n = n, mean = Z%*%ground_truth$gamma, sd = 1)
  Y <- rpois(n = n, lambda = exp(rep(-4, n) + Z%*%ground_truth$beta))
  
  data <- list(X = X, Y = Y, Z = Z)
  data
}


# need to call simulatr_function() to give simulatr a few more pieces of info
generate_data_function <- simulatr_function(
  f = generate_data_f,
  arg_names = formalArgs(generate_data_f),
  loop = TRUE
)


# GCM test with lasso methods
gcm_test <- function(data){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(data$X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  list(p_value <- 2*pnorm(abs(z_stat), lower.tail = FALSE))
}


# dCRT with X sampling and lasso methods
dCRT_X <- function(data, B_small = 1000){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  null_z_stats <- numeric(B_small)
  for(b in 1:B_small){
    X_resampled <- rnorm(n = n, mean = X_on_Z_fit$fitted.values, sd = 1)
    R_resampled <- (X_resampled - X_on_Z_predict)*(Y - Y_on_Z_predict)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
  }
  null_z_stats <- c(z_stat, null_z_stats)
  list(p_value <- 2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat)))
}

# dCRT with Y sampling and lasso methods
dCRT_Y <- function(data, B_small = 1000){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  null_z_stats <- numeric(B_small)
  for(b in 1:B_small){
    Y_resampled <- rpois(n = n, lambda = Y_on_Z_fit$fitted.values)
    R_resampled <- (X - X_on_Z_predict)*(Y_resampled - Y_on_Z_predict)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
  }
  null_z_stats <- c(z_stat, null_z_stats)
  list(p_value <- 2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat)))
}


# create simulatr functions
gcm_test_spec_f <- simulatr_function(f = gcm_test, arg_names = character(0), loop = TRUE)
dCRT_X_spec_f <- simulatr_function(f = dCRT_X, arg_names = character(0), loop = TRUE)
dCRT_Y_spec_f <- simulatr_function(f = dCRT_Y, arg_names = character(0), loop = TRUE)


run_method_functions <- list(gcm_test = gcm_test_spec_f,
                             dCRT_X = dCRT_X_spec_f,
                             dCRT_Y = dCRT_Y_spec_f)


type_I_err_mean <- function(output, ground_truth) {
  mean(output <= 0.05)
}

evaluation_functions <- list(type_I_err_mean = type_I_err_mean)


simulatr_spec <- simulatr_specifier(
  parameter_grid,
  fixed_parameters,
  generate_data_function,
  run_method_functions,
  evaluation_functions
)


check_results <- check_simulatr_specifier_object(simulatr_spec, B_in = 1e4)

# sim_results <- check_simulatr_specifier_object(simulatr_spec)

check_results$metrics

type_I_err <- check_results$metrics |>
  ggplot(aes(x = d,
             y = mean,
             ymin = mean - 2*se,
             ymax = mean + 2*se,
             color = method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(width = 0.5) +
  geom_hline(yintercept=0.05, linetype=2, color = "red") +
  labs(x = "Dimension of covariate Z",
       y = "Type I Error",
       title = "Gaussian-Poisson model") +
  theme(legend.position = "bottom") 


ggsave(filename = "varying-dimension-gaussian-poisson.png",
       plot = type_I_err,
       device = "png",
       width = 5,
       height = 5)

# save the rds result
saveRDS(check_results, "Q2_gaussian_poisson.rds")

