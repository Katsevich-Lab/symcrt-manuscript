
library(simulatr)
library(katlabutils)
library(glmnet)
library(ggplot2)


parameter_grid <- data.frame(
  d = 100,                     # dimension of Z
  s = seq(5, 65, by = 15),     # number of nonzero coefficients for X|Z model
  n = 100                      # sample size
)

parameter_grid

get_ground_truth <- function(n, d, s){
  candidate_gamma <- numeric(d)
  gamma <- rep(1, d)
  candidate_gamma[1:s] <- gamma[1:s]
  
  candidate_beta <- numeric(d)
  beta <- gamma
  candidate_beta[1:s] <- beta[1:s]
  return(list(gamma = candidate_gamma,
              beta = candidate_beta))
}

parameter_grid <- parameter_grid |> add_ground_truth(get_ground_truth)

parameter_grid


fixed_parameters <- list(
  B = 500,                      # number of data realizations
  seed = 4                      # seed to set prior to generating data and running methods
)


# define data-generating model based on the Gaussian linear model
generate_data_f <- function(n, d, ground_truth){
  Z <- katlabutils::fast_generate_mvn(rep(0, d), covariance =  diag(d), num_samples = n)
  X <- rnorm(n = n, mean = Z%*%ground_truth$gamma, sd = 1)
  Y <- rbinom(n = n, size = 1, prob = exp(Z%*%ground_truth$beta) / (1 + exp(Z%*%ground_truth$beta)))
  
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
gcm_test_lasso <- function(data){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run lasso result
  Y_on_Z_fit <- glmnet::cv.glmnet(Z, Y, family = "binomial", thresh = 1e-5)
  X_on_Z_fit <- glmnet::cv.glmnet(Z, X, family = "gaussian")
  Y_on_Z_fit_lasso <- predict(Y_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  X_on_Z_fit_lasso <- predict(X_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  R <- (X - X_on_Z_fit_lasso)*(Y - Y_on_Z_fit_lasso)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  list(p_value <- 2*pnorm(abs(z_stat), lower.tail = FALSE))
}

# GCM test with post-lasso methods
gcm_test_refit <- function(data){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run lasso result
  Y_on_Z_fit <- glmnet::cv.glmnet(Z, Y, family = "binomial", thresh = 1e-5)
  X_on_Z_fit <- glmnet::cv.glmnet(Z, X, family = "gaussian")
  Y_on_Z_fit_lasso <- predict(Y_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  X_on_Z_fit_lasso <- predict(X_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  # run post-lasso on X|Z model
  act_set <- which(as.vector(coef(X_on_Z_fit, s = "lambda.min"))[-1] != 0)
  if(length(act_set) == 0){
    X_on_Z_refit <- rep(mean(X), n)
    X_on_Z_fit_refit <- X_on_Z_refit
  }else{
    if(length(act_set) >= 100){
      act_set <- act_set[1:99]
    }
    Z_subset <- Z[, act_set]
    X_on_Z_refit <- lm(X ~ Z_subset)
    X_on_Z_fit_refit <- unname(X_on_Z_refit$fitted.values)
  }
  
  # test
  R <- (X - X_on_Z_fit_refit)*(Y - Y_on_Z_fit_lasso)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  list(p_value <- 2*pnorm(abs(z_stat), lower.tail = FALSE))
}

# dCRT with X sampling and lasso methods
dCRT_X_lasso <- function(data, B_small = 1000){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run lasso result
  Y_on_Z_fit <- glmnet::cv.glmnet(Z, Y, family = "binomial", thresh = 1e-5)
  X_on_Z_fit <- glmnet::cv.glmnet(Z, X, family = "gaussian")
  Y_on_Z_fit_lasso <- predict(Y_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  X_on_Z_fit_lasso <- predict(X_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  R <- (X - X_on_Z_fit_lasso)*(Y - Y_on_Z_fit_lasso)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  null_z_stats <- numeric(B_small)
  gamma_lasso <- as.vector(coef(X_on_Z_fit, s = "lambda.min"))[-1]
  for(b in 1:B_small){
    X_resampled <- rnorm(n = n, mean = Z%*%gamma_lasso, sd = 1)
    R_resampled <- (X_resampled - X_on_Z_fit_lasso)*(Y - Y_on_Z_fit_lasso)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
  }
  null_z_stats <- c(z_stat, null_z_stats)
  list(p_value <- 2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat)))
}

# dCRT with Y sampling and lasso methods
dCRT_Y_lasso <- function(data, B_small = 1000){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run lasso result
  Y_on_Z_fit <- glmnet::cv.glmnet(Z, Y, family = "binomial", thresh = 1e-5)
  X_on_Z_fit <- glmnet::cv.glmnet(Z, X, family = "gaussian")
  Y_on_Z_fit_lasso <- predict(Y_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  X_on_Z_fit_lasso <- predict(X_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  R <- (X - X_on_Z_fit_lasso)*(Y - Y_on_Z_fit_lasso)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  null_z_stats <- numeric(B_small)
  beta_lasso <- as.vector(coef(Y_on_Z_fit, s = "lambda.min"))[-1]
  for(b in 1:B_small){
    Y_resampled <- rbinom(n = n, size = 1, prob = exp(Z%*%beta_lasso) / (1 + exp(Z%*%beta_lasso)))
    R_resampled <- (X - X_on_Z_fit_lasso)*(Y_resampled - Y_on_Z_fit_lasso)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
  }
  null_z_stats <- c(z_stat, null_z_stats)
  list(p_value <- 2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat)))
}


# dCRT with X sampling and post-lasso methods
dCRT_X_refit <- function(data, B_small = 1000){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run lasso result
  Y_on_Z_fit <- glmnet::cv.glmnet(Z, Y, family = "binomial", thresh = 1e-5)
  X_on_Z_fit <- glmnet::cv.glmnet(Z, X, family = "gaussian")
  Y_on_Z_fit_lasso <- predict(Y_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  X_on_Z_fit_lasso <- predict(X_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  # run post-lasso on X|Z model
  act_set <- which(as.vector(coef(X_on_Z_fit, s = "lambda.min"))[-1] != 0)
  if(length(act_set) == 0){
    X_on_Z_refit <- rep(mean(X), n)
    X_on_Z_fit_refit <- X_on_Z_refit
    # obtain the resample mean
    resampled_mean <- rep(0, n)
  }else{
    if(length(act_set) >= 100){
      act_set <- act_set[1:99]
    }
    Z_subset <- Z[, act_set]
    X_on_Z_refit <- lm(X ~ Z_subset)
    X_on_Z_fit_refit <- unname(X_on_Z_refit$fitted.values)
    # obtain the resample mean
    gamma_refit <- as.vector(coef(X_on_Z_refit))[-1]
    resampled_mean <- Z_subset%*%gamma_refit
  }
  
  # test
  R <- (X - X_on_Z_fit_refit)*(Y - Y_on_Z_fit_lasso)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  null_z_stats <- numeric(B_small)
  for(b in 1:B_small){
    X_resampled <- rnorm(n = n, mean = resampled_mean, sd = 1)
    R_resampled <- (X_resampled - X_on_Z_fit_refit)*(Y - Y_on_Z_fit_lasso)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
  }
  null_z_stats <- c(z_stat, null_z_stats)
  list(p_value <- 2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat)))
}

# dCRT with Y sampling and post-lasso methods
dCRT_Y_refit <- function(data, B_small = 1000){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run lasso result
  Y_on_Z_fit <- glmnet::cv.glmnet(Z, Y, family = "binomial", thresh = 1e-5)
  X_on_Z_fit <- glmnet::cv.glmnet(Z, X, family = "gaussian")
  Y_on_Z_fit_lasso <- predict(Y_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  X_on_Z_fit_lasso <- predict(X_on_Z_fit, Z,
                              s = "lambda.min",
                              type = "response")
  # run post-lasso on X|Z model
  act_set <- which(as.vector(coef(X_on_Z_fit, s = "lambda.min"))[-1] != 0)
  if(length(act_set) == 0){
    X_on_Z_refit <- rep(mean(X), n)
    X_on_Z_fit_refit <- X_on_Z_refit
  }else{
    if(length(act_set) >= 100){
      act_set <- act_set[1:99]
    }
    Z_subset <- Z[, act_set]
    X_on_Z_refit <- lm(X ~ Z_subset)
    X_on_Z_fit_refit <- unname(X_on_Z_refit$fitted.values)
  }
  
  # test
  R <- (X - X_on_Z_fit_refit)*(Y - Y_on_Z_fit_lasso)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  null_z_stats <- numeric(B_small)
  beta_lasso <- as.vector(coef(Y_on_Z_fit, s = "lambda.min"))[-1]
  for(b in 1:B_small){
    Y_resampled <- rbinom(n = n, size = 1, prob = exp(Z%*%beta_lasso) / (1 + exp(Z%*%beta_lasso)))
    R_resampled <- (X - X_on_Z_fit_refit)*(Y_resampled - Y_on_Z_fit_lasso)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
  }
  null_z_stats <- c(z_stat, null_z_stats)
  list(p_value <- 2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat)))
}

# create simulatr functions
gcm_test_lasso_spec_f <- simulatr_function(f = gcm_test_lasso, arg_names = character(0), loop = TRUE)
gcm_test_refit_spec_f <- simulatr_function(f = gcm_test_refit, arg_names = character(0), loop = TRUE)
dCRT_X_lasso_spec_f <- simulatr_function(f = dCRT_X_lasso, arg_names = character(0), loop = TRUE)
dCRT_Y_lasso_spec_f <- simulatr_function(f = dCRT_Y_lasso, arg_names = character(0), loop = TRUE)
dCRT_X_refit_spec_f <- simulatr_function(f = dCRT_X_refit, arg_names = character(0), loop = TRUE)
dCRT_Y_refit_spec_f <- simulatr_function(f = dCRT_Y_refit, arg_names = character(0), loop = TRUE)


run_method_functions <- list(gcm_test_lasso = gcm_test_lasso_spec_f,
                             gcm_test_refit = gcm_test_refit_spec_f,
                             dCRT_X_lasso = dCRT_X_lasso_spec_f,
                             dCRT_Y_lasso = dCRT_Y_lasso_spec_f,
                             dCRT_X_refit = dCRT_X_refit_spec_f,
                             dCRT_Y_refit = dCRT_Y_refit_spec_f)


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


check_results <- check_simulatr_specifier_object(simulatr_spec, B_in = 500)

# sim_results <- check_simulatr_specifier_object(simulatr_spec)

check_results$metrics

type_I_err <- check_results$metrics |>
  ggplot(aes(x = s,
             y = mean,
             ymin = mean - 2*se,
             ymax = mean + 2*se,
             color = method)) +
  geom_point() +
  geom_line() +
  geom_errorbar(width = 2) +
  geom_hline(yintercept=0.05, linetype=2, color = "red") +
  labs(x = "Sparsity of X|Z and Y|Z model",
       y = "Type I Error") +
  theme(legend.position = "bottom") 


ggsave(filename = "varying-both-sparsity-gaussian-binomial.png",
       plot = type_I_err,
       device = "png",
       width = 5,
       height = 5)

# save the rds result
saveRDS(check_results, 
        "~/Dropbox (Penn)/CRT-revision/symcrt-revision/code/result/Q2_varying_both_sparsity.rds")


