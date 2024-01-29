# computation of nu in different settings
library(simulatr)
library(katlabutils)
library(glmnet)
library(ggplot2)
library(dplyr)

## X|Z binomial; Y|Z Gaussian; varying dimension
parameter_grid <- data.frame(
  d = seq(10, 18, by = 2),      # dimension of Z
  n = 200                       # sample size
)

get_ground_truth <- function(d){
  gamma <- rep(0.1, d)
  beta <- rep(0.2, d)*3
  return(list(gamma = gamma,
              beta = beta))
}

parameter_grid <- parameter_grid |> add_ground_truth(get_ground_truth)

fixed_parameters <- list(
  B = 1e4,                      # number of data realizations
  seed = 4                      # seed to set prior to generating data and running methods
)


# define data-generating model based on the Gaussian linear model
generate_data_f <- function(n, d, ground_truth){
  Z <- katlabutils::fast_generate_mvn(rep(0, d), covariance =  diag(d), num_samples = n)
  X <- rpois(n = n, lambda = exp(rep(-1, n) + Z%*%ground_truth$gamma))
  Y <- rpois(n = n, lambda = exp(rep(-1, n) + Z%*%ground_truth$beta))
  
  data <- list(X = X, Y = Y, Z = Z)
  data
}


# need to call simulatr_function() to give simulatr a few more pieces of info
generate_data_function <- simulatr_function(
  f = generate_data_f,
  arg_names = formalArgs(generate_data_f),
  loop = TRUE
)

# dCRT with X sampling and lasso methods
dCRT_X <- function(data, B_small = 1000){
  X <- data$X
  Y <- data$Y
  Z <- data$Z
  n <- length(X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "poisson"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  z_stat <- sum(1/sqrt(n)*sum(R))
  null_z_stats <- numeric(B_small)
  for(b in 1:B_small){
    X_resampled <- rpois(n = n, lambda = X_on_Z_fit$fitted.values)
    R_resampled <- (X_resampled - X_on_Z_predict)*(Y - Y_on_Z_predict)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))
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
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "poisson"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  z_stat <- sum(1/sqrt(n)*sum(R))
  null_z_stats <- numeric(B_small)
  for(b in 1:B_small){
    Y_resampled <- rpois(n = n, lambda = Y_on_Z_fit$fitted.values)
    R_resampled <- (X - X_on_Z_predict)*(Y_resampled - Y_on_Z_predict)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))
  }
  null_z_stats <- c(z_stat, null_z_stats)
  list(p_value <- 2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat)))
}



# create simulatr functions
dCRT_X_spec_f <- simulatr_function(f = dCRT_X, arg_names = character(0), loop = TRUE)
dCRT_Y_spec_f <- simulatr_function(f = dCRT_Y, arg_names = character(0), loop = TRUE)


run_method_functions <- list(dCRT_X = dCRT_X_spec_f,
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


check_results <- simulatr::check_simulatr_specifier_object(simulatr_spec, B_in = 5e3)

# create the result directory and save the results as RDS file
result_dir <- "simulation-results/asymmetry_investigation"
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
  cat("Directory created:", result_dir, "\n")
} else {
  cat("Directory already exists:", result_dir, "\n")
}
saveRDS(check_results, sprintf("%s/check_results_double_poisson.rds", 
                               result_dir))

