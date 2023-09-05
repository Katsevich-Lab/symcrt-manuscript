# computation of nu in different settings
library(symcrt)
library(simulatr)
library(katlabutils)
library(glmnet)
library(ggplot2)
library(symcrt)


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
  oracle_Y <- exp(rep(-4, n) + Z%*%ground_truth$beta)
  oracle_X <- Z %*% ground_truth$gamma
  
  data <- list(X = X, Y = Y, Z = Z, oracle_X = oracle_X, oracle_Y = oracle_Y)
  data
}


# need to call simulatr_function() to give simulatr a few more pieces of info
generate_data_function <- simulatr_function(
  f = generate_data_f,
  arg_names = formalArgs(generate_data_f),
  loop = TRUE
)


# GCM test with lasso methods
X_on_Z_regression <- function(data, ground_truth){
  X <- data$X
  Z <- data$Z
  oracle_X <- data$oracle_X
  n <- length(data$X)
  # run glm result
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  
  list(MSE_X = sum((X_on_Z_predict - oracle_X)^2)/length(X_on_Z_predict))
}


# GCM test with lasso methods
Y_on_Z_regression <- function(data, ground_truth){
  Y <- data$Y
  Z <- data$Z
  oracle_Y <- data$oracle_Y
  n <- length(data$X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  
  list(MSE_Y = sum((Y_on_Z_predict - oracle_Y)^2)/length(Y_on_Z_predict))
}

# create simulatr functions
MSE_X_spec_f <- simulatr_function(f = X_on_Z_regression, arg_names = character(0), loop = TRUE)
MSE_Y_spec_f <- simulatr_function(f = Y_on_Z_regression, arg_names = character(0), loop = TRUE)


run_method_functions <- list(MSE_X_on_Z_regression = MSE_X_spec_f,
                             MSE_Y_on_Z_regression = MSE_Y_spec_f)


MSE_mean <- function(output, ground_truth) {
  mean(as.numeric(output))
}

evaluation_functions <- list(MSE_mean = MSE_mean)


simulatr_spec <- simulatr_specifier(
  parameter_grid,
  fixed_parameters,
  generate_data_function,
  run_method_functions,
  evaluation_functions
)


check_results <- check_simulatr_specifier_object(simulatr_spec, B_in = 1e4)

# sim_results <- check_simulatr_specifier_object(simulatr_spec)

renamed_result <- check_results$metrics %>% rename("Regression" = "method")

MSE <- renamed_result |>
  ggplot(aes(x = d,
             y = mean,
             ymin = mean - 2*se,
             ymax = mean + 2*se,
             color = Regression)) +
  scale_y_log10() +
  geom_point() +
  geom_line() +
  geom_errorbar(width = 0.5) +
  labs(x = "Dimension of covariate Z",
       y = "MSE",
       title = "Gaussian-Poisson model") +
  theme(legend.position = "bottom") 


ggsave(filename = "MSE-gaussian-poisson.png",
       plot = MSE,
       device = "png",
       width = 5,
       height = 5)

