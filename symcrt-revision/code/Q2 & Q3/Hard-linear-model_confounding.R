# computation of nu in different settings
library(symcrt)
library(simulatr)
library(katlabutils)
library(glmnet)
library(ggplot2)
B <- 2e5

## X|Z Gaussian; Y|Z binomial; varying sparsity in X|Z model;
parameter_grid <- data.frame(
  d = 100,                     # dimension of Z
  s_x = seq(5, 65, by = 15),   # number of nonzero coefficients for X|Z model
  n = 100,                     # sample size
  s_y = 5                      # number of nonzero coefficients for Y|Z model
)

parameter_grid

get_ground_truth <- function(n, d, s_x, s_y){
  candidate_gamma <- numeric(d)
  gamma <- rep(1, d)
  candidate_gamma[1:s_x] <- gamma[1:s_x]
  
  candidate_beta <- numeric(d)
  beta <- gamma
  candidate_beta[1:s_y] <- beta[1:s_y]
  return(list(gamma = candidate_gamma,
              beta = candidate_beta))
}

parameter_grid <- parameter_grid |> add_ground_truth(get_ground_truth)


# define data-generating model based on the Gaussian linear model
generate_data_f <- function(n, d, ground_truth){
  Z <- katlabutils::fast_generate_mvn(rep(0, d), covariance =  diag(d), num_samples = n)
  X <- rnorm(n = n, mean = Z%*%ground_truth$gamma, sd = 1)
  Y <- rbinom(n = n, size = 1, prob = exp(Z%*%ground_truth$beta) / (1 + exp(Z%*%ground_truth$beta)))
  
  data <- list(X = X, Y = Y, Z = Z)
  data
}


confounding_list <- numeric(nrow(parameter_grid))
for (i in 1:nrow(parameter_grid)) {
  data <- generate_data_f(B, parameter_grid$d[i], parameter_grid$ground_truth[[i]])
  X <- data$X
  Y <- data$Y
  print(symcrt::simulate_confounding(parameter_grid$n[i], X, Y))
}
