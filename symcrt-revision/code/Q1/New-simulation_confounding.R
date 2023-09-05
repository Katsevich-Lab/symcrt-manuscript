# computation of nu in different settings
library(symcrt)
library(simulatr)
library(katlabutils)
library(glmnet)
library(ggplot2)
B <- 1e6

## X|Z binomial; Y|Z Poisson;

# define expit function
expit <- function(theta)(exp(theta)/(1+exp(theta)))

# data-generating function based on negative binomial distribution
generate_data_pois <- function(n, gamma_0, gamma_1, beta_0, beta_1){
  Z <- rnorm(n = n, mean = 0, sd = 1)
  X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
  Y <- rpois(n = n, lambda = exp(beta_0 + beta_1*Z))
  list(X = X, Y = Y, Z = Z)
}

n <- 1000                  # sample size
gamma_0 <- -4              # intercept in X on Z model
gamma_1 <- 1               # slope in X on Z model
beta_0 <- -3               # intercept in Y on Z model
beta_1 <- 1                # slope in Y on Z model

data <- generate_data_pois(B, gamma_0, gamma_1, beta_0, beta_1)
X <- data$X
Y <- data$Y
print(symcrt::simulate_confounding(n, X, Y))
