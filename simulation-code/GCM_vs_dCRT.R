# load libraries
library(katlabutils)
library(MASS)
library(reshape2)
library(ggplot2)
library(tidyverse)

# define expit function
expit <- function(theta)(exp(theta)/(1+exp(theta)))

# data-generating function based on negative binomial distribution
generate_data_pois <- function(n, gamma_0, gamma_1, beta_0, beta_1){
  Z <- rnorm(n = n, mean = 0, sd = 1)
  X <- rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
  Y <- rpois(n = n, lambda = exp(beta_0 + beta_1*Z))
  list(X = X, Y = Y, Z = Z)
}

# GCM test
gcm_test <- function(X, Y, Z){
  n <- length(X)
  Y_on_Z_fit <- glm(Y ~ Z,
                    family = "poisson")
  X_on_Z_fit <- glm(X ~ Z, family = "binomial")
  R <- (X - X_on_Z_fit$fitted.values)*(Y - Y_on_Z_fit$fitted.values)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  2*pnorm(abs(z_stat), lower.tail = FALSE)
}

# dCRT
dCRT <- function(X, Y, Z, B){
  n <- length(X)
  Y_on_Z_fitted <- glm(Y ~ Z,
                       family = "poisson")$fitted.values
  X_on_Z_fitted <- glm(X ~ Z, family = "binomial")$fitted.values
  R <- (X - X_on_Z_fitted)*(Y - Y_on_Z_fitted)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  null_z_stats <- numeric(B)
  for(b in 1:B){
    X_resampled <- rbinom(n = n, size = 1, prob = X_on_Z_fitted)
    R_resampled <- (X_resampled - X_on_Z_fitted)*(Y - Y_on_Z_fitted)
    null_z_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
  }
  null_z_stats <- c(z_stat, null_z_stats)
  2*min(mean(null_z_stats >= z_stat), mean(null_z_stats <= z_stat))
}

# simulation parameters
n <- 1000                  # sample size
reps_calibration <- 1000   # number of reps used for calibration assessment
reps_timing <- 5           # number of reps used for timing assessment
gamma_0 <- -4              # intercept in X on Z model
gamma_1 <- 1               # slope in X on Z model
beta_0 <- -3               # intercept in Y on Z model
beta_1 <- 1                # slope in Y on Z model
B_small <- 1000            # number of dCRT resamples for calibration assessment
B_large <- 100000          # number of dCRT resamples for timing assessment

# calibration assessment
methods <- c("GCM", "dCRT")
pvals <- matrix(data = NA,
                nrow = reps_calibration,
                ncol = length(methods),
                dimnames = list(rep = 1:reps_calibration, method = methods))
set.seed(1)
for(rep in 1:reps_calibration){
  data <- generate_data_pois(n, gamma_0, gamma_1, beta_0, beta_1)
  pvals[rep, "GCM"] <- gcm_test(data$X, data$Y, data$Z)
  pvals[rep, "dCRT"] <- dCRT(data$X, data$Y, data$Z, B_small)
}

# create the result directory and save the results as RDS file
result_dir <- "simulation-results/GCM_vs_dCRT"
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
  cat("Directory created:", result_dir, "\n")
} else {
  cat("Directory already exists:", result_dir, "\n")
}

# save the p-values
saveRDS(pvals, sprintf("%s/pvals.rds", result_dir))

# compute GCM sampling test statistics
set.seed(1)
result <- data.frame(sampled_stats = rep(0, reps_calibration),
                     resampled_stats = rep(0, reps_calibration))
for (i in 1:reps_calibration) {
  data <- generate_data_pois(n, gamma_0, gamma_1, beta_0, beta_1)
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)
  Y_on_Z_fitted <- glm(Y ~ Z,
                       family = "poisson")$fitted.values
  X_on_Z_fitted <- glm(X ~ Z, family = "binomial")$fitted.values
  R <- (X - X_on_Z_fitted)*(Y - Y_on_Z_fitted)
  z_stat <- sum(1/sqrt(n)*sum(R))/sd(R)
  result$sampled_stats[i] <- z_stat
}
for(b in 1:B_small){
  X_resampled <- rbinom(n = n, size = 1, prob = X_on_Z_fitted)
  R_resampled <- (X_resampled - X_on_Z_fitted)*(Y - Y_on_Z_fitted)
  result$resampled_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))/sd(R_resampled)
}

# save the result
saveRDS(result, sprintf("%s/GCM_test_stat.rds", result_dir))
