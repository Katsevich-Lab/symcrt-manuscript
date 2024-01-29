library(tidyverse)
n <- 200
d <- 16
B_small <- 1e3
ground_truth <- list(beta = rep(0.2, d)*3, gamma = rep(0.1, d))
result_var <- data.frame(sampled_var = rep(0, B_small),
                         resampled_X_var = rep(0, B_small),
                         resampled_Y_var = rep(0, B_small))

result_teststat <- data.frame(sampled_teststat = rep(0, B_small),
                              resampled_X_teststat = rep(0, B_small),
                              resampled_Y_teststat = rep(0, B_small),
                              sampled_nteststat = rep(0, B_small),
                              resampled_X_nteststat = rep(0, B_small),
                              resampled_Y_nteststat = rep(0, B_small))

# generate normal-poisson data
generate_data_f <- function(n, d, ground_truth){
  Z <- katlabutils::fast_generate_mvn(rep(0, d), covariance =  diag(d), num_samples = n)
  X <- rpois(n = n, lambda = exp(rep(-1.5, n) + Z%*%ground_truth$beta))
  Y <- rpois(n = n, lambda = exp(rep(-1.5, n) + Z%*%ground_truth$beta))
  
  data <- list(X = X, Y = Y, Z = Z)
  data
}

# compute the sampled sd
set.seed(1)
for (i in 1:B_small) {
  data <- generate_data_f(n, d, ground_truth)
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "poisson"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict) 
  result_teststat$sampled_teststat[i] <- sum(R) / sqrt(n)
}

# compute the resampled sd
set.seed(1)
data <- generate_data_f(n, d, ground_truth)
X <- data$X; Y <- data$Y; Z <- data$Z
n <- length(X)
# run glm result
Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "poisson"))
Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
for(b in 1:B_small){
  # resample from Y|Z
  Y_resampled <- rpois(n = n, lambda = Y_on_Z_fit$fitted.values)
  R_resampled <- (X - X_on_Z_predict)*(Y_resampled - Y_on_Z_predict)
  result_teststat$resampled_Y_teststat[b] <- sum(R_resampled) / (sqrt(n))

  # resample from X|Z
  X_resampled <- rpois(n = n, lambda = X_on_Z_fit$fitted.values)
  R_resampled <- (X_resampled - X_on_Z_predict)*(Y - Y_on_Z_predict)
  result_teststat$resampled_X_teststat[b] <- sum(R_resampled) / (sqrt(n))
}

# plot the density of sample sd versus resample sd
my_theme <-   theme_bw() + 
  theme(
    legend.position = "bottom",
    axis.line = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(
      color = "black"
    ),
    plot.title = element_blank(),
    legend.title = element_blank()
  )

# create the directory and save the result
result_dir <- "simulation-results/asymmetry_investigation"
if (!dir.exists(result_dir)) {
  dir.create(result_dir)
  cat("Directory created:", result_dir, "\n")
} else {
  cat("Directory already exists:", result_dir, "\n")
}

saveRDS(result_teststat, 
        sprintf("%s/result_teststat.rds", result_dir))
