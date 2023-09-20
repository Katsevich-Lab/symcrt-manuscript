library(ggplot2)
library(tidyverse)
# unnormalized case
n <- 200
d <- 16 
B_small <- 1e3
ground_truth <- list(beta = rep(1, d), gamma = rep(1, d))
result <- data.frame(sampled_unnormalized_stats = rep(0, B_small),
                     resampled_X_unnormalize_stats = rep(0, B_small),
                     resampled_Y_unnormalize_stats = rep(0, B_small))

# generate normal-poisson data
generate_data_f <- function(n, d, ground_truth){
  Z <- katlabutils::fast_generate_mvn(rep(0, d), covariance =  diag(d), num_samples = n)
  X <- rnorm(n = n, mean = Z%*%ground_truth$gamma, sd = 1)
  Y <- rpois(n = n, lambda = exp(rep(-4, n) + Z%*%ground_truth$beta))
  
  data <- list(X = X, Y = Y, Z = Z)
  data
}

# compute the unnormalized sampled statistics
set.seed(1)
for (i in 1:B_small) {
  data <- generate_data_f(n, d, ground_truth)
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  result$sampled_unnormalize_stats[i] <- sum(1/sqrt(n)*sum(R))
}
set.seed(1)
data <- generate_data_f(n, d, ground_truth)
X <- data$X; Y <- data$Y; Z <- data$Z
n <- length(X)
# run glm result
Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")

# compute the unnormalized resampled statistics
for(b in 1:B_small){
  # resample from Y|Z
  Y_resampled <- rpois(n = n, lambda = Y_on_Z_fit$fitted.values)
  R_resampled <- (X - X_on_Z_predict)*(Y_resampled - Y_on_Z_predict)
  result$resampled_Y_unnormalize_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))
  
  # resample form X|Z
  X_resampled <- rnorm(n = n, mean = X_on_Z_fit$fitted.values, sd = 1)
  R_resampled <- (X_resampled - X_on_Z_predict)*(Y - Y_on_Z_predict)
  result$resampled_X_unnormalize_stats[b] <- sum(1/sqrt(n)*sum(R_resampled))
}

# comparison histogram between normalized sampled/resampled statistics
unnormalized_statistics_comparison <- tibble(result) |>
  pivot_longer(cols = c("sampled_unnormalize_stats", "resampled_Y_unnormalize_stats", 
                        "resampled_X_unnormalize_stats"),
               values_to = "statistics",
               names_to = "category") |>
  ggplot() + 
  geom_histogram(aes(x = statistics, color = category), position = "stack") +
  # stat_bin(aes(x = statistics, color = category), bins = 20) +
  labs(x = "Test statistics", 
       y = "Density",
       title = "Resampled versus sampled normalized statistics when resampling Y|Z") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_discrete(name="")

# comparison density curve between normalized sampled/resampled statistics
unnormalized_statistics_density_comparison <- tibble(result) |>
  pivot_longer(cols = c("sampled_unnormalize_stats", "resampled_Y_unnormalize_stats", 
                        "resampled_X_unnormalize_stats"),
               values_to = "statistics",
               names_to = "category") |>
  ggplot() + 
  geom_density(aes(x = statistics, color = category)) +
  labs(x = "Test statistics", 
       y = "Density",
       title = "Resampled versus sampled unnormalized statistics when resampling Y|Z") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_discrete(name="")

# plot the empirical CDF
emp_CDF_unnormalized_statistics <- tibble(result) |>
  pivot_longer(cols = c("sampled_unnormalize_stats", "resampled_Y_unnormalize_stats", 
                        "resampled_X_unnormalize_stats"),
               values_to = "statistics",
               names_to = "category") |> 
  mutate(category = as.factor(category)) |>
  mutate(category = fct_recode(category, 
                               "sampled_stats" = "sampled_unnormalize_stats",
                               "resampled_Y_stats" = "resampled_Y_unnormalize_stats",
                               "resampled_X_stats" = "resampled_X_unnormalize_stats")) |>
  ggplot(aes_string(x = "statistics", colour = "category")) + 
  stat_ecdf() +
  labs(x = "Value of test statistics", 
       y = "Empirical CDF",
       title = "Empirical CDF comparison for unnormalized statistics") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_discrete(name="")

# save the empirical CDF of normalized statistics
ggsave("figures/emp_CDF_unnormalized_statistics.pdf",
  plot = emp_CDF_unnormalized_statistics,
  width = 5,
  height = 4)

# normalized case
result <- data.frame(sampled_normalize_stats = rep(0, B_small),
                     resampled_Y_normalize_stats = rep(0, B_small),
                     resampled_Y_normalize_stats = rep(0, B_small))

# compute the normalized sampled statistics
set.seed(1)
for (i in 1:B_small) {
  data <- generate_data_f(n, d, ground_truth)
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)
  # run glm result
  Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  result$sampled_normalize_stats[i] <- sum(1/sqrt(n)*sum(R)) / sd(R)
}
set.seed(1)
data <- generate_data_f(n, d, ground_truth)
X <- data$X; Y <- data$Y; Z <- data$Z
n <- length(X)
# run glm result
Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")

# compute the normalized resampled statistics
for(b in 1:B_small){
  Y_resampled <- rpois(n = n, lambda = Y_on_Z_fit$fitted.values)
  R_resampled <- (X - X_on_Z_predict)*(Y_resampled - Y_on_Z_predict)
  result$resampled_Y_normalize_stats[b] <- sum(1/sqrt(n)*sum(R_resampled)) / sd(R_resampled)
  
  # resample form X|Z
  X_resampled <- rnorm(n = n, mean = X_on_Z_fit$fitted.values, sd = 1)
  R_resampled <- (X_resampled - X_on_Z_predict)*(Y - Y_on_Z_predict)
  result$resampled_X_normalize_stats[b] <- sum(1/sqrt(n)*sum(R_resampled)) / sd(R_resampled)
}

# comparison histogram between normalized sampled/resampled statistics
normalized_statistics_comparison <- tibble(result) |>
  pivot_longer(cols = c("sampled_normalize_stats", "resampled_Y_normalize_stats", 
                        "resampled_X_normalize_stats"),
               values_to = "statistics",
               names_to = "category") |>
  ggplot() + 
  geom_histogram(aes(x = statistics, color = category), position = "stack") +
  # stat_bin(aes(x = statistics, color = category), bins = 20) +
  labs(x = "Test statistics", 
       y = "Density",
       title = "Resampled versus sampled normalized statistics when resampling Y|Z") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_discrete(name="")

# comparison density curve between normalized sampled/resampled statistics
normalized_statistics_density_comparison <- tibble(result) |>
  pivot_longer(cols = c("sampled_normalize_stats", "resampled_Y_normalize_stats", 
                        "resampled_X_normalize_stats"),
               values_to = "statistics",
               names_to = "category") |> 
  mutate(category = as.factor(category)) |>
  mutate(category = fct_recode(category, 
                               "sampled_stats" = "sampled_normalize_stats",
                               "resampled_Y_stats" = "resampled_Y_normalize_stats",
                               "resampled_X_stats" = "resampled_X_normalize_stats")) |>
  ggplot() + 
  geom_density(aes(x = statistics, color = category)) +
  labs(x = "Test statistics", 
       y = "Density",
       title = "Resampled versus sampled normalized statistics when resampling Y|Z") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_discrete(name="")


# plot the empirical CDF for normalized statistics
emp_CDF_normalized_statistics <- tibble(result) |>
  pivot_longer(cols = c("sampled_normalize_stats", "resampled_Y_normalize_stats", 
                        "resampled_X_normalize_stats"),
               values_to = "statistics",
               names_to = "category") |> 
  mutate(category = as.factor(category)) |>
  mutate(category = fct_recode(category, 
                               "sampled_stats" = "sampled_normalize_stats",
                               "resampled_Y_stats" = "resampled_Y_normalize_stats",
                               "resampled_X_stats" = "resampled_X_normalize_stats")) |>
  ggplot(aes_string(x = "statistics", colour = "category")) + 
  stat_ecdf() +
  labs(x = "Value of test statistics", 
       y = "Empirical CDF",
       title = "Empirical CDF comparison for normalized statistics") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_discrete(name="")

# save the empirical CDF of normalized statistics
ggsave("figures/emp_CDF_normalized_statistics.pdf",
       plot = emp_CDF_normalized_statistics,
       width = 5,
       height = 4)
