library(tidyverse)
n <- 200
d <- 16
B_small <- 1e3
ground_truth <- list(beta = rep(1, d), gamma = rep(1, d))
result <- data.frame(sampled_var = rep(0, B_small),
                     resampled_X_var = rep(0, B_small),
                     resampled_Y_var = rep(0, B_small))

# generate normal-poisson data
generate_data_f <- function(n, d, ground_truth){
  Z <- katlabutils::fast_generate_mvn(rep(0, d), covariance =  diag(d), num_samples = n)
  X <- rnorm(n = n, mean = Z%*%ground_truth$gamma, sd = 1)
  Y <- rpois(n = n, lambda = exp(rep(-4, n) + Z%*%ground_truth$beta))
  
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
  X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
  Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
  X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
  R <- (X - X_on_Z_predict)*(Y - Y_on_Z_predict)
  result$sampled_var[i] <- var(R)
}

# compute the resampled sd
set.seed(1)
data <- generate_data_f(n, d, ground_truth)
X <- data$X; Y <- data$Y; Z <- data$Z
n <- length(X)
# run glm result
Y_on_Z_fit <- suppressWarnings(glm(Y ~ Z, family = "poisson"))
X_on_Z_fit <- suppressWarnings(glm(X ~ Z, family = "gaussian"))
Y_on_Z_predict <- predict(Y_on_Z_fit, as.data.frame(Z), type = "response")
X_on_Z_predict <- predict(X_on_Z_fit, as.data.frame(Z), type = "response")
for(b in 1:B_small){
  # resample from Y|Z
  Y_resampled <- rpois(n = n, lambda = Y_on_Z_fit$fitted.values)
  R_resampled <- (X - X_on_Z_predict)*(Y_resampled - Y_on_Z_predict)
  result$resampled_Y_var[b] <- var(R_resampled)
  
  # resample from X|Z
  X_resampled <- rnorm(n = n, mean = X_on_Z_fit$fitted.values, sd = 1)
  R_resampled <- (X_resampled - X_on_Z_predict)*(Y - Y_on_Z_predict)
  result$resampled_X_var[b] <- var(R_resampled)
}

# plot the density of sample sd versus resample sd
var_comparison <- tibble(result) |>
  mutate(normal_stats = rnorm(B_small)) |>
  pivot_longer(cols = c("sampled_var", "resampled_X_var", "resampled_Y_var"),
               values_to = "var",
               names_to = "category") |>
  ggplot() + 
  geom_density(aes(x = var, color = category)) +
  scale_x_log10() +
  labs(x = "variance", 
       y = "density",
       title = "Resampling variance versus sampling variance") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_color_discrete(name="")

# save the plot
ggsave("figures/var_comparison.pdf",
       plot = var_comparison,
       width = 7,
       height = 5)
