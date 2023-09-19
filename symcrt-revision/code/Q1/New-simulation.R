## question 1
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


# timing assessment
computation_times <- matrix(data = NA,
                            nrow = reps_timing,
                            ncol = length(methods),
                            dimnames = list(rep = 1:reps_timing, method = methods))
for(rep in 1:reps_timing){
  data <- generate_data_pois(n, gamma_0, gamma_1, beta_0, beta_1)
  computation_times[rep, "GCM"] <- system.time(
    gcm_test(data$X, data$Y, data$Z)
  )[["elapsed"]]
  computation_times[rep, "dCRT"] <- system.time(
    dCRT(data$X, data$Y, data$Z, B_large)
  )[["elapsed"]]
}


# create QQ plots of null p-values
qq_plots <- pvals |>
  melt(value.name = "p-value") |>
  ggplot(aes(y = `p-value`, color = method)) +
  stat_qq_points() +
  stat_qq_band() +
  geom_abline() +
  scale_x_continuous(trans = revlog_trans(), breaks = c(1, 0.1, 0.01, 0.001)) +
  scale_y_continuous(trans = revlog_trans(), breaks = c(1, 1e-2, 1e-4, 1e-6, 1e-8)) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "Expected null p-value",
       y = "Observed p-value")

ggsave(filename = "qq-plots.pdf",
       plot = qq_plots,
       device = "pdf",
       width = 3.5,
       height = 3.5)

# Type-I-error comparison across different methods with sig. level = 5e-2
Type_I_error <- data.frame(mean = rep(0, 2),
                           se = rep(0, 2),
                           method = c("GCM", "dCRT"))
Type_I_error$mean <- apply(pvals, 2, function(x) mean(x <= 0.05))
Type_I_error$se <- apply(pvals, 2, function(x) sd(x <= 0.05)) / sqrt(reps_calibration)


type_I_err_2 <- Type_I_error |>
  ggplot(aes(x = method,
             y = mean,
             ymin = mean - 2*se,
             ymax = mean + 2*se)) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  geom_hline(yintercept=0.05, linetype=2, color = "red") +
  labs(x = "Inferential methods",
       y = "Type I Error") +
  theme(legend.position = "bottom") 

ggsave(filename = "type_I_err_5e-2.pdf",
       plot = type_I_err_2,
       device = "pdf",
       width = 3.5,
       height = 3.5)


# Type-I-error comparison across different methods with sig. level = 5e-3
Type_I_error$mean <- apply(pvals, 2, function(x) mean(x <= 5e-3))
Type_I_error$se <- apply(pvals, 2, function(x) sd(x <= 5e-3)) / sqrt(reps_calibration)


type_I_err_3 <- Type_I_error |>
  ggplot(aes(x = method,
             y = mean,
             ymin = mean - 2*se,
             ymax = mean + 2*se)) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  geom_hline(yintercept=5e-3, linetype=2, color = "red") +
  labs(x = "Dimension of covariate Z",
       y = "Type I Error") +
  theme(legend.position = "bottom") 

ggsave(filename = "type_I_err_5e-3.png",
       plot = type_I_err_3,
       device = "png",
       width = 3,
       height = 3)

# create histograms of null p-values
histograms <- pvals |>
  melt(value.name = "p-value") |>
  ggplot(aes(x = `p-value`, fill = method)) +
  geom_histogram(binwidth = 0.1, boundary = 0, color = "black") +
  facet_wrap(~method) +
  scale_x_continuous(limits = c(0,1)) +
  theme(legend.position = "none")

ggsave(filename = "histograms.png",
       plot = histograms,
       device = "png",
       width = 7,
       height = 3)

# create plot of computation times
comp_times <- computation_times |>
  melt(value.name = "time") |>
  group_by(method) |>
  summarise(mean_time = mean(time)) |>
  ggplot(aes(x = method, color = method, y = mean_time)) +
  geom_point() +
  scale_y_continuous(trans = "log10",
                     limits = c(0.001, 10),
                     breaks = c(0.001, 0.01, 0.1, 1, 10),
                     labels = c(0.001, 0.01, 0.1, 1, 10)) +
  labs(x = NULL,
       y = "Time (seconds)") +
  theme(legend.position = "none")

ggsave(filename = "comp-times.png",
       plot = comp_times,
       device = "png",
       width = 3,
       height = 3)

# plot GCM sampling distribution

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


test_stats <- tibble(result) |>
  ggplot(aes(x = sampled_stats)) +
  scale_x_continuous(limits = c(-4.5, 4.5), 
                     breaks = seq(-4.5, 4.5, by = 1.5), 
                     labels = seq(-4.5, 4.5, by = 1.5)) +
  geom_histogram(mapping = aes(x = sampled_stats, y=..density..), color = "black", bins = 20) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red") +
  labs(x = "Sampling distribution of GCM test statistics")


ecdf_comparison <- tibble(result) |>
  mutate(normal_stats = rnorm(reps_calibration)) |>
  pivot_longer(cols = c("sampled_stats", "resampled_stats", "normal_stats"),
               values_to = "test_statistics",
               names_to = "method") |>
  ggplot(aes_string(x = "test_statistics", colour = "method")) + 
  stat_ecdf() +
  labs(x = "Value of test statistics", 
       y = "Empirical CDF") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  scale_color_discrete(name="")


ggsave(filename = "sampled-test-stats.pdf",
       plot = test_stats,
       device = "pdf",
       width = 3.5,
       height = 3.5)

ggsave(filename = "ecdf_comparison.pdf",
       plot = ecdf_comparison,
       device = "pdf",
       width = 5,
       height = 5)
