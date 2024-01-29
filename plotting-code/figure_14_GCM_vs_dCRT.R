library(katlabutils)
library(MASS)
library(reshape2)
library(ggplot2)
library(tidyverse)

# load the result
pvals <- readRDS("simulation-results/GCM_vs_dCRT/pvals.rds")
GCM_result <- readRDS("simulation-results/GCM_vs_dCRT/GCM_test_stat.rds")

# create QQ plots of null p-values
my_theme <-   theme_bw() + 
  theme(
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

# filepaths
filepath_qq <- "qq-plots.pdf"
filepath_type_I_err <- "type_I_err_5e-2.pdf"
filepath_histogram <- "sampled-test-stats.pdf"
filepath_ecdf <- "ecdf_comparison.pdf"
folder_path <- "figures"

# file names
file_names <- c(filepath_qq, filepath_type_I_err, 
                filepath_histogram, filepath_ecdf)

# Construct the full paths
full_paths <- file.path(folder_path, file_names)

# Check if files exist
files_exist <- file.exists(full_paths)

if(all(files_exist)){
  # if figure already exists, don't recompute it
  cat("Figure 14 already exists!\n")
}else{
  # qq-plot 
  qq_plots <- pvals |>
    melt(value.name = "p-value") |>
    ggplot(aes(y = `p-value`, color = method)) +
    stat_qq_points() +
    stat_qq_band() +
    geom_abline() +
    scale_x_continuous(trans = revlog_trans(), 
                       breaks = c(1, 0.1, 0.01, 0.001)) +
    scale_y_continuous(trans = revlog_trans(), 
                       breaks = c(1, 1e-2, 1e-4, 1e-6, 1e-8)) +
    my_theme +
    theme(legend.position = c(0.2, 0.8),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14)) +
    labs(x = "Expected null p-value",
         y = "Observed p-value")
  
  ggsave(filename = sprintf("%s/%s", folder_path, filepath_qq),
         plot = qq_plots,
         device = "pdf",
         width = 5,
         height = 5)
  
  
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
    my_theme + 
    theme(legend.position = c(0.2, 0.8),
          legend.text = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16))
  
  ggsave(filename = sprintf("%s/%s", folder_path, filepath_type_I_err),
         plot = type_I_err_2,
         device = "pdf",
         width = 5,
         height = 5)
  
  # plot the histogram for GCM sampling test stat
  test_stats <- tibble(GCM_result) |>
    ggplot(aes(x = sampled_stats)) +
    scale_x_continuous(limits = c(-4, 4), 
                       breaks = seq(-4.5, 4.5, by = 1.5), 
                       labels = seq(-4.5, 4.5, by = 1.5)) +
    geom_histogram(mapping = aes(x = sampled_stats, y=..density..), color = "black", bins = 20) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red") +
    labs(x = "Sampling distribution of GCM test statistics") + 
    my_theme +
    theme(legend.position = c(0.2, 0.8),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16))
  
  ggsave(filename = sprintf("%s/%s", folder_path, filepath_histogram),
         plot = test_stats,
         device = "pdf",
         width = 5,
         height = 5)
  
  
  # plot the empirical CDF
  ecdf_comparison <- tibble(GCM_result) |>
    mutate(normal_stats = rnorm(reps_calibration)) |>
    pivot_longer(cols = c("sampled_stats", "resampled_stats", "normal_stats"),
                 values_to = "test_statistics",
                 names_to = "method") |>
    mutate(method = factor(method, 
                           levels = c("sampled_stats", "resampled_stats", "normal_stats"),
                           labels = c("GCM test statistics", "resampled test statistics", "normal realizations"))) |>
    ggplot(aes_string(x = "test_statistics", colour = "method")) + 
    stat_ecdf() +
    labs(x = "Value of test statistics", 
         y = "Empirical CDF") +
    scale_color_discrete(name="") +
    my_theme +
    theme(legend.position = c(0.32, 0.8),
          legend.text = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16))
  
  
  ggsave(filename = sprintf("%s/%s", folder_path, filepath_ecdf),
         plot = ecdf_comparison,
         device = "pdf",
         width = 5,
         height = 5)
  
}









