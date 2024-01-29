library(ggplot2)
library(dplyr)
library(tidyr)

# check if the filepath is there or not
figure_path <- "figures/histogram_asymmetry_investigation.pdf"
if(file.exists(figure_path)){
  # if figure already exists, don't recompute it
  cat("Figure 16 already exists!\n")
} else{
  # load the result
  result_teststat <- readRDS("simulation-results/asymmetry_investigation/result_teststat.rds")
  
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
  
  # histogram for resampling test stat versus sampling test stat
  histograms_unnormalized <- tibble(result_teststat) |>
    pivot_longer(cols = c("sampled_teststat", 
                          "resampled_X_teststat", 
                          "resampled_Y_teststat"
    ),
    values_to = "teststat",
    names_to = "category") |>
    mutate(category = case_when(
      category == "sampled_teststat" ~ "sampling test stat",
      category == "resampled_X_teststat" ~ "resampling test stat (X|Z)",
      category == "resampled_Y_teststat" ~ "resampling test stat (Y|Z)"
    )) |>
    ggplot(aes(x = teststat, fill = category)) +
    geom_histogram(bins = 20, boundary = 0, color = "black") +
    xlim(c(-40,40)) +
    facet_wrap(~category) +
    xlab("unnormalized test statistics") +
    my_theme +
    scale_fill_manual(values = c("resampling test stat (X|Z)"="#FF0033", 
                                 "resampling test stat (Y|Z)"="blue",
                                 "sampling test stat"="orange")) +
    theme(legend.position = "none",
          strip.text = element_text(size = 12),
          legend.text = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 16), 
          axis.title.y = element_text(size = 16))
  
  # save the histograms for resampling unnormalized test stat
  ggsave(figure_path,
         plot = histograms_unnormalized,
         width = 9,
         height = 5)
}
