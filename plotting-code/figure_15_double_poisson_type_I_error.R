library(ggplot)
library(dplyr)

# check if the filepath is there or not
figure_path <- "figures/asymmetry_Poisson_type_I_error.pdf"

if(file.exists(figure_path)){
  # if figure already exists, don't recompute it
  cat("Figure 15 already exists!\n")
} else{
  # plot type-I error
  check_results <- readRDS("simulation-results/check_results_double_poisson.rds")
  type_I_err <- check_results$metrics |>
    mutate(method = factor(method, 
                           levels = c("dCRT_X", "dCRT_Y", "dCRT_X_normalize", "dCRT_Y_normalize", "gcm_test"),
                           labels = c("dCRT (resample X|Z)", "dCRT (resample Y|Z)", "ndCRT (resample X|Z)", "ndCRT (resample Y|Z)", "GCM"))) |>
    ggplot(aes(x = d,
               y = mean,
               ymin = mean - 2*se,
               ymax = mean + 2*se,
               color = method)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("dCRT (resample X|Z)"="#FF0033", 
                                  "dCRT (resample Y|Z)"="blue"))+
    geom_errorbar(width = 0.5) +
    geom_hline(yintercept=0.05, linetype=2, color = "black") +
    labs(x = "Dimension of covariate Z",
         y = "Type I Error") +
    theme_bw() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5),
          legend.key.size = unit(0.8, "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 10),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10)) +
    guides(color = guide_legend(nrow = 1))
  
  # save the plot
  ggsave(filename = figure_path,
         plot = type_I_err,
         device = "pdf",
         width = 6,
         height = 4)
}



