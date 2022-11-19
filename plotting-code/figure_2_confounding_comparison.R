######################################################################
#
# Plot Figure 2
#
######################################################################

# set up libraries and settings for plotting
source("plotting-code/plotting_setup.R")

# read the simulation results
Li <- read_csv("simulation-results/confounding/Li.csv")
Liu <- read_csv("simulation-results/confounding/Liu.csv")
rd_Candès <- read_csv("simulation-results/confounding/Candès.csv")

# combine data
simulation_data <- rbind(Liu, rd_Candès) %>%
  mutate(setting = factor(setting, 
                          levels = c("Candès et al.: Randomly Distributed Signals",
                                     "Liu et al.: Equally Spaced Signals",
                                     "Liu et al.: Concentrated Signals"),
                          labels = c("Candès et al. (2018)\nRandomly Distributed Signals",
                                     "Liu et al. (2022)\nEqually Spaced Signals",
                                     "Liu et al. (2022)\nConcentrated Signals")))

# define function to evaluate Type-I error
typeI_err <- function(a, alpha){
  1 - pnorm(qnorm(1-alpha/2) - a) + 
    pnorm(-qnorm(1-alpha/2) - a)
}

# dot + line plot
dotplot <- simulation_data %>%
  filter(class == "Null") %>%
  ggplot(aes(x = index, y = type_I_err)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(size = 1) + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = typeI_err(Li$x, 0.05), color = "dodgerblue") +
  ggh4x::facet_grid2(. ~ setting, scales = "free", independent = "x") +
  ylab("Type I error\n(marginal GCM test)") +
  xlab("Position of variable being tested") +
  geom_rug(data = simulation_data[which(simulation_data$class == "Alternative"), ],
           aes(x = index),
           inherit.aes = FALSE,
           colour = "darkgreen") +
  plotting_theme
  
# histogram plot
histplot <- simulation_data %>%
  filter(class == "Null") %>%
  ggplot(aes(x = type_I_err)) +
  geom_histogram(bins = 50, colour = "black") + 
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = typeI_err(Li$x, 0.05), colour = "dodgerblue") +
  facet_wrap(~setting, scales = "free_y") + 
  scale_x_continuous(limits = c(0,1)) +
  xlab("Type I error (marginal GCM test)") +
  ylab("Frequency") +
  plotting_theme +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# combine plots
plot <- plot_grid(dotplot, histplot, 
                  nrow = 2, 
                  align = "v")

# save the result
ggsave(plot = plot,
       filename = "figures/type_I_Err_inflation_comparison.pdf", 
       device = "pdf",
       width = TEXTWIDTH, 
       height = 0.5*TEXTHEIGHT)
