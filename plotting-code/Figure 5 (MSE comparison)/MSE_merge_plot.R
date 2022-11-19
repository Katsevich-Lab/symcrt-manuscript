######################################################################
#
# Plot Figure 5.
#
######################################################################
library(ggplot2)
library(dplyr)
library(tibble)
library(MASS)
library(cowplot)
library(grid)
library(gridExtra)

TEXTWIDTH = 6.3
TEXTHEIGHT = 8.64

# specify which setting we consider
sim_version <- "diagnostic"
distribution <- "gaussian"
way_to_learn <- "supervised"

# source two files
source("NULL_MSE.R")
source("Alternative_MSE.R")


# summarise MSE according to group
MSE_null <- result_null |>
  filter(parameter %in% c("MSE_shared_X_Z", "MSE_shared_Y_Z")) |>
  group_by_at(c(name, "nu", "reg_method", "parameter", "test_type")) |>
  summarise(
    MSE = mean(value, na.rm = TRUE),
    MSE_se = sd(value, na.rm = TRUE) / sqrt(n())
  ) |>
  ungroup() |>
  dplyr::group_by_at(name) |>
  dplyr::mutate(nu_scale = nu/max(nu))

MSE_alternative <- result_alternative |>
  filter(parameter %in% c("MSE_total_X_Z", "MSE_total_Y_Z")) |>
  group_by_at(c(name, "theta", "reg_method", "parameter", "test_type")) |>
  summarise(
    MSE = mean(value, na.rm = TRUE),
    MSE_se = sd(value, na.rm = TRUE) / sqrt(n())
  ) |>
  ungroup() |>
  dplyr::group_by_at(name) |>
  dplyr::mutate(theta_scale = theta/max(theta))


MSE_null <- MSE_null |> as.data.frame()
MSE_null$parameter <- as.character(MSE_null$parameter)
MSE_alternative <- MSE_alternative |> as.data.frame()
MSE_alternative$parameter <- as.character(MSE_alternative$parameter)

MSE_null$parameter[which(MSE_null$parameter == "MSE_shared_X_Z")] <- "X on Z"
MSE_null$parameter[which(MSE_null$parameter == "MSE_shared_Y_Z")] <- "Y on Z"
MSE_alternative$parameter[which(MSE_alternative$parameter == "MSE_total_X_Z")] <- "X on Z"
MSE_alternative$parameter[which(MSE_alternative$parameter == "MSE_total_Y_Z")] <- "Y on Z"




# create MSE plot
p1 =  MSE_null |>
  ggplot(aes_string(x = "nu_scale", 
                    y = "MSE",
                    colour = "reg_method")) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = c(0,0.2,0.8), minor_breaks = seq(0, 0.8, 0.1), limits = c(0, NA)) +
  ggh4x::facet_grid2(parameter ~ test_type, 
                     scales = "fixed") +
  geom_point(size = 0.5) +
  geom_line() +
  scale_color_manual(values = c("red", "#90ee90"))+
  theme_bw() + 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm"), size = 10),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none")+
  xlab("Confounding strength") 


p2 =  MSE_alternative |>
  ggplot(aes_string(x = "theta_scale", 
                    y = "MSE",
                    colour = "reg_method")) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = c(0, 0.6,0.8), minor_breaks = seq(0, 0.8, 0.2), limits = c(0, NA)) +
  ggh4x::facet_grid2(parameter ~ test_type, 
                     scales = "fixed") +
  geom_point(size = 0.5) +
  geom_line() +
  scale_color_manual(values = c("red", "#90ee90"))+
  theme_bw() + 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
        strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm"), size = 10),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none")+
  xlab("Effect size") 

# get the legend
auxiliary =  MSE_null |>
  ggplot(aes_string(x = "nu_scale", 
                    y = "MSE",
                    colour = "reg_method")) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  ggh4x::facet_grid2(parameter ~ test_type, 
                     scales = "free", 
                     independent = "x") +
  geom_point(size = 0.5) +
  geom_line() +
  scale_color_manual(values = c("red", "#90ee90"))+
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none")+
  xlab("Confounding strength")+
  theme(axis.title.y = element_blank(),
       legend.position = "bottom") + 
  theme(legend.title=element_blank())

legend <- get_legend(auxiliary)

# combine plot
plot <- plot_grid(p1, p2, ncol = 2, align = "v")


# create common x and y labels
y.grob <- textGrob("MSE", 
                   gp=gpar(col="black"), rot=90)


# add to plot
g <- grid.arrange(arrangeGrob(plot, left = y.grob, nrow=1), 
                  plot_grid(legend), nrow=2, heights=c(8, 1))


# save the plot
ggsave(plot = g,
       filename = "MSE_merge_plog_n100.pdf", 
       device = "pdf",
       width = TEXTWIDTH, 
       height = 0.5*TEXTHEIGHT)


