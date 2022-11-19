######################################################################
#
# Plot Figure 4.
#
######################################################################


source("plotting-code/plotting_setup.R")

# specify which setting we consider
distribution <- "gaussian"
way_to_learn <- "supervised"
setting <- "alternative"
source("plotting-code/Figure 3-4 (Type I error & power for Gaussian supervised)/help.R")

# create the method_list 
methods_df <- generate_method_list(method_strings = method_strings,
                                   distribution = distribution,
                                   way_to_learn = way_to_learn)


# extract the parameter_grid from the specifier object
simspec_dir <- sprintf(
  "simulation-code/sim_spec_objects/sim_spec_%s_%s_%s.rds",
  distribution,
  way_to_learn,
  setting
)

simspec_file <- readRDS(simspec_dir)
p_grid <- simspec_file@parameter_grid


# compute adaptive threshold
source("compute_threshold.R")

# change MaxwayCRT to Maxway CRT
methods_df$infer_method[which(methods_df$infer_method == "MaxwayCRT")] <- "Maxway CRT"
threshold$infer_method[which(threshold$infer_method == "MaxwayCRT")] <- "Maxway CRT"

# read the results from disk
simresults_dir <- sprintf(
  "simulation-results/%s_%s_%s_results.rds",
  distribution,
  way_to_learn,
  setting
)


results <- readRDS(simresults_dir)


# rename grid_row_id as grid_id
if("grid_row_id" %in% names(results)){
  results <- results |>
    dplyr::rename(grid_id = grid_row_id) |>
    # NOTE: grid_id is a factor of integers, but potentially in the wrong order.
    # To get the right order, convert to character before converting to integer.
    dplyr::mutate(grid_id = as.integer(as.character(grid_id)))
}


# create the factor matrix
variable_parameters = bind_rows(
  dplyr::tibble(n = 100 * 2^seq(0, 4, 1), d = 400, s = 5, rho = 0.4,  
                variable_setting = paste0("n = ", n),
                fixed_setting = sprintf("d = %d, s = %d, rho = %.1f", d, s, rho)),
  dplyr::tibble(n = 200, d = 100 * 2^seq(0, 4, 1), s = 5, rho = 0.4,
                variable_setting = paste0("d = ", d),
                fixed_setting = sprintf("n = %d, s = %d, rho = %.1f", n, s, rho)),
  dplyr::tibble(n = 200, d = 400, s = 5 * 2^seq(0, 4, 1), rho = 0.4,
                variable_setting = paste0("s = ", s),
                fixed_setting = sprintf("n = %d, d = %d, rho = %.1f", n, d, rho)),
  dplyr::tibble(n = 200, d = 400, s = 5, rho = seq(0, 0.8, 0.2),
                variable_setting = paste0("rho = ", rho),
                fixed_setting = sprintf("n = %d, d = %d, s = %d", n, d, s))) %>%
  dplyr::mutate(setting = row_number(), 
                variable_setting = factor(variable_setting),
                fixed_setting = factor(fixed_setting))


# plot purely the result for dCRT
result <- list()
varying_values <- colnames(p_grid)[1:4]
rejection <- list()
for (k in 1:length(varying_values)) {
  # find the varying index
  name <- varying_values[k]
  assign("name", name)
  var_name <- sprintf("arm_%s", name)
  col_index <- which(colnames(p_grid)==var_name)
  index_name <- which(colnames(p_grid)==name)
  grid_ind <- data.frame(grid_id = p_grid$grid_id[which(p_grid[,col_index]==TRUE)],
                         name_1 = p_grid[which(p_grid[,col_index]==TRUE),index_name],
                         name_2 = p_grid[which(p_grid[,col_index]==TRUE),]$theta)
  
  
  # rename grid_ind
  colnames(grid_ind) <- c("grid_id", name, "theta")
  # join results with p_grid and methods_df
  result[[k]] <- results[which(results$grid_id %in% grid_ind$grid_id),] |>
    dplyr::left_join(grid_ind, by = "grid_id") |>
    dplyr::left_join(methods_df, by = "method")
  
  # join results with methods_df and threshold
  result[[k]] <- result[[k]] |>
    left_join(threshold[which(threshold$grid_id %in% grid_ind$grid_id), ], 
              by = c("infer_method", "reg_method", "grid_id"))
  
  ### compute rejection rate
  rejection_rate <- result[[k]] |>
    filter(parameter == "test_statistic") |>
    group_by_at(c(name, "theta", "infer_method", "reg_method")) |>
    summarise(
      rejection_rate = mean(value < cutoff_lower, na.rm = TRUE) + 
        mean(value > cutoff_upper, na.rm = TRUE)
    ) |>
    dplyr::left_join(variable_parameters[(5*(k-1)+1):(5*k),], by = name) |>
    dplyr::mutate(infer_reg = sprintf("%s (%s)", infer_method, reg_method)) |>
    ungroup() |>
    dplyr::group_by_at(name) |>
    dplyr::mutate(theta_scale = theta/max(theta))
  
  # change Maxway CRT(LASSO) to Maxway CRT
  rejection_rate$infer_reg[which(rejection_rate$infer_reg == "Maxway CRT (LASSO)")] <- "Maxway CRT"
  
  # select only three name values
  grid_ind <- arrange(grid_ind, by = name)
  list_1 <- 1:5
  list_3 <- 11:15
  list_5 <- 21:25
  list <- c(list_1, list_3, list_5)
  
  # select the values
  rejection_rate <- rejection_rate[which(as.data.frame(rejection_rate)[,1] %in% grid_ind[list, 2]), ]
  
  
  # store the rejection rate
  rejection[[k]] <- rejection_rate
}


p1 =  rejection[[1]] |>
  dplyr::mutate(variable_setting = factor(variable_setting, levels=c("n = 100","n = 200","n = 400","n = 800","n = 1600"))) |>
  ggplot(aes_string(x = "theta_scale", 
                    y = "rejection_rate",
                    colour = "infer_reg",
                    linetype = "infer_reg")) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = c(0, 1), minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
  facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
  geom_point(size = 0.5) +
  geom_line() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dashed", "solid")) +
  scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2"))+
  theme_bw() + 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none")


p2 =  rejection[[2]] |>
  dplyr::mutate(variable_setting = factor(variable_setting, levels=c("d = 100","d = 200","d = 400","d = 800","d = 1600"))) |>
  ggplot(aes_string(x = "theta_scale", 
                    y = "rejection_rate",
                    colour = "infer_reg",
                    linetype = "infer_reg")) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
  facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
  geom_point(size = 0.5) +
  geom_line() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dashed", "solid")) +
  scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2"))+
  theme_bw() + 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none")

p3 =  rejection[[3]] |>
  dplyr::mutate(variable_setting = factor(variable_setting, levels=c("s = 5","s = 10","s = 20","s = 40","s = 80"))) |>
  ggplot(aes_string(x = "theta_scale", 
                    y = "rejection_rate",
                    colour = "infer_reg",
                    linetype = "infer_reg")) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
  facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
  geom_point(size = 0.5) +
  geom_line() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dashed", "solid")) +
  scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2"))+
  theme_bw() + 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none")

p4 =  rejection[[4]] |>
  dplyr::mutate(variable_setting = factor(variable_setting, levels=c("rho = 0","rho = 0.2","rho = 0.4","rho = 0.6","rho = 0.8"))) |>
  ggplot(aes_string(x = "theta_scale", 
                    y = "rejection_rate",
                    colour = "infer_reg",
                    linetype = "infer_reg")) + 
  scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
  facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
  geom_point(size = 0.5) +
  geom_line() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dashed", "solid")) +
  scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2"))+
  theme_bw() + 
  theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm")),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="none")


# get legend
auxiliary_1 <- rejection[[4]] |>
  dplyr::mutate(variable_setting = factor(variable_setting, levels=c("rho = 0","rho = 0.2","rho = 0.4","rho = 0.6","rho = 0.8"))) |>
  dplyr::filter(reg_method != "oracle") |>
  ggplot(aes_string(x = "theta", 
                    y = "rejection_rate",
                    colour = "infer_reg",
                    linetype = "infer_reg")) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  ggh4x::facet_grid2(variable_setting ~ fixed_setting, scales = "free", independent = "x") +
  geom_point(size = 0.5) +
  geom_line() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dashed", "solid")) +
  scale_color_manual(values = c(rep("red", 2), "#90ee90", "#90ee90", "#129cf2"))+
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm')) + 
  theme(legend.title=element_blank()) 

auxiliary_2 <- rejection[[4]] |>
  dplyr::mutate(variable_setting = factor(variable_setting, levels=c("rho = 0","rho = 0.2","rho = 0.4","rho = 0.6","rho = 0.8"))) |>
  dplyr::filter(reg_method == "oracle") |>
  ggplot(aes_string(x = "theta", 
                    y = "rejection_rate",
                    colour = "infer_reg",
                    linetype = "infer_reg")) + 
  scale_x_continuous(expand = c(0.01,0.01)) +
  ggh4x::facet_grid2(variable_setting ~ fixed_setting, scales = "free", independent = "x") +
  geom_point(size = 0.5) +
  geom_line() +
  scale_linetype_manual(values = c( "dotted")) +
  scale_color_manual(values = c("grey"))+
  theme_bw() + 
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.4, 'cm')) + 
  theme(legend.title=element_blank()) 


legend_1 <- get_legend(auxiliary_1)
legend_2 <- get_legend(auxiliary_2)


# combine plots
plot <- plot_grid(p1, p2, p3, p4, ncol = 4, align = "v")


# create common x and y labels
y.grob <- textGrob("Power", 
                   gp=gpar(col="black"), rot=90)

x.grob <- textGrob("Effect size", 
                   gp=gpar(col="black"))

# add to plot
g <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow=1), 
                  plot_grid(legend_1, legend_2, nrow = 2), nrow=2,heights=c(7, 1))


# save the plot
ggsave(plot = g,
       filename = sprintf("figures/%s_%s_setting_partial_power.pdf", distribution, way_to_learn), 
       device = "pdf",
       width = TEXTWIDTH, 
       height = 0.6*TEXTHEIGHT)
