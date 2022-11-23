######################################################################
#
# Plot Figure 5.
#
######################################################################

figure_path <- "figures/MSE.pdf"
if(file.exists(figure_path)){
  # if figure already exists, don't recompute it
  cat("Figure 5 already exists!\n")
} else{
  cat("Generating Figure 5...\n")

  source("plotting-code/plotting_setup.R")
  
  # specify which setting we consider
  setting <- c("null", "alternative")
  sim_version <- "diagnostic"
  distribution <- "gaussian"
  way_to_learn <- "supervised"
  source("simulation-code/sim_versions/sim_diagnostic.R")
  
  # compute MSE result for null and alternative settings
  for (w in 1:length(setting)) {
    # create the method_list 
    methods_df <- symcrt::generate_method_df_from_strings(
      method_strings = method_strings[[setting[w]]],
      distribution = distribution,
      way_to_learn = way_to_learn
    )
    
    # extract the parameter_grid from the specifier object
    simspec_dir <- sprintf(
      "simulation-code/sim_spec_objects/%s/sim_spec_%s_%s_%s.rds",
      sim_version,
      distribution,
      way_to_learn,
      setting[w]
    )
    
    simspec_file <- readRDS(simspec_dir)
    p_grid <- simspec_file@parameter_grid
    
    # augment methods_df with reg_method, lambda (for now ignore lambda) , and method
    methods_df <- methods_df |>
      mutate(method_idx = row_number()) |>
      rowwise() |>
      mutate(
        reg_method = X_on_Z_reg$mean_method_type,
        method = sprintf(
          "%s_%s_%s_%s_%d",
          test_type,
          test_hyperparams$way_to_learn,
          X_on_Z_reg$mean_method_type,
          Y_on_Z_reg$mean_method_type,
          method_idx),
        infer_method = test_type,
        way_to_learn = test_hyperparams$way_to_learn
      ) |>
      select(reg_method, method, infer_method)
    
    
    
    # read the results from disk
    simresults_dir <- sprintf(
      "simulation-results/%s/%s_%s_%s_results.rds",
      sim_version,
      distribution,
      way_to_learn,
      setting[w]
    )
    
    results <- readRDS(simresults_dir)
    
    
    
    # rename grid_row_id as grid_id
    if("grid_row_id" %in% names(results)){
      results <- results |>
        rename(grid_id = grid_row_id) |>
        # NOTE: grid_id is a factor of integers, but potentially in the wrong order.
        # To get the right order, convert to character before converting to integer.
        mutate(grid_id = as.integer(as.character(grid_id)))
    }
    
    
    # plot for MSE of Y|Z
    varying_values <- colnames(p_grid)[1:4]
    # find the varying index
    k <- 1
    name <- varying_values[k]
    assign("name", name)
    var_name <- sprintf("arm_%s", name)
    index_name <- which(colnames(p_grid)==name)
    if (setting[w] == "null"){
      grid_ind <- data.frame(grid_id = p_grid$grid_id[which(p_grid$n == 100)],
                             name_1 = p_grid[which(p_grid$n == 100), index_name],
                             name_2 = p_grid[which(p_grid$n == 100),]$nu)
      colnames(grid_ind) <- c("grid_id", name, "nu")
      
      # join results with p_grid and methods_df
      result_null <- results[which(results$grid_id %in% grid_ind$grid_id),] |>
        left_join(grid_ind, by = "grid_id") |>
        left_join(methods_df, by = "method") |>
        mutate(test_type = sprintf("%s setting (shared MSE)", setting[w]))
    }else{
      grid_ind <- data.frame(grid_id = p_grid$grid_id[which(p_grid$n == 100)],
                             name_1 = p_grid[which(p_grid$n == 100), index_name],
                             name_2 = p_grid[which(p_grid$n == 100),]$theta)
      colnames(grid_ind) <- c("grid_id", name, "theta")
      
      # join results with p_grid and methods_df
      result_alternative <- results[which(results$grid_id %in% grid_ind$grid_id),] |>
        left_join(grid_ind, by = "grid_id") |>
        left_join(methods_df, by = "method") |>
        dplyr::mutate(test_type = sprintf("%s setting (total MSE)", setting[w]))
      
    }
  }
  
  
  # summarise MSE according to group
  MSE_null <- result_null |>
    filter(parameter %in% c("MSE_shared_X_Z", "MSE_shared_Y_Z")) |>
    group_by_at(c(name, "nu", "reg_method", "parameter", "test_type")) |>
    summarise(
      MSE = mean(value, na.rm = TRUE),
      MSE_se = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) |>
    group_by_at(name) |>
    mutate(nu_scale = nu/max(nu)) |>
    ungroup()
  
  MSE_alternative <- result_alternative |>
    filter(parameter %in% c("MSE_total_X_Z", "MSE_total_Y_Z")) |>
    group_by_at(c(name, "theta", "reg_method", "parameter", "test_type")) |>
    summarise(
      MSE = mean(value, na.rm = TRUE),
      MSE_se = sd(value, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) |>
    group_by_at(name) |>
    mutate(theta_scale = theta/max(theta)) |>
    ungroup()
  
  
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
    ggplot(aes(x = nu_scale, 
                      y = MSE,
                      colour = reg_method)) + 
    scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 0.4)) +
    # scale_y_continuous(breaks = c(0,0.2,0.8), minor_breaks = seq(0, 0.8, 0.1), limits = c(0, NA)) +
    facet_grid(test_type ~ parameter, scales = "fixed") +
    geom_point(size = 0.5) +
    geom_line() +
    scale_color_manual(values = c("red", "#90ee90"))+
    plotting_theme + 
    theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
          strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm"), size = 10),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
#          axis.ticks.y = element_blank(),
#          axis.text.y = element_blank(),
          legend.position="none")+
    xlab(TeX("Marginal association between $X$ and $Y$ ($\\nu$)")) +
    ylab("Mean-squared estimation error")
  
  p2 =  MSE_alternative |>
    ggplot(aes(x = theta_scale, 
                      y = MSE,
                      colour = reg_method)) + 
    scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(limits = c(0, 0.8)) +
    # scale_y_continuous(breaks = c(0, 0.6,0.8), minor_breaks = seq(0, 0.8, 0.2), limits = c(0, NA)) +
    facet_grid(test_type ~ parameter, scales = "fixed") +
    geom_point(size = 0.5) +
    geom_line() +
    scale_color_manual(values = c("red", "#90ee90"))+
    plotting_theme + 
    theme(strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "cm"), size = 10),
          strip.text.y = element_text(margin = margin(0,0.08,0,0.08, "cm"), size = 10),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          # axis.ticks.y = element_blank(),
          # axis.text.y = element_blank(),
          legend.position="none")+
    xlab(TeX("Effect size ($\\theta$)")) +
    ylab("Mean-squared estimation error")
  
  # get the legend
  auxiliary =  MSE_null |>
    ggplot(aes(x = nu_scale, 
                      y = MSE,
                      colour = reg_method)) + 
    scale_x_continuous(expand = c(0.01,0.01)) +
    ggh4x::facet_grid2(parameter ~ test_type, 
                       scales = "free", 
                       independent = "x") +
    geom_point(size = 0.5) +
    geom_line() +
    scale_color_manual(values = c("red", "#90ee90"))+
    plotting_theme + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position="none")+
    xlab("Confounding strength")+
    theme(axis.title.y = element_blank(),
         legend.position = "bottom",
         legend.title=element_blank())
  
  legend <- get_legend(auxiliary)
  
  # combine plot
  plot <- plot_grid(p1, p2, ncol = 1, align = "v")
  
  
  # create common x and y labels
  y.grob <- textGrob("Mean-squared estimation error", 
                     gp=gpar(col="black"), rot=90)
  
  
  # add to plot
  g <- grid.arrange(arrangeGrob(plot, left = y.grob, nrow=1), 
                    plot_grid(legend), nrow=2, heights=c(8, 1))
  
  
  # save the plot
  ggsave(plot = g,
         filename = figure_path, 
         device = "pdf",
         width = 0.8*TEXTWIDTH, 
         height = 0.65*TEXTHEIGHT)
}
