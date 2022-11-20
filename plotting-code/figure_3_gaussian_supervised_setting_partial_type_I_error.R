######################################################################
#
# Plot Figure 3.
#
######################################################################

# specify which setting we consider
sim_version <- "benchmarking"
distribution <- "gaussian"
way_to_learn <- "supervised"
setting <- "null"
source("simulation-code/sim_versions/sim_benchmarking.R")

figure_path <- sprintf("figures/%s_%s_setting_null_partial.pdf", distribution, way_to_learn)
if(file.exists(figure_path)){
  # if figure already exists, don't recompute it
  cat("Figure 3 already exists!\n")
} else{
  cat("Generating Figure 3...\n")
    
  source("plotting-code/plotting_setup.R")
  
  # create the method_list
  methods_df <- symcrt::generate_method_df_from_strings(
    method_strings = method_strings[[setting]],
    distribution = distribution,
    way_to_learn = way_to_learn
  )
  
  
  # label the normalization type for each setting
  for (i in 1:nrow(methods_df)) {
    test_hyperparams <- methods_df$test_hyperparams[[i]]
    if (methods_df$test_type[i] == "GCM") {
      test_hyperparams$normalize <- "self_normalized"
    } else if (methods_df$test_type[i] == "MaxwayCRT") {
      test_hyperparams$normalize <- "unnormalized"
    } else if (test_hyperparams$normalize) {
      test_hyperparams$normalize <- "normalized"
    } else if (!test_hyperparams$normalize) {
      test_hyperparams$normalize <- "unnormalized"
    }
    methods_df$test_hyperparams[[i]] <- test_hyperparams
  }
  
  # augment methods_df with reg_method, lambda (for now ignore lambda) , and method
  methods_df <- methods_df |>
    mutate(method_idx = row_number()) |>
    rowwise() |>
    mutate(
      method = sprintf(
        "%s_%s_%s_%s_%d",
        test_type,
        test_hyperparams$way_to_learn,
        X_on_Z_reg$mean_method_type,
        Y_on_Z_reg$mean_method_type,
        method_idx
      ),
      infer_method = test_type,
      reg_method = X_on_Z_reg$mean_method_type,
      X_on_Z_reg = X_on_Z_reg$mean_method_type,
      Y_on_Z_reg = Y_on_Z_reg$mean_method_type,
      normalization = test_hyperparams$normalize
    ) |>
    select(method, infer_method, reg_method, X_on_Z_reg, Y_on_Z_reg, normalization)
  
  # extract the parameter_grid from the specifier object
  simspec_dir <- sprintf(
    "simulation-code/sim_spec_objects/%s/sim_spec_%s_%s_%s.rds",
    sim_version,
    distribution,
    way_to_learn,
    setting
  )
  
  simspec_file <- readRDS(simspec_dir)
  p_grid <- simspec_file@parameter_grid
  
  # change MaxwayCRT to Maxway CRT
  methods_df$infer_method[which(methods_df$infer_method == "MaxwayCRT")] <- "Maxway CRT"
  
  
  
  # read the results from disk
  simresults_dir <- sprintf(
    "simulation-results/%s/%s_%s_%s_results.rds",
    sim_version,
    distribution,
    way_to_learn,
    setting
  )
  
  
  results <- readRDS(simresults_dir)
  
  
  # rename grid_row_id as grid_id
  if ("grid_row_id" %in% names(results)) {
    results <- results |>
      rename(grid_id = grid_row_id) |>
      # NOTE: grid_id is a factor of integers, but potentially in the wrong order.
      # To get the right order, convert to character before converting to integer.
      mutate(grid_id = as.integer(as.character(grid_id)))
  }
  
  
  # create the factor matrix
  variable_parameters <- bind_rows(
    tibble(
      n = 100 * 2^seq(0, 4, 1), d = 400, s = 5, rho = 0.4,
      variable_setting = paste0("n = ", n),
      fixed_setting = sprintf("d = %d, s = %d, rho = %.1f", d, s, rho)
    ),
    tibble(
      n = 200, d = 100 * 2^seq(0, 4, 1), s = 5, rho = 0.4,
      variable_setting = paste0("d = ", d),
      fixed_setting = sprintf("n = %d, s = %d, rho = %.1f", n, s, rho)
    ),
    tibble(
      n = 200, d = 400, s = 5 * 2^seq(0, 4, 1), rho = 0.4,
      variable_setting = paste0("s = ", s),
      fixed_setting = sprintf("n = %d, d = %d, rho = %.1f", n, d, rho)
    ),
    tibble(
      n = 200, d = 400, s = 5, rho = seq(0, 0.8, 0.2),
      variable_setting = paste0("rho = ", rho),
      fixed_setting = sprintf("n = %d, d = %d, s = %d", n, d, s)
    )
  ) %>%
    mutate(
      setting = row_number(),
      variable_setting = factor(variable_setting),
      fixed_setting = factor(fixed_setting)
    )
  
  
  # plot purely the result for dCRT
  result <- list()
  varying_values <- colnames(p_grid)[1:4]
  type_I_err <- list()
  for (k in 1:length(varying_values)) {
    # find the varying index
    name <- varying_values[k]
    assign("name", name)
    var_name <- sprintf("arm_%s", name)
    col_index <- which(colnames(p_grid) == var_name)
    index_name <- which(colnames(p_grid) == name)
    grid_ind <- data.frame(
      grid_id = p_grid$grid_id[which(p_grid[, col_index] == TRUE)],
      name_1 = p_grid[which(p_grid[, col_index] == TRUE), index_name],
      name_2 = p_grid[which(p_grid[, col_index] == TRUE), ]$nu
    )
  
  
    # rename grid_ind
    colnames(grid_ind) <- c("grid_id", name, "nu")
    # join results with p_grid and methods_df
    result[[k]] <- results[which(results$grid_id %in% grid_ind$grid_id), ] |>
      left_join(grid_ind, by = "grid_id") |>
      left_join(methods_df, by = "method")
  
    # choose the result that is normalized
    result[[k]] <- result[[k]] |>
      filter(X_on_Z_reg == Y_on_Z_reg) |>
      filter(normalization != "normalized")
  
    # compute Type-I errors and standard errors
    level <- 0.05
    type_I_error <- result[[k]] |>
      filter(parameter == "p_value") |>
      group_by_at(c(name, "nu", "normalization", "reg_method", "X_on_Z_reg", "Y_on_Z_reg", "infer_method")) |>
      summarise(
        type_I_err = mean(value <= level, na.rm = TRUE),
        type_I_err_se = sd(value <= level, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) |>
      left_join(variable_parameters[(5 * (k - 1) + 1):(5 * k), ], by = name) |>
      mutate(infer_reg = sprintf("%s (%s)", infer_method, reg_method)) |>
      group_by_at(name) |>
      mutate(nu_scale = nu / max(nu)) |> 
      ungroup()
  
    # change Maxway CRT(LASSO) to Maxway CRT
    type_I_error$infer_reg[which(type_I_error$infer_reg == "Maxway CRT (LASSO)")] <- "Maxway CRT"
  
  
    # select only three name values
    grid_ind <- arrange(grid_ind, by = name)
    list_1 <- 1:5
    list_3 <- 11:15
    list_5 <- 21:25
    list <- c(list_1, list_3, list_5)
  
    # select the values
    type_I_error <- type_I_error[which(as.data.frame(type_I_error)[, 1] %in% grid_ind[list, 2]), ]
  
  
    # store the type_I_error
    type_I_err[[k]] <- type_I_error
  }
  
  ## plot with ggplot
  
  p1 <- type_I_err[[1]] |>
    mutate(variable_setting = factor(variable_setting, levels = c("n = 100", "n = 200", "n = 400", "n = 800", "n = 1600"))) |>
    ggplot(aes(x = nu_scale, y = type_I_err, colour = infer_reg, linetype = infer_reg)) +
    scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = c(0, 1), minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
    facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
    geom_point(size = 0.5) +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "solid", "dotted", "dashed", "solid")) +
    scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2")) +
    geom_hline(yintercept = level, linetype = "dashed") +
    theme_bw() +
    theme(
      strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    )
  
  
  p2 <- type_I_err[[2]] |>
    mutate(variable_setting = factor(variable_setting, levels = c("d = 100", "d = 200", "d = 400", "d = 800", "d = 1600"))) |>
    ggplot(aes(
      x = nu_scale,
      y = type_I_err,
      colour = infer_reg,
      linetype = infer_reg
    )) +
    scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
    facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
    geom_point(size = 0.5) +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "solid", "dotted", "dashed", "solid")) +
    scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2")) +
    geom_hline(yintercept = level, linetype = "dashed") +
    theme_bw() +
    theme(
      strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none"
    )
  
  p3 <- type_I_err[[3]] |>
    mutate(variable_setting = factor(variable_setting, levels = c("s = 5", "s = 10", "s = 20", "s = 40", "s = 80"))) |>
    ggplot(aes(
      x = nu_scale,
      y = type_I_err,
      colour = infer_reg,
      linetype = infer_reg
    )) +
    scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
    facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
    geom_point(size = 0.5) +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "solid", "dotted", "dashed", "solid")) +
    scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2")) +
    geom_hline(yintercept = level, linetype = "dashed") +
    theme_bw() +
    theme(
      strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none"
    )
  
  p4 <- type_I_err[[4]] |>
    mutate(variable_setting = factor(variable_setting, levels = c("rho = 0", "rho = 0.2", "rho = 0.4", "rho = 0.6", "rho = 0.8"))) |>
    ggplot(aes(
      x = nu_scale,
      y = type_I_err,
      colour = infer_reg,
      linetype = infer_reg
    )) +
    scale_x_continuous(minor_breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.25), limits = c(0.02, NA)) +
    facet_wrap(variable_setting ~ ., scales = "free_y", ncol = 1) +
    geom_point(size = 0.5) +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "solid", "dotted", "dashed", "solid")) +
    scale_color_manual(values = c(rep("red", 2), "#90ee90", "grey", "#90ee90", "#129cf2")) +
    geom_hline(yintercept = level, linetype = "dashed") +
    theme_bw() +
    theme(
      strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "none"
    )
  
  # get legend
  auxiliary_1 <- type_I_err[[4]] |>
    mutate(variable_setting = factor(variable_setting, levels = c("rho = 0", "rho = 0.2", "rho = 0.4", "rho = 0.6", "rho = 0.8"))) |>
    filter(reg_method != "naive") |>
    ggplot(aes(
      x = nu_scale,
      y = type_I_err,
      colour = infer_reg,
      linetype = infer_reg
    )) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    ggh4x::facet_grid2(variable_setting ~ fixed_setting, scales = "free", independent = "x") +
    geom_point(size = 0.5) +
    geom_line() +
    scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed", "solid")) +
    scale_color_manual(values = c(rep("red", 2), "#90ee90", "#90ee90", "#129cf2")) +
    geom_hline(yintercept = level, linetype = "dashed") +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm")
    ) +
    theme(legend.title = element_blank())
  
  auxiliary_2 <- type_I_err[[4]] |>
    mutate(variable_setting = factor(variable_setting, levels = c("rho = 0", "rho = 0.2", "rho = 0.4", "rho = 0.6", "rho = 0.8"))) |>
    filter(reg_method == "naive") |>
    ggplot(aes(
      x = nu_scale,
      y = type_I_err,
      colour = infer_reg,
      linetype = infer_reg
    )) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    ggh4x::facet_grid2(variable_setting ~ fixed_setting, scales = "free", independent = "x") +
    geom_point(size = 0.5) +
    geom_line() +
    scale_linetype_manual(values = c("dotted")) +
    scale_color_manual(values = c("grey")) +
    geom_hline(yintercept = level, linetype = "dashed") +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm")
    ) +
    theme(legend.title = element_blank())
  
  legend_1 <- get_legend(auxiliary_1)
  legend_2 <- get_legend(auxiliary_2)
  
  # combine the plot
  plot <- plot_grid(p1, p2, p3, p4, ncol = 4, align = "v")
  
  
  # create common x and y labels
  y.grob <- textGrob("Type I error",
    gp = gpar(col = "black"), rot = 90
  )
  
  x.grob <- textGrob("Confounding strength",
    gp = gpar(col = "black")
  )
  
  # add to plot
  g <- grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob, nrow = 1),
    plot_grid(legend_1, legend_2, nrow = 2),
    nrow = 2, heights = c(7, 1)
  )
  
  # save the plot
  ggsave(
    plot = g,
    filename = figure_path,
    device = "pdf",
    width = TEXTWIDTH,
    height = 0.6 * TEXTHEIGHT
  )
}