######################################################################
#
# Plot the results of Figure 1
#
######################################################################

# set up libraries and settings for plotting
source("plotting-code/plotting_setup.R")

# specify which setting we consider
sim_version <- "v1"
distribution <- "gaussian"
way_to_learn <- "supervised"
setting <- "null"

# get parameters from this simulation version
source("simulation-code/sim_versions/sim_v1.R")

# create the list of methods applied in this simulation
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
    test_hyperparams$normalize <- "normalized"
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
    X_on_Z_reg = X_on_Z_reg$mean_method_type,
    Y_on_Z_reg = Y_on_Z_reg$mean_method_type,
    normalization = test_hyperparams$normalize
  ) |>
  select(method, infer_method, X_on_Z_reg, Y_on_Z_reg, normalization)

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
    filter(infer_method == "dCRT") |>
    filter(X_on_Z_reg %in% c("LASSO", "PLASSO")) |>
    filter(Y_on_Z_reg %in% c("naive", "LASSO", "PLASSO")) |>
    filter(normalization != "unnormalized")

  # compute Type-I errors and standard errors
  level <- 0.05
  type_I_error <- result[[k]] |>
    filter(parameter == "p_value") |>
    group_by_at(c(name, "nu", "normalization", "X_on_Z_reg", "Y_on_Z_reg")) |>
    summarise(
      type_I_err = mean(value <= level, na.rm = TRUE),
      type_I_err_se = sd(value <= level, na.rm = TRUE) / sqrt(n())
    ) |>
    left_join(variable_parameters[(5 * (k - 1) + 1):(5 * k), ], by = name) |>
    ungroup()

  # store the type_I_error
  type_I_err[[k]] <- type_I_error
}

## extract the case when n=1600
type_I_err[[1]]$Y_on_Z_reg[which(type_I_err[[1]]$Y_on_Z_reg == "naive")] <-
  "Intercept-only"

# plot with ggplot
pn1 <- type_I_err[[1]] |>
  filter(n == 1600) |>
  filter(X_on_Z_reg %in% c("LASSO")) |>
  filter(Y_on_Z_reg %in% c("LASSO", "Intercept-only")) |>
  ggplot(aes_string(
    x = "nu",
    y = "type_I_err",
    linetype = "Y_on_Z_reg"
  )) +
  geom_point() +
  geom_line(aes_string(linetype = "Y_on_Z_reg")) +
  scale_y_continuous(limits = c(0, 0.42)) +
  scale_x_continuous(expand = c(0.005, 0.005)) +
  geom_errorbar(aes(
    ymin = pmax(0, type_I_err - 2 * type_I_err_se),
    ymax = pmin(1, type_I_err + 2 * type_I_err_se)
  ), width = 0.02) +
  geom_hline(yintercept = level, linetype = "dashed") +
  labs(
    x = TeX("Scaled norm of $\\beta$"),
    y = "Type I error",
    linetype = "Estimate of E[Y|Z]"
  ) +
  plotting_theme +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 12)
  )

# save the result
ggsave(
  plot = pn1,
  filename = "figures/negative_result_dCRT_Gaussian_supervised.pdf",
  device = "pdf",
  width = 0.62 * TEXTWIDTH,
  height = 0.6 * 0.7 * TEXTHEIGHT
)