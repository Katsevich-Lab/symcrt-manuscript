######################################################################
#
# Plot the results of Figure 1
#
######################################################################
## 0.4233418 cm for font size
## we use 12 pt in main text. It is about 0.16667 inches.

library(ggplot2)
library(dplyr)
library(tibble)
library(grid)
library(latex2exp)


TEXTWIDTH = 6.3
TEXTHEIGHT = 8.64

# specify which setting we consider
sim_version <- "v1"
distribution <- "gaussian"
way_to_learn <- "supervised"
setting <- "null"
source("help.R")

# create the method_list 
methods_df <- generate_method_list(method_strings = method_strings,
                                   distribution = distribution,
                                   way_to_learn = way_to_learn)

# label the normalization type for each setting
for (i in 1:nrow(methods_df)) {
  test_hyperparams <- methods_df$test_hyperparams[[i]]
  if(methods_df$test_type[i] == "GCM"){
    test_hyperparams$normalize <- "self_normalized"
  }else if(methods_df$test_type[i] == "MaxwayCRT"){
    test_hyperparams$normalize <- "normalized"
  }else if(test_hyperparams$normalize){
    test_hyperparams$normalize <- "normalized"
  }else if(!test_hyperparams$normalize){
    test_hyperparams$normalize <- "unnormalized"
  }
  methods_df$test_hyperparams[[i]] <- test_hyperparams
}

# augment methods_df with reg_method, lambda (for now ignore lambda) , and method
methods_df <- methods_df |>
  dplyr::mutate(method_idx = row_number()) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    method = sprintf(
      "%s_%s_%s_%s_%d",
      test_type,
      test_hyperparams$way_to_learn,
      X_on_Z_reg$mean_method_type,
      Y_on_Z_reg$mean_method_type,
      method_idx),
    infer_method = test_type,
    X_on_Z_reg = X_on_Z_reg$mean_method_type,
    Y_on_Z_reg = Y_on_Z_reg$mean_method_type,
    normalization = test_hyperparams$normalize
  ) |>
  dplyr::select(method, infer_method, X_on_Z_reg, Y_on_Z_reg, normalization)

# extract the parameter_grid from the specifier object
simspec_dir <- sprintf("%s/private/spec_objects/%s/sim_spec_%s_%s_%s.rds",
                       .get_config_path("LOCAL_SYMCRT_DATA_DIR"),
                       sim_version,
                       distribution,
                       way_to_learn,
                       setting)

simspec_file <- readRDS(simspec_dir)
p_grid <- simspec_file@parameter_grid


# read the results from disk
simresults_dir <- sprintf("%s/private/results/%s/%s_%s_%s_results.rds",
                          .get_config_path("LOCAL_SYMCRT_DATA_DIR"),
                          sim_version,
                          distribution,
                          way_to_learn,
                          setting)

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
type_I_err <- list()
for (k in 1:length(varying_values)) {
  # find the varying index
  name <- varying_values[k]
  assign("name", name)
  var_name <- sprintf("arm_%s", name)
  col_index <- which(colnames(p_grid)==var_name)
  index_name <- which(colnames(p_grid)==name)
  grid_ind <- data.frame(grid_id = p_grid$grid_id[which(p_grid[,col_index]==TRUE)],
                         name_1 = p_grid[which(p_grid[,col_index]==TRUE),index_name],
                         name_2 = p_grid[which(p_grid[,col_index]==TRUE),]$nu)
  
  
  # rename grid_ind
  colnames(grid_ind) <- c("grid_id", name, "nu")
  # join results with p_grid and methods_df
  result[[k]] <- results[which(results$grid_id %in% grid_ind$grid_id),] |>
    dplyr::left_join(grid_ind, by = "grid_id") |>
    dplyr::left_join(methods_df, by = "method")
  
  # choose the result that is normalized
  result[[k]] <- result[[k]] |> 
    dplyr::filter(infer_method == "dCRT") |>
    dplyr::filter(X_on_Z_reg %in% c("LASSO", "PLASSO")) |>
    dplyr::filter(Y_on_Z_reg %in% c("naive", "LASSO", "PLASSO")) |>
    dplyr::filter(normalization != "unnormalized")
  
  # compute Type-I errors and standard errors
  level <- 0.05
  type_I_error <- result[[k]] |>
    dplyr::filter(parameter == "p_value") |>
    dplyr::group_by_at(c(name, "nu", "normalization", "X_on_Z_reg", "Y_on_Z_reg")) |>
    dplyr::summarise(
      type_I_err = mean(value <= level, na.rm = TRUE),
      type_I_err_se = sd(value <= level, na.rm = TRUE) / sqrt(n())
    ) |>
    dplyr::left_join(variable_parameters[(5*(k-1)+1):(5*k),], by = name) |>
    ungroup()
  
  # store the type_I_error
  type_I_err[[k]] <- type_I_error
}

## extract the case when n=1600
type_I_err[[1]]$Y_on_Z_reg[which(type_I_err[[1]]$Y_on_Z_reg == "naive")] <- 
  "Intercept-only"

# plot with ggplot
pn1 =  type_I_err[[1]] |>
  dplyr::filter(n == 1600) |>
  dplyr::filter(X_on_Z_reg %in% c("LASSO")) |>
  dplyr::filter(Y_on_Z_reg %in% c("LASSO", "Intercept-only")) |>
  ggplot(aes_string(x = "nu", 
                    y = "type_I_err",
                    linetype = "Y_on_Z_reg")) + 
  geom_point() +
  geom_line(aes_string(linetype = "Y_on_Z_reg")) +
  scale_y_continuous(limits = c(0, 0.42)) +
  scale_x_continuous(expand = c(0.005,0.005)) +
  geom_errorbar(aes(
    ymin = pmax(0, type_I_err - 2 * type_I_err_se),
    ymax = pmin(1, type_I_err + 2 * type_I_err_se)
  ), width = 0.02) +
  geom_hline(yintercept = level, linetype = "dashed") +
  labs(x = TeX("Scaled norm of $\\beta$"),
       y = "Type I error",
       linetype = "Estimate of E[Y|Z]"
  )+
  theme_bw() + 
  theme(legend.position = "bottom",
        axis.text=element_text(size=8),
        axis.title=element_text(size=12))

# save the result
ggsave(plot = pn1,
       filename = "negative_result_dCRT_Gaussian_supervised.pdf", 
       device = "pdf",
       width = 0.62*TEXTWIDTH, 
       height = 0.6*0.7*TEXTHEIGHT)


