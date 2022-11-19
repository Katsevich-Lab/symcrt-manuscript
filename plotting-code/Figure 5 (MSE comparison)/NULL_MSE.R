setting <- "null"
source("help.R")

# create the method_list 
methods_df <- generate_method_list(method_strings = method_strings,
                                   distribution = distribution,
                                   way_to_learn = way_to_learn)

# extract the parameter_grid from the specifier object
simspec_dir <- sprintf("%s/private/spec_objects/%s/sim_spec_%s_%s_%s.rds",
                       .get_config_path("LOCAL_SYMCRT_DATA_DIR"),
                       sim_version,
                       distribution,
                       way_to_learn,
                       setting)

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
grid_ind <- data.frame(grid_id = p_grid$grid_id[which(p_grid$n == 100)],
                       name_1 = p_grid[which(p_grid$n == 100), index_name],
                       name_2 = p_grid[which(p_grid$n == 100),]$nu)
colnames(grid_ind) <- c("grid_id", name, "nu")

# join results with p_grid and methods_df
result_null <- results[which(results$grid_id %in% grid_ind$grid_id),] |>
  left_join(grid_ind, by = "grid_id") |>
  left_join(methods_df, by = "method") |>
  mutate(test_type = sprintf("%s setting (shared MSE)", setting))
