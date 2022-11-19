library(dplyr)
library(tibble)
library(reshape)
library(MASS)

# choose a way to find the quantile: fit a Gaussian or use the empirical quantile
method <- "fit"

# read the results from disk
simsresults_dir <- sprintf("%s/private/results/%s/%s_%s_%s_results.rds",
                           .get_config_path("LOCAL_SYMCRT_DATA_DIR"),
                           sim_version,
                           distribution,
                           way_to_learn,
                           "calibration")

results <- readRDS(simsresults_dir)




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
  dplyr::select(reg_method, method, infer_method)

# rename grid_row_id as grid_id
if("grid_row_id" %in% names(results)){
  results <- results |>
    dplyr::rename(grid_id = grid_row_id) |>
    # NOTE: grid_id is a factor of integers, but potentially in the wrong order.
    # To get the right order, convert to character before converting to integer.
    dplyr::mutate(grid_id = as.integer(as.character(grid_id)))
}

# find the results that do not use naive fitting; 
# might later add no_nu as a parameter

result <- rep(list(), nrow(p_grid))
adaptive_threshold <- rep(list(), nrow(p_grid))


for (k in 1:nrow(p_grid)) {
  # extract the results from the corresponding result list
  result[[k]] <- results[which(results$grid_id == k), ] |> left_join(methods_df, by = "method")
  result[[k]] <- result[[k]][which(result[[k]]$parameter == "test_statistic"), ]
  
  # compute the adaptive threshold
  if(method == "fit"){
    set.seed(1)
    adaptive_threshold[[k]] <- result[[k]][complete.cases(result[[k]]), ] |> 
      group_by_at(c("infer_method", "reg_method")) |>
      summarise(
        cutoff_lower = qnorm(0.025, mean = fitdistr(value, "normal")$estimate[1], sd = fitdistr(value, "normal")$estimate[2]),
        cutoff_upper = qnorm(0.975, mean = fitdistr(value, "normal")$estimate[1], sd = fitdistr(value, "normal")$estimate[2])
      ) |>
      ungroup()
  }else{
    adaptive_threshold[[k]] <- result[[k]] |> 
      group_by_at(c("infer_method", "reg_method")) |>
      summarise(
        cutoff_lower = quantile(value, 0.025, na.rm = TRUE),
        cutoff_upper = quantile(value, 0.975, na.rm = TRUE)
      ) |>
      ungroup()
  }
  
  # add an extra line to adaptive_threshold[[k]]
  adaptive_threshold[[k]][nrow(adaptive_threshold[[k]])+1, 1] = "GCM"
  adaptive_threshold[[k]][nrow(adaptive_threshold[[k]]), 2] = as.character("oracle")
  adaptive_threshold[[k]][nrow(adaptive_threshold[[k]]), 3] = qnorm(0.025)
  adaptive_threshold[[k]][nrow(adaptive_threshold[[k]]), 4] = qnorm(0.975)
  
  # add a new column indexing the grid id
  adaptive_threshold[[k]]$grid_id <- k
}

threshold <- merge_all(adaptive_threshold)

