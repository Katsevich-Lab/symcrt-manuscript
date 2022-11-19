if(setting == "null"){
  method_strings <- c("GCM naive naive",
                      "GCM LASSO LASSO",
                      "GCM PLASSO PLASSO",
                      "dCRT LASSO naive normalized",
                      "dCRT LASSO LASSO normalized",
                      "dCRT PLASSO naive normalized",
                      "dCRT PLASSO PLASSO normalized",
                      "dCRT LASSO naive unnormalized",
                      "dCRT LASSO LASSO unnormalized",
                      "dCRT PLASSO naive unnormalized",
                      "dCRT PLASSO PLASSO unnormalized",
                      "MaxwayCRT LASSO LASSO")
} else if(setting %in% c("calibration", "alternative")){
  method_strings <- c("GCM oracle oracle",
                      "GCM LASSO LASSO",
                      "GCM PLASSO PLASSO",
                      "dCRT LASSO LASSO unnormalized",
                      "dCRT PLASSO PLASSO unnormalized",
                      "MaxwayCRT LASSO LASSO")
} else{
  stop("setting must be one of null, calibration, alternative")
}


generate_method_list <- function(method_strings, distribution, way_to_learn) {
  set.seed(1)
  num_methods <- length(method_strings)
  
  # prepare lists for the four columns of the methods_df
  test_type_list <- vector("character", num_methods)
  X_on_Z_reg_list <- vector("list", num_methods)
  Y_on_Z_reg_list <- vector("list", num_methods)
  test_hyperparams_list <- vector("list", num_methods)
  
  # iterate over methods
  for(method_idx in 1:num_methods){
    method_string <- method_strings[method_idx]
    
    # split the string into pieces
    method_string_split <- strsplit(method_string, split = " ")[[1]]
    test_type <- method_string_split[1]
    
    # add regression methods
    if(test_type %in% c("GCM", "dCRT", "MaxwayCRT")){
      X_on_Z_method_type <- method_string_split[2]
      X_on_Z_reg <- list(mean_method_type = X_on_Z_method_type,
                         mean_method_hyperparams = list(family = distribution))
      Y_on_Z_method_type <- method_string_split[3]
      Y_on_Z_reg <- list(mean_method_type = Y_on_Z_method_type,
                         mean_method_hyperparams = list(family = distribution))
    }
    
    # add way_to_learn
    test_hyperparams <- list(way_to_learn = way_to_learn)
    
    # add the resample family
    if(test_type %in% c("dCRT", "MaxwayCRT")){
      test_hyperparams$resample_family <- distribution
    }
    
    # add the normalization
    if(test_type == "dCRT"){
      normalization <- method_string_split[4]
      if(normalization == "normalized"){
        test_hyperparams$normalize <- TRUE
      } else if(normalization == "unnormalized"){
        test_hyperparams$normalize <- FALSE
      } else{
        stop("normalization variable should either be normalized or unnormalized")
      }
    }
    
    # add the MSE option only for gaussian response with GCM statistic using LASSO/PLASSO
    test_hyperparams$MSE <- FALSE
    if(test_type == "GCM"){
      if(distribution == "gaussian"){
        if(X_on_Z_reg$mean_method_type %in% c("LASSO", "PLASSO")){
          if(Y_on_Z_reg$mean_method_type %in% c("LASSO", "PLASSO")){
            test_hyperparams$MSE <- TRUE
          }
        }
      }
    }
    
    
    # update the lists
    test_type_list[method_idx] <- test_type
    X_on_Z_reg_list[[method_idx]] <- X_on_Z_reg
    Y_on_Z_reg_list[[method_idx]] <- Y_on_Z_reg
    test_hyperparams_list[[method_idx]] <-test_hyperparams
  }
  
  # assemble the lists into a tibble, pass to generate_method_list_from_df, return
  tibble::tibble(test_type = test_type_list,
                 X_on_Z_reg = X_on_Z_reg_list,
                 Y_on_Z_reg = Y_on_Z_reg_list,
                 test_hyperparams = test_hyperparams_list) 
}