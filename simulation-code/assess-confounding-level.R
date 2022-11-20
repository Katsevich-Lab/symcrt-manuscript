# define type-I error computation function
typeI_err <- function(a, alpha){
  1 - pnorm(qnorm(1-alpha/2) - a) + 
    pnorm(-qnorm(1-alpha/2) - a)
}

# Li and Liu (2022) setting
if(!file.exists("simulation-results/confounding/Li.csv")){
  cat("Assessing confounding level of Li et al (2022)...\n")
  MC <- 100
  confouding_level_Li <- numeric(MC)
  set.seed(1)
  for (r in 1:MC) {
    set.seed(r)
    n <- 250
    p <- 500
    beta <- numeric(p)
    magnitude <- 0.3
    beta[1:5] <- magnitude*(2*rbinom(5, 1, 0.5) - 1)
    gamma <- beta
    rho <- 0.5
    sig_Z <- katlabutils::generate_cov_ar1(rho = rho, d = p)
    Z <- katlabutils::fast_generate_mvn(mean = numeric(p),
                                        covariance = sig_Z,
                                        num_samples = n)
    res_X_Z <- rnorm(n)
    res_Y_Z <- rnorm(n)
    X <- Z %*% gamma + res_X_Z
    Y <- Z %*% beta + res_Y_Z
    c <- symcrt::simulate_confounding(n, X, Y)
    confouding_level_Li[r] <- c
  }
  saveRDS(mean(confouding_level_Li), "simulation-results/confounding/Li.rds")
} else{
  cat("Li (2022) confounding level already computed!\n")
}

# Liu et al (2021) adjacent setting

if(!file.exists("simulation-results/confounding/Liu.csv")){
  cat("Assessing confounding level of Liu et al (2022)...\n")
  set.seed(1)
  n <- 800
  MC <- 100*n
  p <- 800
  s <- 50
  magnitude <- 0.175
  beta_s <- numeric(p)
  beta_eqsp <- numeric(p)
  
  # first s coordinates are nonzero
  beta_s[1:s] <- magnitude*(2*rbinom(s, 1, 0.5) - 1)
  
  # nonzero components are equally p/s spaced.
  nonzero_seq_L <- seq(1, p, p/s)
  zero_seq_L <- setdiff(1:p, nonzero_seq_L)
  beta_eqsp[nonzero_seq_L] <- magnitude*(2*rbinom(s, 1, 0.5) - 1)
  
  # generate covariate
  cat("Generating Monte Carlo samples from (X,Y,Z)...\n")
  rho <- 0.5
  sig_Z <- katlabutils::generate_cov_ar1(rho = rho, d = p)
  covariate <- katlabutils::fast_generate_mvn(mean = numeric(p),
                                              covariance = sig_Z,
                                              num_samples = MC)
  
  # compute the response in first setting (first s are nonzero)
  Y_s <- covariate%*%beta_s + rnorm(MC)
  
  # compute the response in equally spaced setting
  Y_eqsp <- covariate%*%beta_eqsp + rnorm(MC)
  
  
  # compute the confounding level for each test setting
  cat("Computing Type-I error for marginal GCM test...\n")
  confounding_level_s_Liu <- numeric(p)
  confounding_level_eqsp_Liu <- numeric(p)
  for (r in 1:p) {
    confounding_level_eqsp_Liu[r] <- symcrt::simulate_confounding(n =n,
                                                                  X = covariate[, r],
                                                                  Y = Y_eqsp)
    confounding_level_s_Liu[r] <- symcrt::simulate_confounding(n = n,
                                                               X = covariate[, r],
                                                               Y = Y_s)
  }
  
  first_s_Liu <- data.frame(confounding = confounding_level_s_Liu, index = c(1:(p)))
  eqsp_Liu <- data.frame(confounding = confounding_level_eqsp_Liu, index = c(1:(p)))
  
  # Liu et al setting
  eqsp_Liu$setting <- "Liu et al.: Equally Spaced Signals"
  first_s_Liu$setting <- "Liu et al.: Concentrated Signals"
  Liu <- rbind(eqsp_Liu, first_s_Liu)
  Liu$class <- "Null"
  Liu$class[which((Liu$index %in% nonzero_seq_L) & (Liu$setting == "Liu et al.: Equally Spaced Signals"))] <- "Alternative"
  Liu$class[which((Liu$index %in% c(1:s)) & (Liu$setting == "Liu et al.: Concentrated Signals"))] <- "Alternative"
  
  # compute type-I error for Liu et al.
  Liu$type_I_err <- typeI_err(Liu$confounding, 0.05)
  
  write_csv(Liu, "simulation-results/confounding/Liu.csv")
} else{
  cat("Liu et al (2022) confounding level already computed!\n")
}

# Candès et al. (2018) setting
if(!file.exists("simulation-results/confounding/Candès.csv")){
  cat("Assessing confounding level of Candes et al (2018)...\n")
  set.seed(1)
  n <- 800
  MC <- 100*n
  p <- 1500
  s <- 50
  magnitude <- 20
  beta_rd <- numeric(p)
  # beta_eqsp <- numeric(p)
  
  # first s coordinates are nonzero
  beta_rd[sample(1:p, s)] <- magnitude*(2*rbinom(s, 1, 0.5) - 1)
  
  # nonzero components are equally p/s spaced.
  nonzero_seq_C <- which(beta_rd != 0)
  zero_seq_C <- setdiff(1:p, nonzero_seq_C)
  beta_rd[nonzero_seq_C] <- magnitude*(2*rbinom(s, 1, 0.5) - 1)
  
  cat("Generating Monte Carlo samples from (X,Y,Z)...\n")
  # generate covariate
  rho <- 0.3
  sig_Z <- katlabutils::generate_cov_ar1(rho = rho, d = p)
  covariate <- katlabutils::fast_generate_mvn(mean = numeric(p),
                                              covariance = sig_Z,
                                              num_samples = MC)
  
  # scale the matrix as performed in the orginal paper
  covariate_scale <- covariate
  for (b in 1:100) {
    covariate_scale[((b-1)*n + 1):(b*n), ] <- scale(covariate[((b-1)*n + 1):(b*n), ]) / sqrt(n - 1)
  }
  
  # compute the response in first setting (first s are nonzero)
  predictor_rd <- covariate_scale%*%beta_rd
  Y_s <- rbinom(MC, 1, prob = exp(predictor_rd)/(1+exp(predictor_rd)))
  
  # compute the confounding level for each test setting
  confounding_level_rd_Candès <- numeric(p)
  cat("Computing Type-I error for marginal GCM test...\n")
  for (r in 1:p) {
    if(r %% 10 == 0){
      cat(sprintf("Working of variable %d out of %d...\n", r, p))
    }
    confounding_level_rd_Candès[r] <- symcrt::simulate_confounding(n = n, 
                                                                   X = covariate_scale[, r], 
                                                                   Y = Y_s)
  }
  
  rd_Candès <- data.frame(confounding = confounding_level_rd_Candès, index = c(1:(p)))
  
  
  # Candes et al setting
  rd_Candès$class <- "Null"
  rd_Candès$class[which(rd_Candès$index %in% nonzero_seq_C)] <- "Alternative"
  rd_Candès$setting <- "Candès et al.: Randomly Distributed Signals"
  
  # compute type-I error for Candes et al.
  rd_Candès$type_I_err <- typeI_err(rd_Candès$confounding, 0.05)
  
  # save the csv file
  write.csv(rd_Candès, "simulation-results/confounding/Candès.csv")
} else{
  cat("Candes et al (2018) confounding level already computed!\n")
}