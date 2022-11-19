######################################################################
#
# Plot Figure 2
#
######################################################################

library(ggplot2)
library(gridExtra)
library(readr)
library(dplyr)
library(cowplot)

TEXTWIDTH = 6.3
TEXTHEIGHT = 8.64

# define type-I error computation function
typeI_err <- function(a, alpha){
  1 - pnorm(qnorm(1-alpha/2) - a) + 
    pnorm(-qnorm(1-alpha/2) - a)
}


########## reproduce the data using following code ##########

# Li and Liu (2022) setting
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

print(mean(confouding_level_Li))
# 4.432801
print(sd(confouding_level_Li))
# 1.508544

write.csv(mean(confouding_level_Li), "Li.csv")


# Liu et al (2021) adjacent setting
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


# Candès et al. (2018) setting

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
# confounding_level_eqsp_Candès <- numeric(p)
for (r in 1:p) {
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
write.csv(rd_Candès, "Candès.csv")

# Liu et al setting
eqsp_Liu$setting <- "Liu et al.: Equally Spaced Signals"
first_s_Liu$setting <- "Liu et al.: Concentrated Signals"
Liu <- rbind(eqsp_Liu, first_s_Liu)
Liu$class <- "Null"
Liu$class[which((Liu$index %in% nonzero_seq_L) & (Liu$setting == "Liu et al.: Equally Spaced Signals"))] <- "Alternative"
Liu$class[which((Liu$index %in% c(1:s)) & (Liu$setting == "Liu et al.: Concentrated Signals"))] <- "Alternative"



# compute type-I error for Liu et al.
Liu$type_I_err <- typeI_err(Liu$confounding, 0.05)

write.csv(Liu, "Liu.csv")

########## reproduce the plot using saved data (or data generated above) ##########

Li <- read_csv("Li.csv")
Liu <- read_csv("Liu.csv")
rd_Candès <- read_csv("Candès.csv")

# combine data
simulation_data <- rbind(Liu, rd_Candès) %>%
  mutate(setting = factor(setting, 
                          levels = c("Candès et al.: Randomly Distributed Signals",
                                     "Liu et al.: Equally Spaced Signals",
                                     "Liu et al.: Concentrated Signals"),
                          labels = c("Candès et al. (2018)\nRandomly Distributed Signals",
                                     "Liu et al. (2022)\nEqually Spaced Signals",
                                     "Liu et al. (2022)\nConcentrated Signals")))

# dot + line plot

dotplot <- simulation_data %>%
  filter(class == "Null") %>%
  ggplot(aes(x = index, y = type_I_err)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_point(size = 1) + 
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  geom_hline(yintercept = typeI_err(Li$x, 0.05), color = "dodgerblue") +
  ggh4x::facet_grid2(. ~ setting, scales = "free", independent = "x") +
  ylab("Type I error\n(marginal GCM test)") +
  xlab("Position of variable being tested") +
  geom_rug(data = simulation_data[which(simulation_data$class == "Alternative"), ],
           aes(x = index),
           inherit.aes = FALSE,
           colour = "darkgreen") +
  theme_bw()
  
# histogram plot
histplot <- simulation_data %>%
  filter(class == "Null") %>%
  ggplot(aes(x = type_I_err)) +
  geom_histogram(bins = 50, colour = "black") + 
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = typeI_err(Li$x, 0.05), colour = "dodgerblue") +
  facet_wrap(~setting, scales = "free_y") + 
  scale_x_continuous(limits = c(0,1)) +
  xlab("Type I error (marginal GCM test)") +
  ylab("Frequency") +
  theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# combine plots
plot <- plot_grid(dotplot, histplot, 
                  nrow = 2, 
                  align = "v")

# save the result
ggsave(plot = plot,
       filename = "type_I_Err_inflation_comparison.pdf", 
       device = "pdf",
       width = TEXTWIDTH, 
       height = 0.5*TEXTHEIGHT)
