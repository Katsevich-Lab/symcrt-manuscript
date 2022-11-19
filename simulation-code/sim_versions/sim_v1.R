### define global parameters ###
alpha <- 0.05                # nominal level
maxType_I_null <- 0.99       # maximum Type-I error for naive GCM in null simulation
Type_I_alt <- 0.5            # Type-I error for naive GCM in alternative simulation
maxPower <- 0.99             # maximum power for oracle GCM in alternative simulation
test_type <- "two_side"      # sidedness of the tests
no_nu <- 5                   # number of values for nu parameter (effect of Z on Y)
no_theta <- 5                # number of values for theta parameter (effect of X on Y)
seed <- 4                    # seed to control random number generation
B_nu <- 20000                # number of replicates for nu computation
B_theta <- 20000             # number of replicates for theta computation
B_reps <- 400                # number of repetitions for each parameter setting
varying_values <- list(
  n = 100 * 2^seq(0, 4, 1),  # sample size
  d = 100 * 2^seq(0, 4, 1),  # dimension of Z
  s = 5 * 2^seq(0, 4, 1),    # number of variables in Z impacting Y
  rho = seq(0, 0.8, 0.2)     # auto-correlation parameter for Z
)
baseline_values <- list(n = 200,
                        d = 400,
                        s = 5,
                        rho = 0.4)
distributions <- c("binomial", "gaussian")
ways_to_learn <- c("supervised", "semi_supervised")
settings <- c("null", "calibration", "alternative")

simulations <- tidyr::crossing(
  distribution = distributions,
  way_to_learn = ways_to_learn,
  setting = settings
) |>
  as.data.frame()

method_strings <- list()
method_strings[["null"]] <- c(
  "GCM naive naive",
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
  "MaxwayCRT LASSO LASSO"
)
method_strings[["alternative"]] <- c(
  "GCM oracle oracle",
  "GCM LASSO LASSO",
  "GCM PLASSO PLASSO",
  "dCRT LASSO LASSO unnormalized",
  "dCRT PLASSO PLASSO unnormalized",
  "MaxwayCRT LASSO LASSO"
)
method_strings[["calibration"]] <- method_strings[["alternative"]]
