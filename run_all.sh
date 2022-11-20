#############################################################
# REPRODUCE THE SIMULATIONS
#############################################################

### REPRODUCE THE V1 SIMULATION ###
# TBD: Gene

### REPRODUCE THE DIAGNOSTIC SIMULATION ###
# TBD: Gene

### REPRODUCE THE CONFOUNDING SIMULATION ###
# TBD: Gene

#############################################################
# PRODUCE THE PLOTS
#############################################################

# Create Figure 1
Rscript plotting-code/figure_1_negative_results.R

# Create Figure 2
Rscript plotting-code/figure_2_confounding_comparison.R

# Create Figures 3-4 
Rscript plotting-code/figure_3_gaussian_supervised_setting_partial_type_I_error.R
Rscript plotting-code/figure_4_gaussian_supervised_setting_partial_power.R

# Create Figure 5
Rscript plotting-code/figure_5_MSE.R

# Create Figures 6-13
Rscript plotting-code/figures_6_8_10_12_full_type_I_error.R
Rscript plotting-code/figures_7_9_11_13_power.R