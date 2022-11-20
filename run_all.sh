#############################################################
# REPRODUCE THE SIMULATIONS
#############################################################

max_gb=7.5
profile="standard"

### INSTALL ALL R PACKAGES NECESSARY ###
Rscript -e 'install.packages("packrat")'
Rscript -e 'packrat::restore()'

### REPRODUCE THE V1 SIMULATION ###
sim_version_fp=$(pwd)"/simulation-code/sim_versions/sim_v1.R"
sim_spec_dir=$(pwd)"/simulation-code/sim_spec_objects/v1"
output_dir=$(pwd)"/simulation-results/v1"
bash simulation-code/run_simulation.sh $sim_version_fp $sim_spec_dir $output_dir $max_gb $profile

### REPRODUCE THE DIAGNOSTIC SIMULATION ###
sim_version_fp=$(pwd)"/simulation-code/sim_versions/sim_diagnostic.R"
sim_spec_dir=$(pwd)"/simulation-code/sim_spec_objects/diagnostic"
output_dir=$(pwd)"/simulation-results/diagnostic"
bash simulation-code/run_simulation.sh $sim_version_fp $sim_spec_dir $output_dir $max_gb $profile

### REPRODUCE THE CONFOUNDING SIMULATION ###
Rscript simulation-code/assess-confounding-level.R

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