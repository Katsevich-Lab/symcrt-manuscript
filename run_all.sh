#############################################################
# REPRODUCE ALL SIMULATIONS
#
# Note: This code took about 24 hours to run on a computer 
# cluster with 500 cores, with 8GB/core. The code will rerun
# computations to regenerate only those files that are not
# already present in this directory. 
#
#############################################################

# READ COMMAND LINE ARGUMENTS
# max_gb=$1    # GB available per core (default 8)
# profile=$2   # Nextflow profile (default "standard")
max_gb=7.5
profile="standard"

### INSTALL ALL R PACKAGES NECESSARY VIA RENV ###
Rscript -e 'if(!require(renv)) install.packages("renv", repos = "http://cran.us.r-project.org")'
Rscript -e 'renv::activate(); renv::restore()'  # This step may take up to 15 minutes
OLD_R_LIBS_USER=$R_LIBS_USER
export R_LIBS_USER=$(Rscript -e 'cat(.libPaths()[1])')

### REPRODUCE THE BENCHMARKING SIMULATION ###
sim_version_fp=$(pwd)"/simulation-code/sim_versions/sim_benchmarking.R"
sim_spec_dir=$(pwd)"/simulation-code/sim_spec_objects/benchmarking"
output_dir=$(pwd)"/simulation-results/benchmarking"
bash simulation-code/run_simulation.sh $sim_version_fp $sim_spec_dir $output_dir $max_gb $profile

### REPRODUCE THE DIAGNOSTIC SIMULATION ###
sim_version_fp=$(pwd)"/simulation-code/sim_versions/sim_diagnostic.R"
sim_spec_dir=$(pwd)"/simulation-code/sim_spec_objects/diagnostic"
output_dir=$(pwd)"/simulation-results/diagnostic"
bash simulation-code/run_simulation.sh $sim_version_fp $sim_spec_dir $output_dir $max_gb $profile

export R_LIBS_USER=$OLD_R_LIBS_USER

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
