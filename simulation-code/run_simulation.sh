# 0. Parse command line arguments
sim_version_fp=$1  # path to file specifying the parameters of the simulation
sim_spec_dir=$2    # path to directory to store simulatr specifier objects
output_dir=$3      # path to directory to store output files
max_gb=$4          # GB available per core (default 8)
profile=$5         # Nextflow profile (default "standard")
# set defaults
if [ -z "$profile" ]; then
  profile="standard"
fi
if [ -z "$max_gb" ]; then
  max_gb=8
fi

echo "Launching the simulation..."

# 1. Create output directory if it does not exist
mkdir -p $output_dir

# 2. Make sure pipelines/packages are up to date
Rscript -e 'devtools::install_github("Katsevich-Lab/simulatr")'
Rscript -e 'devtools::install_git("git@github.com:Katsevich-Lab/symcrt.git")'
nextflow pull katsevich-lab/simulatr-pipeline

# 3. Create simulatr specifier objects
echo "Creating simulatr specifier objects"
if [ ! -d "$sim_spec_dir" ]; then
  mkdir $sim_spec_dir
  nextflow run create_simspec_objects.nf \
      --sim_version_fp $sim_version_fp \
      --sim_spec_dir $sim_spec_dir \
      -profile $profile
else
  echo "Simulatr specifier objects already created"
fi

# 4. Run all of the simulations
for simspec_filename in $sim_spec_dir/*
do
  # obtain simulation string from simspec_filename
  sim_string=$(Rscript -e '
  simspec_filename <- commandArgs(trailingOnly=TRUE)[1];
  sim_string = strsplit(simspec_filename, split = "sim_spec_|.rds")[[1]][2]
  cat(sim_string)' $simspec_filename)
  output_filename=$sim_string"_results.rds"
  if [ ! -f "$output_dir/$output_filename" ]; then
    echo "Running the "$sim_string" simulation..."
    nextflow run katsevich-lab/simulatr-pipeline \
      --simulatr_specifier_fp $simspec_filename \
      --result_dir $output_dir \
      --result_file_name $output_filename \
      --max_gb $max_gb \
      -resume \
      -profile $profile \
      -with-trace $output_dir/$sim_string"_trace.txt"
    wait
    if [ -f "$output_dir/$output_filename" ]; then
      if [ -z "$NXF_WORK" ]; then
        rm -r $NXF_WORK/*
      fi
    else
        exit
    fi
  else
    echo $output_filename" already exists"
  fi
done

echo "Done."
