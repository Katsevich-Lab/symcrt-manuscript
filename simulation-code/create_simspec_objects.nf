// First, extract (distribution, way_to_learn, setting) for each simulation
process extract_sim_indices {
  memory '2GB'
  time '15m'

  input:
  val sim_version_fp from params.sim_version_fp

  output:
  path "simulations.txt" into simulations_raw_ch

  """
  #!/usr/bin/env Rscript

  source("$sim_version_fp")
  sapply(1:nrow(simulations), function(row)(as.list(simulations[row,]))) |> 
    unlist() |>
    writeLines("simulations.txt")
  """
}

simulations_ch = simulations_raw_ch.splitText().map{it.trim()}.collate(3)

// Second, create simulatr specifier objects for each simulation
process create_simspec_object {
  memory '2GB'
  time '1h'

  publishDir params.sim_spec_dir, mode: "copy"

  input:
  tuple val(distribution), val(way_to_learn), val(setting) from simulations_ch
  val sim_version_fp from params.sim_version_fp

  output:
  path 'sim_spec_*.rds' into output_ch

  """
  #!/usr/bin/env Rscript

  source('$sim_version_fp')
  sim_spec_obj <- symcrt::create_simspec_object(
    alpha, 
    maxType_I_null, 
    Type_I_alt, 
    maxPower, 
    test_type, 
    no_nu, 
    no_theta, 
    seed, 
    B_nu, 
    B_theta, 
    B_reps, 
    baseline_values,
    varying_values, 
    '$distribution',
    '$way_to_learn', 
    '$setting', 
    method_strings)
  simspec_filename <- sprintf(
    "sim_spec_%s_%s_%s.rds",
    '$distribution',
    '$way_to_learn',
    '$setting')
  saveRDS(sim_spec_obj, simspec_filename)
  """
}