Rscript -e 'renv::activate(); renv::restore()'  # This step may take up to 15 minutes
Rscript -e 'devtools::package_info("symcrt")'
