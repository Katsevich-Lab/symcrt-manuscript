# Reconciling model-X and doubly robust approaches to conditional independence testing

This repository reproduces the results reported in arXiv version 1 the following paper:

Z. Niu, A. Chakraborty, O. Dukes, E. Katsevich.
"Reconciling model-X and doubly robust approaches to conditional independence testing." 
([arXiv]())

## Repository structure

The code to reproduce the numerical simulations is in `simulation-code`, as well as the accompanying R package [symcrt](https://github.com/katsevich-lab/symcrt). Some of the simulations are based on the `simulatr` [R package](https://github.com/timothy-barry/simulatr) and associated [Nextflow pipeline](https://github.com/katsevich-lab/simulatr-pipeline). The `simulatr` specifier objects are stored in `simulation-code/sim_spec_objects`. The results outputted by the simulation code are stored in `simulation-results`. These results are split into three parts: the primary simulations of the paper (underlying all figures but Figures 2 and 5), the simulation to assess the confounding levels of previous studies (underlying Figure 2), and the simulation to compare the estimation performances of the lasso and post lasso (underlying Figure 5). The results for these parts are stored in subdirectories `benchmarking`, `confounding`, and `diagnostic` of `simulation-results`, respectively. Code to produce figures based on these results is stored in `plotting-code`, and the corresponding figures are stored in `figures`. The manuscript LaTeX source is stored in `manuscript`.

## Instructions for reproducing (a subset of) the results in the paper

- Identify a computing environment (e.g. a laptop, high performance computing cluster, or cloud computing platform) you would like to use. If you would like to reproduce a subset of Figures 3-13, a distributed computing environment (i.e. cluster or cloud) is strongly recommended.
- Check that you have R version 4.1.3 or higher installed. Copies of all necessary R packages will be installed for you via [`renv`](https://rstudio.github.io/renv/articles/renv.html). 
- [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) or upgrade to the latest version, and [configure](https://www.nextflow.io/docs/latest/config.html#) it for use in your computing environment.
- Clone this repository into your computing environment.
- Optional: Determine how much RAM is available to each core (code defaults to 8GB).
- Optional/advanced: Determine which [Nextflow profile](https://www.nextflow.io/docs/latest/config.html#config-profiles) you would like to use (code defaults to `standard`). 
- Delete any pre-computed objects you would like to re-compute, including `simulatr` specifier objects under `simulation-code/sim_spec_objects`, simulation results under `simulation-results`, or figures under `figures`.
- Open a shell and navigate to the `symcrt-manuscript` directory. 
- Launch the computations via `bash run_all.sh <MAXGB> <PROFILE>`, where `<MAXGB>` should be replaced by the available RAM per core in GB (or omitted if you would like to use the default of 8GB) and `<PROFILE>` should be replaced by the Nextflow profile (or omitted if not using profiles). This will re-compute (only) the deleted pre-computed objects. If no pre-computed objects are deleted, then the code will not do anything.

## Acknowledgments

ZN was partially supported by the grant "Statistical Software for Single Cell CRISPR Screens" awarded to EK by Analytics at Wharton. OD was partially supported by FWO grant 1222522N and NIH grant AG065276. EK was partially supported by NSF DMS-2113072.  We acknowledge help from Timothy Barry with our simulation studies and the underlying computational infrastructure, including his `simulatr` R package and Nextflow pipeline. We acknowledge dedicated support from the staff at the Wharton High Performance Computing Cluster. We acknowledge Lucas Janson for providing details about the simulation setting in Candes et al (2018). We acknowledge Eric Tchetgen Tchetgen for helpful discussions on hypothesis testing in the semiparametric models.
