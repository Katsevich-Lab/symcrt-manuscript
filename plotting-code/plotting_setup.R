# load libraries
suppressPackageStartupMessages({
  library(RColorBrewer)
  library(reshape)    # for merge_all()
  library(MASS)       # for fitdistr()
  library(gridExtra)
  library(grid)
  library(latex2exp)
  library(cowplot)    # for plot_grid()
  library(simulatr)
  library(symcrt)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

# dimensions of the manuscript page
TEXTWIDTH = 6.3
TEXTHEIGHT = 8.64

# theme for plotting
plotting_theme <- theme_bw()