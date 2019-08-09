################################################################################
#
# Install Missing Packages on Node (if running in parallel on cluster)
#
################################################################################
source("scripts/util/__Util_MiscFunctions.R")
needed_packages <- c("reshape2", "igraph", "ggplot2", "msm", "dplyr", "tidyr", "gtools", "parallel", "snowfall", "rlecuyer")
install_missing_packages(needed_packages)