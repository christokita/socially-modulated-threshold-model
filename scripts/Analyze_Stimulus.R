################################################################################
#
# Analyze stimulus data
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")

####################
# Load data and summarise
####################
files <- list.files("output/Rdata/_ProcessedData/Stim/Sigma0-Epsilon0.1-Beta1.1/", full.names = TRUE)
stim_data <- lapply(files, function(file) {
  # Load group size data
  load(file)
  # 
})


