################################################################################
#
# Entropy across parameter
#
################################################################################

rm(list = ls())
source("scripts/util/__Util__MASTER.R")
library(RColorBrewer)
library(scales)


####################
# Load Data: epsilon
####################
directory <- "output/Rdata/GroupSizeEpsilonSweep_Sigma0-Beta1.1/"

# List files 
files <- list.files(directory, full.names = TRUE)
for (file in files) {
  load(file)
  epsilon <- as.numeric(gsub(".*epsilon([\\.0-9]+).Rdata", "\\1", file, perl = TRUE))
  entropy_sum$epsilon <- epsilon
  if (!exists("entropy")) {
    entropy <- entropy_sum
  } else {
    entropy <- rbind(entropy, entropy_sum)
  }
}

save(entropy, file = "output/ParameterSpace/GroupSizeEpsilonSweep_Sigma0-Beta1.1.Rdata")

rm(entropy)

####################
# Load Data: beta
####################
directory <- "output/Rdata/GroupSizeBetaSweep_Sigma0-Epsilon0.1/" 
# List files 
files <- list.files(directory, full.names = TRUE)
for (file in files) {
  load(file)
  if (!exists("entropy")) {
    entropy <- entropy_sum
  } else {
    entropy <- rbind(entropy, entropy_sum)
  }
}

save(entropy, file = "output/ParameterSpace/GroupSizeBetaSweep_Sigma0-Epsilon0.1.Rdata")

