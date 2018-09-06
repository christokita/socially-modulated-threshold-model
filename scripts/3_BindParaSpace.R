################################################################################
#
# Entropy across parameter
#
################################################################################

####################
# Load Data: epsilon
####################
rm(list = ls())
directory <- "output/Rdata/GroupSizeEpsilonSweep_Sigma0-Beta0.9/"

# List files 
files <- list.files(directory, full.names = TRUE)
for (file in files) {
  load(file)
  epsilon <- as.numeric(gsub(".*-epsilon([0-9-\\.]+).Rdata", "\\1", file, perl = TRUE))
  entropy_sum$epsilon <- epsilon
  if (!exists("entropy")) {
    entropy <- entropy_sum
  } else {
    entropy <- rbind(entropy, entropy_sum)
  }
}

save(entropy, file = "output/ParameterSpace/GroupSizeEpsilonSweep_Sigma0-Beta0.9.Rdata")

rm(entropy)

####################
# Load Data: beta
####################
rm(list = ls())
directory <- "output/Rdata/GroupSizeBetaSweep_Sigma0-Epsilon-0.1/" 
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

save(entropy, file = "output/ParameterSpace/GroupSizeBetaSweep_Sigma0-Epsilon-0.1.Rdata")

