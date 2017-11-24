################################################################################
#
# Bind together extra group sizes from same parameter combo
#
################################################################################

load("output/Rdata/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1_LargerGroups.Rdata")

groups_taskDist1       <- groups_taskDist
groups_taskCorr1       <- groups_taskCorr
groups_taskStep1       <- groups_taskStep
groups_taskTally1      <- groups_taskTally
groups_stim1           <- groups_stim
groups_entropy1        <- groups_entropy
groups_graphs1         <- groups_graphs
groups_specialization1 <- groups_specialization
groups_thresh1         <- groups_thresh

load("output/Rdata/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1_SmallerGroups.Rdata")

groups_taskDist       <- c(groups_taskDist, groups_taskDist1)
groups_taskCorr       <- c(groups_taskCorr, groups_taskCorr1)
groups_taskStep       <- c(groups_taskStep, groups_taskStep1)
groups_taskTally      <- c(groups_taskTally, groups_taskTally1)
groups_stim           <- c(groups_stim, groups_stim1)
groups_entropy        <- c(groups_entropy, groups_entropy1)
groups_graphs         <- c(groups_graphs, groups_graphs1)
groups_specialization <- c(groups_specialization, groups_specialization1)
groups_thresh         <- c(groups_thresh, groups_thresh1)

save(groups_taskDist, groups_taskCorr, groups_taskStep, groups_taskTally, groups_stim, groups_entropy, groups_graphs,
     groups_specialization, groups_thresh,
     file = "output/Rdata/Sigma001-Eps001-Phi001-ConnectP01-Bias1.1.Rdata")
