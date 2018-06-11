
rm(list = ls())
setwd('.')

####################
# Install missing packages before everything
####################
source("scripts/util/__Util_MiscFunctions.R")
needed_packages <- c("reshape2", "igraph", "ggplot2", "msm", "dplyr", "tidyr", "gtools", "parallel", "snowfall")
install_missing_packages(needed_packages)

####################
# Source necessary scripts/libraries
####################
source("scripts/util/__Util__MASTER.R")
library(parallel)
library(snowfall)

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(5, 10, 20) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 10000 #number of generations to run simulation 
reps           <- 20 #number of replications per simulation (for ensemble)
chunk_size     <- 5 #number of simulations sent to single core 

# Threshold Parameters
ThreshM        <- rep(10, m) #population threshold means 
ThreshSD       <- ThreshM * 0.1 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.6, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
p              <- 0.5 #baseline probablity of initiating an interaction per time step
epsilon        <- 0.01 #relative weighting of social interactions for adjusting thresholds
beta           <- 1.1 #probability of interacting with individual in same state relative to others


####################
# Prep for Parallelization
####################
# Break up parameter replications into smaller batches

chunk_run  <- 1:(reps / chunk_size)
run_in_parallel <- expand.grid(n = Ns, run = chunk_run)
run_in_parallel <- run_in_parallel %>% 
  arrange(n)

# Prepare for parallel
no_cores <- detectCores() - 1
sfInit(parallel = TRUE, cpus = no_cores)
sfExportAll()
sfLibrary(dplyr)
sfLibrary(reshape2)
sfLibrary(igraph)
sfLibrary(ggplot2)
sfLibrary(msm)
sfLibrary(gtools)
sfLibrary(snowfall)
sfClusterSetupRNGstream(seed = 323)

####################
# Run ensemble simulation
####################
# Loop through group size (and chucnks)
improveSpec <- sfLapply(1:nrow(run_in_parallel), function(k) {
  # Set group size 
  n <- run_in_parallel[k, 1]
  chunk <- run_in_parallel[k, 2]
  # Prep lists for collection of simulation outputs from this group size
  ens_thresh      <- list()
  # Run Simulations
  for (sim in 1:chunk_size) {
    # Seed internal thresholds
    threshMat <- seed_thresholds(n = n, 
                                 m = m, 
                                 threshold_means = ThreshM,
                                 threshold_sds = ThreshSD)
    ens_thresh[[sim]] <- threshMat
  }
  return(ens_thresh)
})
sfStop()

test <- lapply(improveSpec, function(i) {
  bound <- do.call("rbind", i)
})

test <- do.call("rbind", test)
qplot(data = as.data.frame(test), x =Thresh1, y = Thresh2)
table(duplicated(test))