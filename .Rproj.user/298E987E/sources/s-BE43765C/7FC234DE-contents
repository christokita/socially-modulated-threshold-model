################################################################################
#
# Social interaction model: set for running on Della cluster
#
################################################################################

rm(list = ls())

####################
# Install missing packages before everything
####################
source("scripts/util/__Util_MiscFunctions.R")
needed_packages <- c("reshape2", "igraph", "ggplot2", "msm", "dplyr", "tidyr", "gtools", "parallel", "snowfall", "rlecuyer")
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
reps           <- 10 #number of replications per simulation (for ensemble)
chunk_size     <- 5 #number of simulations sent to single core 

# Threshold Parameters
ThreshM        <- rep(10, m) #population threshold means 
ThreshSD       <- ThreshM * 0.0 #population threshold standard deviations
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
no_cores <- detectCores() 
sfInit(parallel = TRUE, cpus = no_cores)
sfExportAll()
sfLibrary(dplyr)
sfLibrary(reshape2)
sfLibrary(igraph)
sfLibrary(ggplot2)
sfLibrary(msm)
sfLibrary(gtools)
sfLibrary(snowfall)
sfLibrary(tidyr)
sfClusterSetupRNGstream(seed = 323)

####################
# Run ensemble simulation
####################
# Loop through group size (and chucnks)
parallel_simulations <- sfLapply(1:nrow(run_in_parallel), function(k) {
  # Set group size 
  n <- run_in_parallel[k, 1]
  chunk <- run_in_parallel[k, 2]
  # Prep lists for collection of simulation outputs from this group size
  ens_taskDist    <- list()
  ens_taskTally   <- list()
  ens_entropy     <- list()
  ens_stim        <- list()
  ens_thresh      <- list()
  ens_thresh1Time <- list()
  ens_thresh2Time <- list()
  ens_graphs      <- list()
  # Run Simulations
  for (sim in 1:chunk_size) {
    ####################
    # Seed structures and intial matrices
    ####################
    # Set initial probability matrix (P_g)
    P_g <- matrix(data = rep(0, n * m), ncol = m)
    # Seed task (external) stimuli
    stimMat <- seed_stimuls(intitial_stim = InitialStim, 
                            gens = gens)
    # Seed internal thresholds
    threshMat <- seed_thresholds(n = n, 
                                 m = m, 
                                 threshold_means = ThreshM,
                                 threshold_sds = ThreshSD)
    # Start task performance
    X_g <- matrix(data = rep(0, length(P_g)), ncol = ncol(P_g))
    # Create cumulative task performance matrix
    X_tot <- X_g
    # Create cumulative adjacency matrix
    g_tot <-  matrix(data = rep(0, n * n), ncol = n)
    colnames(g_tot) <- paste0("v-", 1:n)
    rownames(g_tot) <- paste0("v-", 1:n)
    # Prep lists for data collection within simulation
    taskTally <- list()
    thresh1time <- list()
    thresh2time <- list()
    thresh1time[[1]] <- threshMat[ ,1]
    thresh2time[[1]] <- threshMat[ ,2]
    
    ####################
    # Simulate individual run
    ####################
    # Run simulation
    for (t in 1:gens) { 
      # Current timestep is actually t+1 in this formulation, because first row is timestep 0
      # Update stimuli
      stimMat <- update_stim(stim_matrix = stimMat, 
                             deltas = deltas, 
                             alpha = alpha, 
                             state_matrix = X_g, 
                             time_step = t)
      # Calculate task demand based on global stimuli
      P_g <- calc_determ_thresh(time_step        = t + 1, # first row is generation 0
                                threshold_matrix = threshMat, 
                                stimulus_matrix  = stimMat)
      # Update task performance
      X_g <- update_task_performance(task_probs   = P_g,
                                     state_matrix = X_g,
                                     quit_prob    = quitP)
      # Update social network (previously this was before probability/task update)
      g_adj <- temporalNetwork(X_sub_g = X_g,
                               prob_interact = p,
                               bias = beta)
      g_tot <- g_tot + g_adj
      # Adjust thresholds
      threshMat <- adjust_thresholds_social_capped(social_network = g_adj,
                                                   threshold_matrix = threshMat,
                                                   state_matrix = X_g,
                                                   epsilon = epsilon,
                                                   threshold_max = 2 * ThreshM[1])
      # Capture threshold values
      thresh1time[[t + 1]] <- threshMat[,1]
      thresh2time[[t + 1]] <- threshMat[,2]
      # Update total task performance profile
      X_tot <- X_tot + X_g
      # Capture current task performance tally
      tally <- matrix(c(t, colSums(X_g)), ncol = ncol(X_g) + 1)
      colnames(tally) <- c("t", colnames(X_g))
      taskTally[[t]] <- tally
    }
    
    ####################
    # Post run calculations
    ####################
    # Bind together task tally
    col_names <- colnames(taskTally[[1]])
    taskTally <- matrix(unlist(taskTally), 
                        ncol = length(taskTally[[1]]), 
                        byrow = TRUE, 
                        dimnames = list(c(NULL), c(col_names)))
    test <- label_parallel_runs(matrix = taskTally, n = n, simulation = sim, chunk = chunk)
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    entropy <- label_parallel_runs(matrix = entropy, n = n, simulation = sim, chunk = chunk)
    # Calculate total task distribution
    totalTaskDist <- X_tot / gens
    totalTaskDist <- label_parallel_runs(matrix = totalTaskDist, n = n, simulation = sim, chunk = chunk)
    # Create tasktally table
    stimMat <- label_parallel_runs(matrix = stimMat, n = n, simulation = sim, chunk = chunk)
    # Thresh tracking matrices
    thresh1time <- summarise_threshold_tracking(tracked_threshold = thresh1time, n = n, time_steps = gens)
    thresh2time <- summarise_threshold_tracking(tracked_threshold = thresh2time, n = n, time_steps = gens)
    # Add total task distributions, entropy values, and graphs to lists
    ens_taskDist[[sim]]    <- totalTaskDist
    ens_entropy[[sim]]     <- entropy
    ens_taskTally[[sim]]   <- taskTally
    ens_stim[[sim]]        <- stimMat
    ens_thresh[[sim]]      <- threshMat 
    ens_thresh1Time[[sim]] <- thresh1time
    ens_thresh2Time[[sim]] <- thresh2time
    ens_graphs[[sim]]      <- g_tot / gens
  }
  # Join together
  all_data <- list()
  all_data[["task_distribution"]] <- ens_taskDist
  all_data[["entropy"]] <- ens_entropy
  all_data[["task_tally"]] <- ens_taskTally
  all_data[["stimulus"]] <- ens_stim
  all_data[["threshold_matrix"]] <- ens_thresh
  all_data[["threshold1_time"]] <- ens_thresh1Time
  all_data[["threshold2_time"]] <- ens_thresh2Time
  all_data[["socail_network"]] <- ens_graphs
  # Return
  return(all_data)
})

sfStop()

####################
# Save data
####################
filename <- paste0("Sigma", ThreshSD[1], "-Epsilon", epsilon, "-Beta", beta, ".Rdata")
save(parallel_simulations, file = filename)
