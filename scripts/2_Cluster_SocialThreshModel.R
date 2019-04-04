################################################################################
#
# Social interaction model: set for running on Della cluster
#
################################################################################

rm(list = ls())

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
Ns             <- seq(5, 100, 5) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 200000 #number of generations to run simulation 
reps           <- 100 #number of replications per simulation (for ensemble)
chunk_size     <- 5 #number of simulations sent to single core 

# Threshold Parameters
ThreshM        <- rep(50, m) #population threshold means 
ThreshSD       <- ThreshM * 0 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.8, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
p              <- 1 #baseline probablity of initiating an interaction per time step
epsilon        <- 0.4 #relative weighting of social interactions for adjusting thresholds
beta           <- 1.1 #probability of interacting with individual in same state relative to others


####################
# Prep for Parallelization
####################
# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
dir_name <- paste0("Sigma", (ThreshSD/ThreshM)[1], "-Epsilon", epsilon, "-Beta", beta, "-LongRun")
full_path <- paste0(storage_path, dir_name)
dir.create(full_path)
sub_dirs <- c("TaskDist", "Entropy", "TaskTally", "Stim", 
              "Thresh", "Thresh1Time", "Thresh2Time", "Graphs")
for (sub_dir in sub_dirs) {
  dir.create(paste0(full_path, "/", sub_dir), showWarnings = FALSE)
}

# Break up parameter replications into smaller batches\
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
sfLibrary(stringr)
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
  ens_entropy     <- list()
  ens_thresh      <- list()
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
                                                   threshold_max = 100)
      # Update total task performance profile
      X_tot <- X_tot + X_g
    }
    
    ####################
    # Post run calculations
    ####################
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    entropy <- label_parallel_runs(matrix = entropy, n = n, simulation = sim, chunk = chunk)
    # Calculate total task distribution
    totalTaskDist <- X_tot / gens
    totalTaskDist <- label_parallel_runs(matrix = totalTaskDist, n = n, simulation = sim, chunk = chunk)
    # Create thresh table
    threshMat <- label_parallel_runs(matrix = threshMat, n = n, simulation = sim, chunk = chunk)
    # Add total task distributions, entropy values, and graphs to lists
    ens_taskDist[[sim]]    <- totalTaskDist
    ens_entropy[[sim]]     <- entropy
    ens_thresh[[sim]]      <- threshMat 
    ens_graphs[[sim]]      <- g_tot / gens
  }
  # Bind and write
  save_parallel_data(data = ens_taskDist, 
                     path = full_path, 
                     sub_directory = "TaskDist",
                     n = n, 
                     chunk = chunk)
  save_parallel_data(data = ens_entropy, 
                     path = full_path, 
                     sub_directory = "Entropy",
                     n = n, 
                     chunk = chunk)
  save_parallel_data(data = ens_thresh, 
                     path = full_path, 
                     sub_directory = "Thresh",
                     n = n, 
                     chunk = chunk)
  save(ens_graphs, 
       file = paste0(full_path,
                     "/Graphs/", 
                     str_pad(string = n, width = 3, pad = "0"), 
                     "-", 
                     str_pad(string = chunk, width = 2, pad = "0"), 
                     ".Rdata"))
  # Return all_clear
  rm(ens_taskDist, ens_entropy, ens_thresh, ens_graphs)
  return(paste0("DONE: n = ", n, 
                ", Sims ", ((chunk-1) * chunk_size)+1, "-", chunk * chunk_size))
  sys.sleep(1)
})

sfStop()

