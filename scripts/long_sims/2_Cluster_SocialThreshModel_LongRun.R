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
n              <- 80 #vector of number of individuals to simulate
m              <- 2 #number of tasks
Tsteps         <- 500000 #number of time steps to run simulation 
reps           <- 100 #number of replications per simulation (for ensemble)

# Threshold Parameters
ThreshM        <- rep(50, m) #population threshold means 
ThreshSD       <- ThreshM * 0 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.8, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
p              <- 1 #baseline probablity of initiating an interaction per time step
epsilon        <- 0.1 #relative weighting of social interactions for adjusting thresholds
beta           <- 1.075 #probability of interacting with individual in same state relative to others


####################
# Prep for Parallelization
####################
# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
dir_name <- paste0("LongSim-Sigma", (ThreshSD/ThreshM)[1], "-Epsilon", epsilon, "-Beta", beta)
full_path <- paste0(storage_path, dir_name)
dir.create(full_path)
sub_dirs <- c("TaskDist", "Entropy", "Stim", "Thresh", "Graphs")
for (sub_dir in sub_dirs) {
  dir.create(paste0(full_path, "/", sub_dir), showWarnings = FALSE)
}


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
parallel_simulations <- sfLapply(1:reps, function(sim) {
  
    ####################
    # Seed structures and intial matrices
    ####################
    # Set initial probability matrix (P_g)
    P_g <- matrix(data = rep(0, n * m), ncol = m)
    # Seed task (external) stimuli
    stimMat <- seed_stimuls(intitial_stim = InitialStim, 
                            Tsteps = Tsteps)
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
    for (t in 1:Tsteps) { 
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
    # Calculate total task distribution
    totalTaskDist <- X_tot / Tsteps
    # Normalize interaciton matrix
    g_tot <- g_tot / Tsteps
    
    ####################
    # Save data
    ####################
    save(entropy, file = paste0(full_path, "/Entropy/Sim_", sim, ".Rdata"))
    save(totalTaskDist, file = paste0(full_path, "/TaskDist/Sim_", sim, ".Rdata"))
    save(stimMat, file = paste0(full_path, "/Stim/Sim_", sim, ".Rdata"))
    save(threshMat, file = paste0(full_path, "/Thresh/Sim_", sim, ".Rdata"))
    save(g_tot, file = paste0(full_path, "/Graphs/Sim_", sim, ".Rdata"))
    
})

sfStop()

