################################################################################
#
# Social interaction model: sweep across beta values with set order of simtulus encounter
# designed for running in parallel on cluster
#
# **This script only focuses on measuring resulint DOL and task distributions**
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
# Choose task probabilistically 
####################
update_task_performance <- function(task_probs, state_matrix, quit_prob) {
  # Create possible task space
  tasks <- seq(1:ncol(task_probs))
  # Loop through individuals
  for(row in 1:nrow(state_matrix)) {
    # Inactive workers randomly sample one stimulus
    if (sum(state_matrix[row, ]) == 0) {
      # Sample task probability
      #tasks_order <- sample(x = tasks, size = length(tasks), replace = FALSE)
      tasks_order <- tasks
      # Loop through tasks and go with first one that results in activity
      for (task in tasks_order) {
        prob <- task_probs[row, task]
        activity <- sample(x = c(0, 1), size = 1, prob = c(1 - prob, prob))
        if (activity == 1) {
          state_matrix[row, task] <- activity
          break
        }
      }
    } 
    else { #active workers quit with certain probability
      quit_now <- sample(x = c("yes", "no"), size = 1, prob = c(quit_prob, (1 - quit_prob)))
      if (quit_now == "yes") {
        state_matrix[row, ] <- 0
      }
    }
  }
  # Return
  colnames(state_matrix) <- paste0("Task", 1:ncol(task_probs))
  rownames(state_matrix) <- paste0("v-", 1:nrow(task_probs))
  return(state_matrix) 
}

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
n              <- 80 #group size
m              <- 2 #number of tasks
Tsteps         <- 50000 #number of time steps to run simulation 
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
epsilon        <- 0.1 #relative weighting of social interactions for adjusting thresholds
betas          <- seq(1, 1.25, 0.01) #probability of interacting with individual in same state relative to others


####################
# Prep for Parallelization
####################
# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
dir_name <- paste0("n", n,  "-Sigma", (ThreshSD/ThreshM)[1], "-Epsilon", epsilon, "_BetaSweep_stimorder")
full_path <- paste0(storage_path, dir_name)
dir.create(full_path)
sub_dirs <- c("TaskDist", "Entropy")
for (sub_dir in sub_dirs) {
  dir.create(paste0(full_path, "/", sub_dir), showWarnings = FALSE)
}

# Break up parameter replications into smaller batches\
chunk_run  <- 1:(reps / chunk_size)
run_in_parallel <- expand.grid(beta =  round(betas, digits = 3), run = chunk_run) #rounding to make sure numbers are what they appear
run_in_parallel <- run_in_parallel %>% 
  arrange(beta)

# Read files already simulated and filter out of needed parameters
already_ran <- list.files(path = paste0(full_path, "/Entropy/"))
ran_files <- lapply(already_ran, function(f) {
  beta_val <- as.numeric(gsub("([\\.0-9]+)-[0-9]+\\.Rdata", "\\1", x = f, perl = TRUE))
  chunk_num <- as.numeric(gsub("[\\.0-9]+-([0-9]+)\\.Rdata", "\\1", x = f, perl = TRUE))
  to_return <- data.frame(beta = beta_val, run = chunk_num)
  return(to_return)
})
ran_files <- do.call('rbind', ran_files)
run_in_parallel <- run_in_parallel %>% 
  anti_join(ran_files)

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
  beta <- run_in_parallel[k, 1]
  chunk <- run_in_parallel[k, 2]
  # Prep lists for collection of simulation outputs from this group size
  ens_taskDist    <- list()
  ens_entropy     <- list()
  # Run Simulations
  for (sim in 1:chunk_size) {
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
                                                   threshold_max = Inf)
      # Update total task performance profile
      X_tot <- X_tot + X_g
    }
    
    ####################
    # Post run calculations
    ####################
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    entropy <- label_parallel_runs_parmeters(matrix = entropy, beta = beta, epsilon = epsilon, simulation = sim, chunk = chunk)
    # Calculate total task distribution
    totalTaskDist <- X_tot / Tsteps
    totalTaskDist <- label_parallel_runs_parmeters(matrix = totalTaskDist, beta = beta, epsilon = epsilon, simulation = sim, chunk = chunk)
    # Add total task distributions, entropy values, and graphs to lists
    ens_taskDist[[sim]]    <- totalTaskDist
    ens_entropy[[sim]]     <- entropy
  }
  # Bind and write
  save_parallel_data_parameter(data = ens_taskDist, 
                               path = full_path, 
                               sub_directory = "TaskDist",
                               parameter_value = beta, 
                               chunk = chunk)
  save_parallel_data_parameter(data = ens_entropy, 
                               path = full_path, 
                               sub_directory = "Entropy",
                               parameter_value = beta, 
                               chunk = chunk)
  # Return all_clear
  rm(ens_taskDist, ens_entropy)
  return(paste0("DONE: beta = ", beta, 
                ", Sims ", ((chunk-1) * chunk_size)+1, "-", chunk * chunk_size))
  sys.sleep(1)
})

sfStop()

