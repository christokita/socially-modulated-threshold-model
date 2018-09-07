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
Ns             <- c(5, 10, seq(20, 100, 20)) #vector of number of individuals to simulate
Ns             <- c(100) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 50000 #number of generations to run simulation
reps           <- 100 #number of replications per simulation (for ensemble)

# Threshold Parameters
ThreshM        <- rep(50, m) #population threshold means 
InitialStim    <- rep(0, m) #intital vector of stimuli
Sigmas         <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15)
delta_values   <- seq(0.5, 0.8, 0.1) #vector of stimuli increase rates
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
# p              <- 1 #baseline probablity of initiating an interaction per time step
# epsilon        <- 0 #relative weighting of social interactions for adjusting thresholds
# beta           <- 1.1 #probability of interacting with individual in same state relative to others


####################
# Prep for Parallelization
####################
# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
dir_name <-"FixedThreshold-SigmaDeltaSweep"
full_path <- paste0(storage_path, dir_name, "/")
dir.create(full_path, showWarnings = FALSE)

# Prepare table for values to be iterated over
parameter_values <- expand.grid(n = Ns, D = delta_values, S = Sigmas)

# Check if file already exists
files <- list.files(full_path)
completed_runs <- data.frame(n = as.numeric(gsub(x = files, "^.*-n([0-9]+)\\.Rdata$", "\\1", perl = T)))
completed_runs$D <- as.numeric(gsub(x = files, "^Delta([\\.0-9]+)-.*$", "\\1", perl = T))
completed_runs$S <- as.numeric(gsub(x = files, "^.*-Sigma([\\.0-9]+)-.*$", "\\1", perl = T))
parameter_values <- anti_join(parameter_values, completed_runs, by = c("n", "D", "S"))
parameter_values <- parameter_values %>% 
  select(D, S) %>% 
  unique()

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
# LOOP 1: loop through threshold variation and deltas
parallel_simulations <- sfLapply(1:nrow(parameter_values), function(k) {
  # Set threshold variation
  ThreshSD <- ThreshM * parameter_values[k, "S"]
  # Set deltas
  deltas <- rep(parameter_values[k, "D"], m)
  # LOOP 2: loop through group sizes
  for (i in 1:length(Ns)) {
    # Set group size
    n <- Ns[i]
    # Prep lists for collection of simulation outputs from this group size
    ens_taskDist    <- list()
    ens_entropy     <- list()
    ens_thresh      <- list()
    ens_graphs      <- list()
    # LOOP 3: Loop through replicates
    for (sim in 1:reps) {
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
      # g_tot <-  matrix(data = rep(0, n * n), ncol = n)
      # colnames(g_tot) <- paste0("v-", 1:n)
      # rownames(g_tot) <- paste0("v-", 1:n)
      
      ####################
      # Simulate individual run
      ####################
      # LOOP 4: Run simulation
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
        # g_adj <- temporalNetwork(X_sub_g = X_g,
        #                          prob_interact = p,
        #                          bias = beta)
        # g_tot <- g_tot + g_adj
        # # Adjust thresholds
        # threshMat <- adjust_thresholds_social_capped(social_network = g_adj,
        #                                              threshold_matrix = threshMat,
        #                                              state_matrix = X_g,
        #                                              epsilon = epsilon,
        #                                              threshold_max = 2 * ThreshM[1])
        # Update total task performance profile
        X_tot <- X_tot + X_g
      }
      
      ####################
      # Post run calculations
      ####################
      # Calculate Entropy
      entropy <- as.data.frame(mutualEntropy(TotalStateMat = X_tot))
      entropy$n <- n
      entropy$replicate <- sim
      entropy$sigma <- ThreshSD[1] / ThreshM[1]
      entropy$delta <- deltas[1]
      # Calculate total task distribution
      total_task_dist <- as.data.frame(X_tot / gens)
      total_task_dist$n <- n
      total_task_dist$replicate <- sim
      total_task_dist$sigma <- ThreshSD[1] / ThreshM[1]
      total_task_dist$delta <- deltas[1]
      # Add total task distributions, entropy values, and graphs to lists
      ens_taskDist[[sim]]    <- total_task_dist
      ens_entropy[[sim]]     <- entropy
    }
    # Bind and save
    task_distributions <- do.call("rbind", ens_taskDist)
    entropies <- do.call("rbind", ens_entropy)
    save(task_distributions, entropies, file = paste0(full_path, "Delta", deltas[1],
                                                      "-Sigma", ThreshSD[1]/ThreshM[1],
                                                      "-n", n, ".Rdata"))
    # Return all_clear
    rm(ens_taskDist, ens_entropy)
    Sys.sleep(1)
  }
})

sfStop()

