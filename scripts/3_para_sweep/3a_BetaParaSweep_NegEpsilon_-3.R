################################################################################
#
# Social interaction model: Sweep beta and group size parameter space
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
Tsteps         <- 50000 #number of time steps to run simulation
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
epsilon        <- -0.1 #relative weighting of social interactions for adjusting thresholds
betas          <- seq(0.85, 0.89, 0.01) #probability of interacting with individual in same state relative to others


####################
# Prep for Parallelization
####################
# Create parameter combinations for parallelization
run_in_parallel <- expand.grid(n = Ns, beta = betas)
run_in_parallel <- run_in_parallel %>% 
  arrange(n)

# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
file_name <- paste0("GroupSizeBetaSweep_Sigma", ThreshSD[1], "-Epsilon", epsilon)
full_path <- paste0(storage_path, file_name, '/')
dir.create(full_path, showWarnings = FALSE)

# Check if there is already some runs done
files <- list.files(full_path)
completed_runs <- data.frame(n = as.numeric(gsub(x = files, "n([0-9]+)-.*", "\\1", perl = T)))
completed_runs$beta <- as.numeric(gsub(x = files, ".*-beta([\\.0-9]+).Rdata$", "\\1", perl = T))
run_in_parallel <- anti_join(run_in_parallel, completed_runs, by = c("n", "beta"))

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
# sfClusterSetupRNGstream(seed = 100)

####################
# Run ensemble simulation
####################
# Loop through group size (and chucnks)
parallel_simulations <- sfLapply(1:nrow(run_in_parallel), function(k) {
  # Set group size 
  n <- run_in_parallel[k, 1]
  beta <- run_in_parallel[k, 2]
  # Prep lists for collection of simulation outputs from this group size
  ens_entropy     <- list()
  # Run Simulations
  for (sim in 1:reps) {
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
                                                   threshold_max = 2 * ThreshM[1])
      # Update total task performance profile
      X_tot <- X_tot + X_g
    }
    ####################
    # Post run calculations
    ####################
    # Calculate Entropy
    entropy <- as.data.frame(mutualEntropy(TotalStateMat = X_tot))
    entropy$n <- n
    entropy$beta <- beta
    # Add entropy values to list
    ens_entropy[[sim]] <- entropy
    # Clean
    rm(X_tot, stimMat, threshMat, g_tot, g_adj, P_g, X_g)
  }
  # Bind together and summarise
  entropy_sum <- do.call("rbind", ens_entropy)
  entropy_sum <- entropy_sum %>% 
    group_by(n, beta) %>% 
    summarise(Dsym_mean = mean(Dsym),
              Dysm_SD = sd(Dsym),
              Dtask_mean = mean(Dtask),
              Dtask_SD = sd(Dtask),
              Dind_mean = mean(Dind),
              Dind_SD = sd(Dind))
  entropy_sum <- as.data.frame(entropy_sum)
  save(entropy_sum, file = paste0(full_path, 
                                  "n", 
                                  str_pad(string = n, width =  3, pad =  "0"),
                                  "-beta",
                                  beta, 
                                  ".Rdata"))
  Sys.sleep(1)
})

sfStop()

# Bind and save
# parallel_data <- do.call('rbind', parallel_simulations)

# Create directory for depositing data
# storage_path <- "/scratch/gpfs/ctokita/"
# file_name <- paste0("GroupSizeBetaSweep_Sigma", ThreshSD[1], "-Epsilon", epsilon)
# full_path <- paste0(storage_path, file_name, '.Rdata')
# save(parallel_data, file = full_path)
