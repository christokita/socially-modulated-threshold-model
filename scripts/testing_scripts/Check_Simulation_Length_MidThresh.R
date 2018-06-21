################################################################################
#
# Social interaction model: test simulation length time for affect on results
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
Ns             <- c(5, 10, 20, 30, 40, 50, 100) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 70000 #number of generations to run simulation 
reps           <- 100 #number of replications per simulation (for ensemble)
chunk_size     <- 10 #number of simulations sent to single core 

# Threshold Parameters
ThreshM        <- rep(25, m) #population threshold means 
ThreshSD       <- ThreshM * 0 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.8, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
p              <- 0.5 #baseline probablity of initiating an interaction per time step
epsilon        <- 0.05 #relative weighting of social interactions for adjusting thresholds
beta           <- 1.1 #probability of interacting with individual in same state relative to others

# Where to check time values
avg_window <- 100
times <- c(1000, seq(5000, 70000, 5000))

####################
# Prep for Parallelization
####################
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
  check_info <- list()
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
                                                   threshold_max = 2 * ThreshM[1])
      # Update total task performance profile
      X_tot <- X_tot + X_g
      # Capture stats if it is appropriate window
      if (t %in% times) {
        # Get DOL
        window_entropy <- mutualEntropy(X_tot)
        # Get threshold information
        thresh_info <- matrix(data = c(sd(threshMat[, 1]),
                             max(threshMat[, 1]),
                             min(threshMat[, 1]),
                             sd(threshMat[, 2]),
                             max(threshMat[, 2]),
                             min(threshMat[, 2])), 
                             nrow = 1, 
                             dimnames = list(NULL, 
                                             c("Thresh1SD", "Thresh1Max", "Thresh1Min", 
                                               "Thresh2SD", "Thresh2Max", "Thresh2Min")))
        # Get stim information
        stim_info <- matrix(data = c(mean(stimMat[(t+1-avg_window):t+1 , 1]),
                                    sd(stimMat[(t+1-avg_window):t+1 , 1]),
                                    mean(stimMat[(t+1-avg_window):t+1 , 2]),
                                    sd(stimMat[(t+1-avg_window):t+1 , 2])),
                           nrow = 1,
                           dimnames = list(NULL,
                                           c("Stim1Avg", "Stim1SD",
                                             "Stim2Avg", "Stim2SD")))
        # Bind
        all_info <- matrix(data = c(n, sim, chunk, t),
                           nrow = 1,
                           dimnames = list(NULL,
                                           c("n", "sim", "chunk", "t")))
        all_info <- cbind(all_info, window_entropy)
        all_info <- cbind(all_info, thresh_info)
        all_info <- cbind(all_info, stim_info)
        if (!exists("sim_info")) {
          sim_info <- all_info
        } else {
          sim_info <- rbind(sim_info, all_info)
        }
        rm(all_info, window_entropy, thresh_info, stim_info)
      }
    }
    check_info[[sim]] <- sim_info
    rm(sim_info)
  }
  check_info <- do.call('rbind', check_info)
  return(check_info)
  rm(check_info)
  sys.sleep(1)
})

sfStop()

parallel_simulations <- do.call("rbind", parallel_simulations)
# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
save(parallel_simulations, file =  paste0(storage_path, "CheckSimLength_Sigma", ThreshSD[1], "-Epsilon", epsilon, "-Beta", beta, ".Rdata"))

