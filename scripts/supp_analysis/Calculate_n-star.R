################################################################################
#
# Social interaction model: calculate n* via simulation (to compare to analytical results)
#
################################################################################
# We will compare DOL at 1,000 time steps to DOL at 6,000 time steps. Where the change in DOL is 0, we can assume that is n*
# as this is where we expect probability of interacting with each behavioral type to be equal

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
m              <- 2 #number of tasks
gens           <- 6000 #number of generations to run simulation 
reps           <- 100 #number of replications per simulation (for ensemble)
chunk_size     <- 10 #number of simulations sent to single core 

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
betas          <- seq(1.025, 1.25, 0.025) #probability of interacting with individual in same state relative to others

# Where to check time values
times <- c(1000, 6000)

####################
# Prep for Parallelization
####################

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
parallel_simulations <- sfLapply(1:length(betas), function(k) {
  # Set beta
  beta <- betas[k]
  # Set group sizes 
  n_star_calculated <- round((m * beta) / (deltas[1] * (beta - m + 1)))
  Ns <- seq(n_star_calculated - 5, n_star_calculated + 5, 1)
  # Prep lists for collection of simulation outputs from this group size
  group_size_data <- lapply(1:length(Ns), function(x) {
    # set group size 
    n <- Ns[x]
    # set list for collection of sim_data 
    group_size_list <- list()
    # Run Simulations
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
          # Bind
          all_info <- matrix(data = c(n, beta, sim, t),
                             nrow = 1,
                             dimnames = list(NULL,
                                             c("n", "beta", "sim", "t")))
          all_info <- cbind(all_info, window_entropy)
          if (!exists("sim_info")) {
            sim_info <- all_info
          } else {
            sim_info <- rbind(sim_info, all_info)
          }
          rm(all_info, window_entropy)
        }
      }
      group_size_list[[sim]] <- sim_info
      print(paste0("Done ", sim))
      rm(sim_info)
    } #end of replicate simulations loop
    to_return <- do.call("rbind", group_size_list)
    return(to_return)
  }) # end of group size loop
  sys.sleep(1)
})

sfStop()

parallel_simulations <- do.call("rbind", parallel_simulations)
# Create directory for depositing data
storage_path <- "/scratch/gpfs/ctokita/"
save(parallel_simulations, file =  paste0(storage_path, "CalculateNstar-", ThreshM[1], "_Sigma", ThreshSD[1], "-Epsilon", epsilon, ".Rdata"))

