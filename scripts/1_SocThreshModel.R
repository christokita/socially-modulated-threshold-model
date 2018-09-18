
rm(list = ls())

####################
# Source necessary scripts/libraries
####################
source("scripts/util/__Util__MASTER.R")

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(5, 10, 20, 30, 50, 70) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 50000 #number of generations to run simulation 
reps           <- 20 #number of replications per simulation (for ensemble)

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
beta           <- 1.1 #probability of interacting with individual in same state relative to others


####################
# Run ensemble simulation
####################
# Prep meta-lists for collection of group size simulations
groups_taskDist    <- list()
groups_stim        <- list()
groups_thresh      <- list()
groups_entropy     <- list()
groups_graphs      <- list()

# Loop through group sizes
for (i in 1:length(Ns)) {
  # Set group size
  n <- Ns[i]
  # Prep lists for collection of simulation outputs from this group size
  ens_taskDist    <- list()
  ens_entropy     <- list()
  ens_stim        <- list()
  ens_thresh      <- list()
  ens_graphs      <- list()
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
    }
    
    ####################
    # Post run calculations
    ####################
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    # Calculate total task distribution
    totalTaskDist <- X_tot / gens
    # Create tasktally table
    stimMat <- cbind(stimMat, 0:(nrow(stimMat) - 1))
    colnames(stimMat)[ncol(stimMat)] <- "t"
    # Add total task distributions, entropy values, and graphs to lists
    ens_taskDist[[sim]]    <- totalTaskDist
    ens_entropy[[sim]]     <- entropy
    ens_stim[[sim]]        <- stimMat
    ens_thresh[[sim]]      <- threshMat 
    ens_graphs[[sim]]      <- g_tot / gens
  }
  # Add to list of lists
  groups_taskDist[[i]]    <- ens_taskDist
  groups_stim[[i]]        <- ens_stim
  groups_thresh[[i]]      <- ens_thresh
  groups_entropy[[i]]     <- ens_entropy
  groups_graphs[[i]]      <- ens_graphs
}


