################################################################################
#
# Social interaction model: non-parallelized simulation to also track thresholds
#
################################################################################

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
Ns             <- c(60) #vector of number of individuals to simulate
m              <- 2 #number of tasks
Tsteps         <- 50000 #number of time steps to run simulation 
reps           <- 1 #number of replications per simulation (for ensemble)

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
# Run ensemble simulation
####################
# Prep meta-lists for collection of group size simulations
groups_taskDist    <- list()
groups_taskTally   <- list()
groups_stim        <- list()
groups_thresh      <- list()
groups_entropy     <- list()
groups_thresh1Time <- list()
groups_thresh2Time <- list()
groups_graphs      <- list()

# Loop through group sizes
for (i in 1:length(Ns)) {
  # Set group size
  n <- Ns[i]
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
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    # Calculate total task distribution
    totalTaskDist <- X_tot / Tsteps
    # Create tasktally table
    stimMat <- cbind(stimMat, 0:(nrow(stimMat) - 1))
    colnames(stimMat)[ncol(stimMat)] <- "t"
    # Add total task distributions, entropy values, and graphs to lists
    ens_taskDist[[sim]]    <- totalTaskDist
    ens_entropy[[sim]]     <- entropy
    ens_taskTally[[sim]]   <- taskTally
    ens_stim[[sim]]        <- stimMat
    ens_thresh[[sim]]      <- threshMat 
    ens_thresh1Time[[sim]] <- thresh1time
    ens_thresh2Time[[sim]] <- thresh2time
    ens_graphs[[sim]]      <- g_tot / Tsteps
  }
  # Add to list of lists
  groups_taskDist[[i]]    <- ens_taskDist
  groups_taskTally[[i]]   <- ens_taskTally
  groups_stim[[i]]        <- ens_stim
  groups_thresh[[i]]      <- ens_thresh
  groups_entropy[[i]]     <- ens_entropy
  groups_thresh1Time[[i]] <- ens_thresh1Time
  groups_thresh2Time[[i]] <- ens_thresh2Time
  groups_graphs[[i]]      <- ens_graphs
}

thresh_time <- do.call('rbind', thresh1time)
thresh_time <- as.data.frame(thresh_time)
thresh_time <- thresh_time %>% 
  mutate(t = 1:nrow(.)) %>% 
  gather(., Id, Threshold, -t)

ggplot(thresh_time, aes(x = t, y = Threshold, group = Id)) +
  geom_line(size = 0.1, alpha = 0.1) +
  scale_x_continuous(label = comma) +
  theme_ctokita()
