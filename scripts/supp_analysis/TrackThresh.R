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

####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(80) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 50000 #number of generations to run simulation 
reps           <- 1 #number of replications per simulation (for ensemble)

# Threshold Parameters
ThreshM        <- rep(50, m) #population threshold means 
ThreshSD       <- ThreshM * 0 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.8, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active
thresh_max     <- 100

# Social Network Parameters
p              <- 1 #baseline probablity of initiating an interaction per time step
epsilon        <- 0.1 #relative weighting of social interactions for adjusting thresholds
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
                                                   threshold_max = thresh_max)
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
    totalTaskDist <- X_tot / gens
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
    ens_graphs[[sim]]      <- g_tot / gens
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

thresh_time1 <- do.call('rbind', thresh1time)
thresh_time1 <- as.data.frame(thresh_time1)
thresh_time1 <- thresh_time1 %>% 
  mutate(t = 0:(nrow(.)-1)) %>% 
  gather(., Id, Threshold, -t)

thresh_time2 <- do.call('rbind', thresh2time)
thresh_time2 <- as.data.frame(thresh_time2)
thresh_time2 <- thresh_time2 %>% 
  mutate(t = 0:(nrow(.)-1)) %>% 
  gather(., Id, Threshold, -t)

#Save
save(thresh_time1, thresh_time2, X_tot, file = paste0("output/ThresholdTime/Examples/n", n, "-eps", epsilon, "-beta", beta, ".Rdata"))

####################
# Analuyze data
####################
# Load
load('output/ThresholdTime/Examples/n80-eps0.1-beta1.1.Rdata')

# Plot time series
gg_threshtime <- ggplot(thresh_time1, aes(x = t, y = Threshold, group = Id)) +
  geom_line(size = 0.1, alpha = 0.1, colour = "#1f78b4") +
  scale_x_continuous(name = expression(paste("Time step (", italic(t), ")")),
                     breaks = c(1, seq(10000, 50000, 10000)),
                     labels = c("0", "", "", "", "", "50,000"),
                     expand = c(0, 0)) +
  scale_y_continuous(name = expression(paste("Task 1 threshold (", italic(theta[i1,t]), ")")),
                     limits = c(0, 100),
                     breaks = seq(0, 100, 50)) +
  theme_ctokita() +
  theme(axis.title = element_blank(),
    axis.text.x = element_text(hjust = 0.7))
gg_threshtime

ggsave(gg_threshtime, filename = paste0("output/ThresholdTime/Examples/TimeSeries_n", n, "-beta", beta, "-epsilon", epsilon, ".png"), 
                                        width = 31, height = 19.5, units = "mm", dpi = 400)

# Plot histogram of task performance
task_data <- as.data.frame(X_tot) / gens
task_data$bias <- task_data$Task2 - task_data$Task1
gg_taskdist <- ggplot(data = task_data, aes(x = Task1)) +
  geom_histogram(bins = 25, size = 0.2, color = "white", fill = "#1f78b4") +
  theme_ctokita() +
  scale_y_continuous(breaks = seq(0, 40, 20)) +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(hjust = 0.7))
gg_taskdist

ggsave(gg_taskdist, filename = paste0("output/ThresholdTime/Examples/TaskDist_n", n, "-beta", beta, "-epsilon", epsilon, ".svg"), 
       width = 30.5, height = 19.5, units = "mm", dpi = 400)
