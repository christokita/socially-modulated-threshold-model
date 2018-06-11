################################################################################
#
# Model incorporating both thresholds and network dynamics
#
################################################################################

rm(list = ls())
source("scripts/__Util__MASTER.R")


####################
# Set global variables
####################
# Initial paramters: Free to change
# Base parameters
Ns             <- c(40) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 20000 #number of generations to run simulation 
corrStep       <- 200 #number of time steps for calculation of correlation 
reps           <- 1 #number of replications per simulation (for ensemble)

# Threshold Parameters
ThreshM        <- rep(50, m) #population threshold means 
ThreshSD       <- ThreshM * 0.05 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.8, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
p              <- 0.5 #baseline probablity of initiating an interaction per time step
epsilon        <- 0.1 #relative weighting of social interactions for adjusting thresholds
beta           <- 1.1 #probability of interacting with individual in same state relative to others



####################
# Run simulation multiple times
####################

# Loop through group sizes
for (i in 1:length(Ns)) {
  # Set group size
  n <- Ns[i]
  
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
    
    # Prep correlation tracking matrix
    thresh1time <- list()
    thresh2time <- list()
    thresh1time[[1]] <- threshMat[,1]
    thresh2time[[2]] <- threshMat[,2]
    
    ####################
    # Simulate
    ####################
    # Run simulation
    for (t in 1:gens) { 
      # Current timestep is actually t+1 in this formulation, because first row is timestep 0
      # Update stimuli
      for (j in 1:ncol(stimMat)) {
        # update stim
        stimMat[t + 1, j] <- update_stim(stimulus = stimMat[t, j],
                                         delta = deltas[j], 
                                         alpha = alpha, 
                                         Ni = sum(X_g[ , j]), 
                                         n = n)
      }
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
      # threshMat <- adjust_thresh_social(social_network = g_adj,
      #                                   threshold_matrix = threshMat,
      #                                   state_matrix = X_g,
      #                                   epsilon = epsilon)
      threshMat <- adjust_thresholds_social_capped(social_network = g_adj,
                                                   threshold_matrix = threshMat,
                                                   state_matrix = X_g,
                                                   epsilon = epsilon,
                                                   threshold_max = 2 * ThreshM[1])
      thresh1time[[t + 1]] <- threshMat[,1]
      thresh2time[[t + 1]] <- threshMat[,2]
      
      # Update total task performance profile
      X_tot <- X_tot + X_g
      
    }

    # Print simulation completed
    print(paste0("DONE: N = ", n, ", Simulation ", sim))
  }
  
}

library(RColorBrewer)
library(scales)
library(tidyr)
library(ggthemes)

thresh1time <- do.call("rbind", thresh1time)
row.names(thresh1time) <- NULL
thresh1time <- as.data.frame(thresh1time)
thresh1time <- thresh1time %>% 
  mutate(t = 0:(nrow(.)-1)) %>% 
  gather("Id", "Threshold", 1:40)

threshMat <- threshMat %>% 
  as.data.frame(.) %>% 
  mutate(ThreshRatio = log(Thresh1 / Thresh2),
         Id = row.names(.))
threshMat$ThreshRatio[threshMat$ThreshRatio > 10] <- 10
threshMat$ThreshRatio[threshMat$ThreshRatio < -10] <- -10

thresh1time <- merge(thresh1time, threshMat, by = "Id")


gg_thresh <- ggplot(data = thresh1time, 
                    aes(x = t, y = Threshold)) +
  theme_classic(base_size = 10) +
  geom_line(aes(group = Id, colour = ThreshRatio), size = 0.2) +
  scale_colour_gradient2(name = "ln(Threshold Ratio)",
                         high = "#d7191c",
                         mid = "#ffffbf",
                         # mid = "#cecece",
                         low = "#2c7bb6", 
                         midpoint = 0, 
                         limits = c(-5, 5),
                         oob = squish) +
  scale_y_continuous(expand = c(0.0, 0), limits = c(0, 100)) +
  theme(aspect.ratio = 1,
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(4, "pt"),
        axis.text = element_text(color = "black"))
gg_thresh


entropy <- mutualEntropy(TotalStateMat = X_tot)

# ggsave("output/ThresholdTime/Size100_Sigma0.0_TripleTimeLength.png", scale = 0.6, dpi = 600)
