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
Ns             <- c(20) #vector of number of individuals to simulate
m              <- 2 #number of tasks
gens           <- 10000 #number of generations to run simulation 
corrStep       <- 200 #number of time steps for calculation of correlation 
reps           <- 1 #number of replications per simulation (for ensemble)

# Threshold Parameters
ThreshM        <- rep(10, m) #population threshold means 
ThreshSD       <- ThreshM * 0.01 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.6, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
epsilon        <- 0.01 #relative weighting of social interactions for lowering thresholds #0.01 = epsilon = phi
phi            <- 0.01 #default forgetting rate of thresholds
p              <- 0.1 #probability of interacting with individual in other states
q              <- 1.1 #probability of interacting with individual in same state relative to others



####################
# Run simulation multiple times
####################
# Prep meta-lists for collection of group size simulations
groups_taskDist  <- list()
groups_taskCorr  <- list()
groups_taskStep  <- list()
groups_taskTally <- list()
groups_stim      <- list()
groups_thresh    <- list()
groups_entropy   <- list()
groups_graphs    <- list()

# Loop through group sizes
for (i in 1:length(Ns)) {
  # Set group size
  n <- Ns[i]
  
  # Prep lists for collection of simulation outputs
  ens_taskDist  <- list()
  ens_taskCorr  <- list()
  ens_taskStep  <- list()
  ens_taskTally <- list()
  ens_entropy   <- list()
  ens_stim      <- list()
  ens_thresh    <- list()
  ens_graphs    <- list()
  
  # Run Simulations
  for (sim in 1:reps) {
    
    ####################
    # Seed structures and intial matrices
    ####################

    # Set initial probability matrix (P_g)
    P_g <- matrix(data = rep(0, n * m), ncol = m)
    
    # Seed task (external) stimuli
    stimMat <- seedStimuls(InitialSVector = InitialStim, 
                           gens = gens)
    
    # Seed internal thresholds
    threshMat <- seedThresholds(n = n, 
                                m = m, 
                                ThresholdMeans = ThreshM, 
                                ThresholdSDs = ThreshSD)
    
    # Start task performance
    X_g <- matrix(data = rep(0, length(P_g)), ncol = ncol(P_g))
    
    # Create cumulative task performance matrix
    X_tot <- X_g
    
    # Create cumulative adjacency matrix
    g_tot <-  matrix(data = rep(0, n * n), ncol = n)
    colnames(g_tot) <- paste0("v-", 1:n)
    rownames(g_tot) <- paste0("v-", 1:n)
    
    # Prep correlation step matrix
    X_prev <- matrix(data = rep(0, n * m), ncol = m)
    X_prevTot <- matrix(data = rep(0, n * m), ncol = m)
    taskCorr <- list()
    taskStep <- list()
    taskTally <- list()
    
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
      # Update stimuli
      for (j in 1:ncol(stimMat)) {
        # update stim
        stimMat[t + 1, j] <- globalStimUpdate(stimulus = stimMat[t, j],
                                              delta = deltas[j], 
                                              alpha = alpha, 
                                              Ni = sum(X_g[ , j]), 
                                              n = n)
      }
      # Update social network
      g_adj <- temporalNetwork(X_sub_g = X_g,
                               p = p,
                               bias = q)
      g_tot <- g_tot + g_adj
      # Calculate task demand based on global stimuli
      P_g <- calcThresholdDetermMat(TimeStep = t + 1, # first row is generation 0
                                    ThresholdMatrix = threshMat, 
                                    StimulusMatrix = stimMat)
      # Update task performance
      X_g <- updateTaskPerformance(P_sub_g    = P_g,
                                   TaskMat    = X_g,
                                   QuitProb   = quitP)
      # Adjust thresholds
      threshMat <- adjustThresholdsSocial(SocialNetwork = g_adj,
                                          ThresholdMatrix = threshMat, 
                                          X_sub_g = X_g, 
                                          epsilon = epsilon, 
                                          phi = phi)
      thresh1time[[t + 1]] <- threshMat[,1]
      thresh2time[[t + 1]] <- threshMat[,2]
      
      
      # Capture current task performance tally
      tally <- matrix(c(t, colSums(X_g)), ncol = ncol(X_g) + 1)
      colnames(tally) <- c("t", colnames(X_g))
      tally <- transform(tally, Inactive = n - sum(X_g), n = n, replicate = sim)
      taskTally[[t]] <- tally
      
      # Update total task performance profile
      X_tot <- X_tot + X_g
      
      # Create time step for correlation
      if (t %% corrStep == 0) {
        # Get tasks performance in correlation step
        X_step <- X_tot - X_prevTot
        # Add to ensemble list of task steps
        taskStep[[t / corrStep]] <- X_step
        # Calculate rank correlation if it is not the first step
        if(sum(X_prev) != 0) {
          # Normalize
          stepNorm <- X_step / rowSums(X_step)
          prevNorm <- X_prev / rowSums(X_prev)
          # Calculate ranks
          step_ranks <- calculateTaskRank(TaskStepMat = X_step)
          prev_ranks <- calculateTaskRank(TaskStepMat = X_prev)
          # Calculate Correlation
          rankCorr <- cor(prev_ranks, step_ranks, method = "spearman")
          # Put in list
          taskCorr[[(t / corrStep) - 1]] <- diag(rankCorr)
          names(taskCorr)[(t / corrStep) - 1] <- paste0("Gen", t)
        }
        # Update previous step total matrix
        X_prevTot <- X_tot
        # Update previous step total matrix
        X_prev <- X_step
      }
    }
    
    # Calculate Entropy
    entropy <- mutualEntropy(TotalStateMat = X_tot)
    entropy <- transform(entropy, n = n, replicate = sim)
    
    # Calculate total task distribution
    # totalTaskDist <- X_tot / rowSums(X_tot)
    totalTaskDist <- X_tot / gens
    totalTaskDist <- transform(totalTaskDist, Inactive = gens - rowSums(X_tot), n = n, replicate = sim)
    
    # Create tasktally table
    taskTally <- do.call("rbind", taskTally)
    
    # Create tasktally table
    stimMat <- transform(stimMat, n = n, replicate = sim)
    
    # Create tasktally table
    taskCorr <- transform(taskCorr, replicate = sim)
    
    # Add total task distributions, entropy values, and graphs to lists
    ens_taskDist[[sim]]  <- totalTaskDist
    ens_entropy[[sim]]   <- entropy
    ens_taskCorr[[sim]]  <- taskCorr
    ens_taskTally[[sim]] <- taskTally
    ens_taskStep[[sim]]  <- taskStep
    ens_stim[[sim]]      <- stimMat
    ens_thresh[[sim]]    <- threshMat 
    ens_graphs[[sim]]    <- g_tot / gens
    
    # Print simulation completed
    print(paste0("DONE: N = ", n, ", Simulation ", sim))
  }
  
  # Calculate mean correlation for each n
  runCorrs <- lapply(ens_taskCorr, function(x) {
    # Unlist
    runs <- do.call("rbind", x)
    replicate <- runs[nrow(runs), ]
    replicate <- unique(replicate)
    runs <- runs[-nrow(runs), ]
    # Calculate mean
    runMean <- matrix(data = rep(NA, m), ncol =  m)
    for (column in 1:m) {
      runMean[ , column] <- mean(runs[ , column], na.rm = TRUE)
    }
    runMean <- cbind(runMean, replicate)
    colnames(runMean) <- c(paste0("Task", 1:m), "replicate")
    return(runMean)
  })
  runCorrs <- do.call("rbind", runCorrs)
  runCorrs <- transform(runCorrs, n = n)
  
  # Add to list of lists
  groups_taskDist[[i]]  <- ens_taskDist
  groups_taskCorr[[i]]  <- runCorrs
  groups_taskStep[[i]]  <- ens_taskStep
  groups_taskTally[[i]] <- ens_taskTally
  groups_stim[[i]]      <- ens_stim
  groups_thresh[[i]]    <- ens_thresh
  groups_entropy[[i]]   <- ens_entropy
  groups_graphs[[i]]    <- ens_graphs
  
}

# trim out correlations for group size 1
if(1 %in% Ns) {
  groups_taskCorr <- groups_taskCorr[-1]
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
  gather("Id", "Threshold", 1:70)

threshMat <- threshMat %>% 
  as.data.frame(.) %>% 
  mutate(ThreshRatio = log(Thresh1 / Thresh2),
         Id = row.names(.)) %>% 
  select(-Thresh1, -Thresh2)

thresh1time <- merge(thresh1time, threshMat, by = "Id")


base_breaks_x <- function(x){
  b <- pretty(x)
  adjust_for_ticks <- max(b) * 0.0025
  d <- data.frame(y=-Inf, yend=-Inf, x=min(b) - adjust_for_ticks, xend=max(b) + adjust_for_ticks)
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend),size = 1, inherit.aes=FALSE),
       scale_x_continuous(breaks=b))
}
base_breaks_y <- function(x){
  b <- pretty(x)
  adjust_for_ticks <- max(b) * 0.0025
  d <- data.frame(x = -Inf, xend = -Inf, y=min(b) - adjust_for_ticks, yend = max(b) + adjust_for_ticks)
  list(geom_segment(data=d, aes(x = x, y = y, xend = xend, yend = yend), size = 1, inherit.aes = FALSE),
       scale_y_continuous(breaks=b))
}

test <- thresh1time[thresh1time$Id %in% c("v-22", "v-20", "v-25", "v-40"),]
gg_thresh <- ggplot(data = thresh1time, 
                    aes(x = t, y = Threshold)) +
  theme_bw(base_size = 10) +
  geom_line(aes(group = Id, colour = ThreshRatio), size = 0.5) +
  scale_colour_gradient2(name = "ln(Threshold Ratio)",
                         high = "#d7191c",
                         mid = "#f2ec79", 
                         low = "#2c7bb6", 
                         midpoint = 0, 
                         limits = c(-2, 2),
                         oob = squish) +
  base_breaks_x(thresh1time$t) +
  base_breaks_y(thresh1time$Threshold) +
  theme(aspect.ratio = 1,
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks.length = unit(4, "pts")
        axis.text = element_text(color = "black"))
gg_thresh

ggsave("output/ThresholdTime/Size70.png", scale = 0.7)
