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
Ns             <- c(50) #vector of number of individuals to simulate
m              <- 3 #number of tasks
gens           <- 4000 #number of generations to run simulation 
corrStep       <- 200 #number of time steps for calculation of correlation 
reps           <- 1 #number of replications per simulation (for ensemble)

# Threshold Parameters
ThreshM        <- rep(10, m) #population threshold means 
ThreshSD       <- ThreshM * 0.1 #population threshold standard deviations
InitialStim    <- rep(0, m) #intital vector of stimuli
deltas         <- rep(0.6, m) #vector of stimuli increase rates  
alpha          <- m #efficiency of task performance
quitP          <- 0.2 #probability of quitting task once active

# Social Network Parameters
epsilon        <- 0.005 #relative weighting of social interactions for lowering thresholds
phi            <- 0.005 #default forgetting rate of thresholds
p              <- 0.2 #probability of interacting with individual in other states
q              <- 1 #probability of interacting with individual in same state relative to others



####################
# Run simulation multiple times
####################
# Prep meta-lists for collection of group size simulations
groups_taskDist  <- list()
groups_taskCorr  <- list()
groups_taskStep  <- list()
groups_taskTally <- list()
groups_stim      <- list()
groups_entropy   <- list()
groups_graphs    <- list()
groups_specialization <- data.frame(NULL)

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
    taskOverTime  <- matrix(nrow = 0, ncol = n)
    
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
      # Calculate information on neighbors performing tasks
      L_g <- calcSocialInfo(SocialNetwork = g_adj,
                            X_sub_g = X_g)
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
      # Note which task is being peformed
      taskPerf <- matrix(nrow = 1, ncol = n)
      for (i in 1:nrow(X_g)) {
        task <- unname(which(X_g[i, ] == 1))
        if (length(task) == 0) {
          task <- 0
        }
        taskPerf[i] <- task
      }
      colnames(taskPerf) <- row.names(X_g)
      taskOverTime <- rbind(taskOverTime, taskPerf)
      
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
    
    # Calculate specialization of task performance 
    # from Gautrais et al. (2002)
    for (col in 1:ncol(taskOverTime)) {
      # Grab column of individual
      t_prof <- taskOverTime[ , col ]
      # Remove inactivity
      t_prof <- paste(t_prof, collapse = "")
      # Calculate transitions
      t_prof <- gsub("1+", "1", t_prof)
      t_prof <- gsub("2+", "2", t_prof)
      t_prof <- gsub("0+", "", t_prof)
      t_prof <- as.numeric(unlist(strsplit(as.character(t_prof), "")))
      transitions <- lapply(2:length(t_prof), function(entry) {
        a <- t_prof[entry] != t_prof[entry - 1]
      })
      C_i <- sum(unlist(transitions))
      C_i <- C_i / (length(t_prof) - 1)
      # Calulate specialization
      F_i <- 1 - m * C_i
      to_return <- data.frame(individual = paste0("v-", col), 
                              n = n,
                              replicate = sim,
                              TransSpec = F_i)
      groups_specialization <- rbind(groups_specialization, to_return)
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
    colnames(runMean) <- c("Task1", "Task2", "replicate")
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
  groups_entropy[[i]]   <- ens_entropy
  groups_graphs[[i]]    <- ens_graphs
  
}

# trim out correlations for group size 1
if(1 %in% Ns) {
  groups_taskCorr <- groups_taskCorr[-1]
}

# filename <- "Sigma01-Eps2-ConnectP03"
# 
# save(groups_entropy, groups_stim, groups_taskCorr, groups_taskDist,
#      groups_taskStep, groups_taskTally, groups_specialization,
#      file = paste0("output/", filename, ".Rdata"))
